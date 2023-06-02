# A latent-competing risks model of the
# feto-infant mortality trajectory over
# age of gestation

# Init ------------------------------------------------------------

library(tidyverse)
library(maxLik)
library(cowplot)

paths <- list()
paths$input <- list(
  fetoinfant_lifetables = 'out/30-fetoinfant_lifetables.RData',
  figure_specs = 'src/00-figure_specifications.R',
  lifetable_functions = 'src/00-fnct-feto_infant_lt.R'
)
paths$output <- list(
  fetoinfant = 'tmp/21-fetoinfant.RData'
)

# figure specs
source(paths$input$figure_specs)
# fetoinfant lifetable functions
source(paths$input$lifetable_functions)

# constants
cnst <-
  list(
    gestage_brk = seq(24, 77, by = 4),
    lifetable_breaks = 24:77,
    left_truncation_gestage = 24,
    right_censoring_gestage = 77
  )

# figures
fig <- list()
# tables
tab <- list()
# model fits
fit <- list()

# Functions -------------------------------------------------------

# [a,b] -> [-Inf, Inf]
ScaleLogit <- function(x, a, b) {
  x_norm <- (x - a) / (b - a)
  logit_x <- log(x_norm / (1 - x_norm))
  
  return(logit_x)
}

# [-Inf,Inf] -> [a,b]
ScaleInverseLogit <- function(logit_x, a, b) {
  x_norm <- 1 / (1 + exp(-logit_x))
  # scale x back to [a, b]
  x <- a + x_norm * (b - a)
  
  return(x)
}

# Data ------------------------------------------------------------

# a list of feto-infant lifetable as FILT objects
load(paths$input$fetoinfant_lifetables)

# Feto-infant parametric survival -------------------------------

# baseline characteristics
NegGompertzHzrd <-
  function (x, a, b, alpha, beta, split = 21) {
    I1 = ifelse(x < split, 1, 0)
    I2 = ifelse(x >= split, 1, 0)
    a*exp(-b*x)*I1 + alpha*exp(-beta*(x-split))*I2
  }
NegGompertzCumHzrd <-
  function (x, a, b, alpha, beta, split = 21) {
    I1 = ifelse(x < split, 1, 0)
    I2 = ifelse(x >= split, 1, 0)
    Ssplit = ((a - a*exp(-b*split)) / b)
    I1*((a - a*exp(-b*x)) / b) +
      I2*((alpha - alpha*exp(-beta*(x-split))) / beta + Ssplit)
  }
NegGompertzSurv <-
  function (x, a, b, alpha, beta, split = 21) {
    exp(-NegGompertzCumHzrd(x, a, b, alpha, beta, split))
  }

# birth hump characteristics
GaussianHzrd <-
  function (x, a, c, s) {
    a*exp(-(x-c)^2/(2*s^2))
  }
GaussianCumHzrd <-
  function (x, a, c, s) {
    denom <- sqrt(2)*s
    Hx <- a*sqrt(pi/2)*s*(pracma::erf(c/denom) + pracma::erf((x-c)/denom))
    return(Hx)
  }
GaussianSurv <-
  function (x, a, c, s) {
    exp(-GaussianCumHzrd(x, a, c, s))
  }

# feto-infant competing risks survival
FetoinfantSurv <-
  function (x, pars, split = 15, component = 'total') {
    a1 = exp(pars[1]); b = exp(pars[2])
    alpha = exp(pars[3]); beta = exp(pars[4]);
    a2 = exp(pars[5]); c = ScaleInverseLogit(pars[6], a = 10, b = 20); s = exp(pars[7])

    # ontogenescent survival
    ontogen_surv <-
      NegGompertzSurv(
        x = x,
        a = a1, b = b,
        alpha = alpha, beta = beta,
        split = split
      )
    # birth survival
    birth_surv <-
      GaussianSurv(
        x = x,
        a = a2, c = c, s = s
      )

    Sx <-
      switch (component,
              'total' = ontogen_surv * birth_surv,
              'ontogen' = ontogen_surv,
              'birth' = birth_surv,
              stop('Component must be one of "total", "ontogen" or "birth"')
      )
    return(Sx)
  }

# feto-infant competing risks hazard
FetoinfantHzrd <-
  function (x, pars, split = 15, component = 'total') {
    a1 = exp(pars[1]); b = exp(pars[2])
    alpha = exp(pars[3]); beta = exp(pars[4]);
    a2 = exp(pars[5]); c = ScaleInverseLogit(pars[6], a = 10, b = 20); s = exp(pars[7])
    
    # ontogenescent hazard
    ontogen_hzrd <-
      NegGompertzHzrd(
        x = x,
        a = a1, b = b,
        alpha = alpha, beta = beta, split = split
      )

    # birth hazard
    birth_hzrd <-
      GaussianHzrd(
        x = x,
        a = a2, c = c, s = s
      )

    hx <-
      switch (component,
              'total' = ontogen_hzrd + birth_hzrd,
              'ontogen' = ontogen_hzrd,
              'birth' = birth_hzrd,
              stop('Component must be one of "total", "ontogen" or "birth"')
      )
    return(hx)

  }

# Objective function --------------------------------------------

# interval censored likelihood
IntervalCensoredLogLike <-
  function (pars, age, width, obsDx, obsCx, SurvFnct, lambda = 0, split, ...) {

    # predict log-hazard on basis of parameter estimates
    predSurvL <- SurvFnct(x = age, pars, split, ...)
    predSurvR <- SurvFnct(x = age+width, pars, split, ...)

    loglike <- obsDx*log(predSurvL-predSurvR) + obsCx*log(predSurvR) -
      # penalize discontinuities between the two ontogenescent segments
      # penalize birth hump magnitude
      lambda*(sum(abs(pars[5])) + abs(pars[3] - (pars[1]-exp(pars[2])*21)))
    return(loglike)

  }

# Fit fetoinfant survival ---------------------------------------

# Fit a parametric model of fetoinfant survival to
# the fetoinfant lifetable assuming left truncation age
# at 24 and right censoring at 77. Derive various predictions.
# Sample parameters from multivariate normal distribution
# derived from hessian in order to calculate CIs.
FitFetoinfantSurvival <-
  function (filt, lambda = 0, split = 15, simulate = TRUE) {

    stopifnot(any(class(filt) == 'FILT'))

    fit_lifetable <-
      filt %>%
      left_join(FILTMortalityRates(.)) %>%
      left_join(FILTSurvival(.)) %>%
      group_by(stratum) %>%
      group_modify(
        ~ {

          # initialize model parameters via heuristics
          init_pars <-
            .x %>%
            summarise(
              # log hazard at t=0
              loga1 = log(ifelse(first(m) == 0, mean(m), first(m))),
              # log relative rate of ontogenescent mortality decline
              logb =
                log(abs((log(median(m[1:5]))-
                           log(median(m[40:44])))/(x[1]-x[40]))),
              logb = ifelse(is.infinite(logb), -5, logb),
              # log hazard at t=45-24=21
              logalpha = loga1-exp(logb)*21,
              # log relative rate of decline after week 39
              logbeta =
                logb,
              # log hazard contribution of birth component at mode
              loga2 = log(exp(loga1)/2),
              # location of birth component (0 gets rescaled to 39 weeks)
              logc = 0,
              # spread of birth component
              logs = log(1)
            ) %>%
            unlist()

          model <-
            maxLik(
              logLik = IntervalCensoredLogLike,
              start  = init_pars,
              method = 'CG',
              # data and arguments to objective function
              age =
                pull(.x, x)-24,
              width =
                pull(.x, n),
              obsDx =
                pull(.x, D),
              obsCx =
                pull(.x, C),
              # hazard function
              SurvFnct =
                FetoinfantSurv,
              lambda = lambda,
              split = split,
              # options
              iterlim = 1e4
            )

          # 1000 draws from the posterior parameter distribution
          # assuming multivariate normal derived from hessian
          par_draw <-
            expand_grid(
            draw = 1:1000,
            name = c('1:a1', '2:b', '3:alpha', '4:beta', '5:a2', '6:c', '7:s')
            )
          if (isTRUE(simulate)) {
            par_draw$value <- MASS::mvrnorm(
              n = 1e3,
              mu = model$estimate,
              Sigma = chol2inv(-model$hessian)
            ) %>% t() %>% c()
          } else {
            par_draw$value <- rep(model$estimate, 1e3)
          }
          
          # summarise parameter draws into
          # point estimates and credible intervals
          pars <-
            par_draw %>%
            group_by(name) %>%
            summarise(
              avg = mean(value),
              se = sd(value),
              ci025 = quantile(value, 0.025),
              ci975 = quantile(value, 0.975)
            )

          # number of weeks from observation
          # start to end
          omega <-
            cnst$right_censoring_gestage -
            cnst$left_truncation_gestage

          # predictions by posterior draw
          pred_draw <-
            par_draw %>%
            group_by(draw) %>%
            group_modify(~{
              tibble(
                # weeks since left truncation age
                x =
                  seq(0, omega, length.out = 1000),
                # width of age interval
                n = omega/1000,
                # probability of fetoinfant survival until x
                total_lx =
                  FetoinfantSurv(
                    x,
                    pars = .x$value,
                    component = 'total'
                  ),
                # probability of fetoinfant death until x
                total_Fx =
                  1-total_lx,
                # probability of fetoinfant death until x
                # (one in x won't survive)
                total_iFx =
                  1/total_Fx,
                # hazard of fetoinfant death at x
                total_hx =
                  FetoinfantHzrd(
                    x,
                    pars = .x$value,
                    component = 'total'
                  ),
                # hazard of birth component at x
                birth_hx =
                  FetoinfantHzrd(
                    x,
                    pars = .x$value,
                    component = 'birth'
                  ),
                # hazard of ontogenescent component at x
                ontogen_hx =
                  FetoinfantHzrd(
                    x,
                    pars = .x$value,
                    component = 'ontogen'
                  ),
                # cumulative probability of death due
                # to birth component
                birth_Fx =
                  cumsum(total_lx*birth_hx*n),
                # cumulative probability of death due
                # to ontogenescent component
                ontogen_Fx =
                  cumsum(total_lx*ontogen_hx*n),
                # share of deaths due to birth component
                p_birth =
                  birth_Fx/total_Fx
              ) %>%
                mutate(
                  # transform back to gestational age
                  x = x + cnst$left_truncation_gestage
                )
            }) %>%
            ungroup()
          
          # summarise predictions over draws
          pred <-
            pred_draw %>%
            pivot_longer(
              -c('draw', 'x', 'n')
            ) %>%
            group_by(name, x, n) %>%
            summarise(
              avg = mean(value, na.rm = TRUE),
              se = sd(value, na.rm = TRUE),
              q025 = quantile(value, 0.025, na.rm = TRUE),
              q975 = quantile(value, 0.975, na.rm = TRUE)
            ) %>%
            pivot_wider(
              names_from = 'name',
              values_from = c('avg', 'se', 'q025', 'q975')
            ) %>%
            ungroup()

          tibble(
            convergence = model$code,
            loglike = maxValue(model),
            par_draws = list(par_draw),
            par_summary = list(pars),
            hessian = list(model$hessian),
            model = list(model),
            pred_draws = list(pred_draw),
            pred_summary = list(pred),
            lifetable = list(.x)
          )

        }) %>%
      ungroup()

    fit_lifetable

  }

# by sex
fit$sex09 <-
  FitFetoinfantSurvival(filt$sex09)

# by cohort
fit$cohort <-
  FitFetoinfantSurvival(filt$cohort)

# by origin
fit$origin <-
  FitFetoinfantSurvival(filt$origin09)

# by cause of death
fit$maternal <- FitFetoinfantSurvival(filt$maternal, lambda = 0.01)
fit$TreatableNeoplasms <- FitFetoinfantSurvival(filt$TreatableNeoplasms, lambda = 0.01)
fit$UntreatableNeoplasms <- FitFetoinfantSurvival(filt$UntreatableNeoplasms, lambda = 0.01)
fit$InfectionsParacitesOperations <- FitFetoinfantSurvival(filt$InfectionsParacitesOperations, simulate = FALSE,
                                      lambda = 0.01)
fit$ViolenceAccidents <- FitFetoinfantSurvival(filt$ViolenceAccidents, simulate = FALSE,
                                      lambda = 0.1) 

fit$PCML <- FitFetoinfantSurvival(filt$PCML, lambda = 0.01)
fit$Prematurity <- FitFetoinfantSurvival(filt$Prematurity, lambda = 0.01)
fit$Other <- FitFetoinfantSurvival(filt$Other, lambda = 0.01)

fit$UnspecificStillbirth <- FitFetoinfantSurvival(filt$UnspecificStillbirth)
fit$SuddenInfantDeath <- FitFetoinfantSurvival(filt$SuddenInfantDeath, lambda = 0.1) 

# Plot hazards --------------------------------------------------

# This function requires the output of FitFetoinfantSurvival()
# and plots lifetable estimates of fetoinfant survival and mortality
# versus the parametric model estimates.
PlotHazards <-
  function(
    filt_fit,
    ylab,
    xbrk = cnst$gestage_brk,
    xlim =
      c(cnst$left_truncation_gestage,
      cnst$right_censoring_gestage-0.1),
    scaler = 1e5
  ) {

    # aspect ratio
    ar <- 0.85

    lifetables <-
      filt_fit %>% unnest_legacy(lifetable)
    predictions <-
      filt_fit %>% unnest_legacy(pred_summary)

    plot_hzrd <-
      lifetables %>%
      ggplot() +
      geom_point(
        aes(
          x = x+0.5*n, y = m*scaler, group = stratum,
          color = as.character(stratum)
        ),
        size = fig_spec$point_size_m,
        alpha = 0.2
      ) +
      geom_line(
        aes(
          x = x, y = avg_total_hx*scaler, group = stratum,
          color = as.character(stratum)
        ),
        size = fig_spec$line_size_m,
        data = predictions
      ) +
      geom_ribbon(
        aes(
          x = x,
          ymin = q025_total_hx*scaler,
          ymax = q975_total_hx*scaler,
          group = stratum,
          fill = as.character(stratum)
        ),
        alpha = 0.2,
        data = predictions
      ) +
      geom_vline(
        aes(xintercept = 40),
        lty = 3,
        size = fig_spec$line_size_m
      ) +
      scale_y_continuous(
        paste0('Feto-infant deaths per\n',
               formatC(scaler, format = 'd', big.mark = ','),
               ' person-weeks at risk'),
        breaks = c(seq(2, 10, 2), seq(20, 80, 20)),
        trans = 'log10'
      ) +
      scale_x_continuous(
        'Week of gestation',
        breaks = xbrk,
        limits = xlim
      ) +
      scale_color_manual(
        '',
        values = fig_spec$discrete_colors
      ) +
      scale_fill_manual(
        '',
        values = fig_spec$discrete_colors
      ) +
      fig_spec$MyGGplotTheme(ar = ar) +
      theme(
        legend.position = c(0.8, 0.8)
      )

    plot_surv <-
      lifetables %>%
      ggplot() +
      geom_vline(
        aes(xintercept = 40),
        lty = 3,
        size = fig_spec$line_size_m
      ) +
      geom_point(
        aes(
          x = x, y = F_empirical*scaler,
          color = as.character(stratum)
        ),
        alpha = 0.2,
        size = fig_spec$point_size_m
      ) +
      geom_ribbon(
        aes(
          x = x,
          ymin = (1-q025_total_lx)*scaler,
          ymax = (1-q975_total_lx)*scaler,
          fill = as.character(stratum),
        ),
        alpha = 0.2,
        data = predictions
      ) +
      geom_line(
        aes(
          x = x, y = (1-avg_total_lx)*scaler,
          color = as.character(stratum),
        ),
        size = fig_spec$line_size_m,
        data = predictions
      ) +
      scale_y_continuous(
        paste0('Cumulative feto-infant deaths\n',
               'out of ', formatC(scaler, format = 'd',
                                  big.mark = ','),
               ' cohort members')
      ) +
      scale_color_manual(values = fig_spec$discrete_colors) +
      scale_x_continuous(
        '', breaks = xbrk
      ) +
      scale_fill_manual(
        '',
        values = fig_spec$discrete_colors
      ) +
      fig_spec$MyGGplotTheme(ar = ar) +
      theme(
        legend.position = 'none'
      )

    plot_grid(plot_hzrd, plot_surv, ncol = 2, align = 'h')

  }

# by sex
fig$hzrd_sex09 <- PlotHazards(fit$sex09)
fig_spec$ExportPDF(
  fig$hzrd_sex09,
  '50-hzrd-sex09',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# by cohort
fig$hzrd_cohort <- PlotHazards(fit$cohort)
fig_spec$ExportPDF(
  fig$hzrd_cohort,
  '50-hzrd-cohort',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# by origin
fig$hzrd_origin <- PlotHazards(fit$origin)
fig_spec$ExportPDF(
  fig$hzrd_origin,
  '50-hzrd-origin',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# by cause of death

fig$hzrd_maternal <- PlotHazards(fit$maternal)
fig_spec$ExportPDF(
  fig$hzrd_maternal,
  '50-hzrd-maternal_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig$hzrd_TreatableNeoplasms <- PlotHazards(fit$TreatableNeoplasms)
fig_spec$ExportPDF(
  fig$hzrd_TreatableNeoplasms,
  '50-hzrd-TreatableNeoplasms_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig$hzrd_UntreatableNeoplasms <- PlotHazards(fit$UntreatableNeoplasms)
fig_spec$ExportPDF(
  fig$hzrd_UntreatableNeoplasms,
  '50-hzrd-UntreatableNeoplasms_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)


fig$hzrd_InfectionsParacitesOperations <- PlotHazards(fit$InfectionsParacitesOperations)
fig_spec$ExportPDF(
  fig$hzrd_InfectionsParacitesOperations,
  '50-hzrd-InfectionsParacitesOperations_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig$hzrd_ViolenceAccidents <- PlotHazards(fit$ViolenceAccidents) #notworking
fig_spec$ExportPDF(
  fig$hzrd_ViolenceAccidents,
  '50-hzrd-ViolenceAccidents_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig$hzrd_UnspecificStillbirth <- PlotHazards(fit$UnspecificStillbirth)
fig_spec$ExportPDF(
  fig$hzrd_UnspecificStillbirth,
  '50-hzrd-UnspecificStillbirth_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig$hzrd_SuddenInfantDeath <- PlotHazards(fit$SuddenInfantDeath) #notworking
fig_spec$ExportPDF(
  fig$hzrd_SuddenInfantDeath,
  '50-hzrd-SuddenInfantDeath_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig$hzrd_pcml <- PlotHazards(fit$PCML)
fig_spec$ExportPDF(
  fig$hzrd_pcml,
  '50-hzrd-pcml_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# fig$hzrd_respiratory <- PlotHazards(fit$Respiratory)
# fig_spec$ExportPDF(
#   fig$hzrd_respiratory,
#   '50-hzrd-respiratory',
#   'out',
#   width = fig_spec$width,
#   height = fig_spec$width*0.4
# )

fig$hzrd_other <- PlotHazards(fit$Other)
fig_spec$ExportPDF(
  fig$hzrd_other,
  '50-hzrd-other_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

fig$hzrd_prematurity <- PlotHazards(fit$Prematurity)
fig_spec$ExportPDF(
  fig$hzrd_prematurity,
  '50-hzrd-prematurity_v2',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.4
)

# Share of deaths due to birth hump -----------------------------

BirthHumpDeaths <- function (filt_fit, x) {

  filt_fit %>%
    unnest_legacy(pred_summary) %>%
    filter(x == {{x}}) %>%
    select(stratum, avg_p_birth,
           se_p_birth,
           q025_p_birth, q975_p_birth) %>%
    ungroup()
}

# by sex
BirthHumpDeaths(fit$sex09, x = 77)
# by cohort
BirthHumpDeaths(fit$cohort, x = 77)
# by origin
BirthHumpDeaths(fit$origin, x = 77)

# by cause of death
BirthHumpDeaths(fit$maternal, x = 77)
BirthHumpDeaths(fit$TreatableNeoplasms, x = 77)
BirthHumpDeaths(fit$UntreatableNeoplasms, x = 77)
BirthHumpDeaths(fit$ViolenceAccidents, x = 77)
BirthHumpDeaths(fit$UnspecificStillbirth, x = 77)
BirthHumpDeaths(fit$Other, x = 77)
BirthHumpDeaths(fit$SuddenInfantDeath, x = 77)
BirthHumpDeaths(fit$InfectionsParacitesOperations, x = 77)
BirthHumpDeaths(fit$PCML, x = 77)
BirthHumpDeaths(fit$Prematurity, x = 77)

# Plot transitional component -------------------------------------

cod <- c('maternal', 'TreatableNeoplasms','UntreatableNeoplasms', 'ViolenceAccidents', 'UnspecificStillbirth',
         'PCML', 'Prematurity', 'SuddenInfantDeath','InfectionsParacitesOperations', 'Other')

filt_cod <- map(cod, ~{
  fit[[.x]]$pred_summary[[1]] %>%
    mutate(cod = .x)
}) %>% bind_rows() %>%
  group_by(cod) %>%
  mutate(
    cod = as_factor(cod),
    max = max(avg_birth_hx, na.rm = TRUE)
  ) %>%
  ungroup()

fig$birthhump_cod_joint <-
  filt_cod %>%
  mutate(cod = fct_reorder(cod, .x = max, .desc = TRUE)) %>%
  ggplot(aes(x = x)) +
  geom_area(aes(y = avg_birth_hx*1e5, fill = cod)) +
  scale_fill_brewer(type = 'q', palette = 3) +
  #facet_wrap(~cod, ncol = 1) +
  fig_spec$MyGGplotTheme(panel_border = TRUE) +
  labs(y = 'Feto-infant deaths per 100K person-weeks',
       x = 'Week of gestation',
       fill = 'Cause of death')

fig_spec$ExportPDF(
  fig$birthhump_cod_joint,
  '50-birthhump_cod_joint',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.6
)

fig$birthhump_cod_separate <-
  filt_cod %>%
  mutate(cod = fct_reorder(cod, .x = max, .desc = TRUE)) %>%
  ggplot(aes(x = x)) +
  geom_area(aes(y = avg_birth_hx*1e5, fill = cod)) +
  scale_fill_brewer(type = 'q', palette = 3) +
  facet_wrap(~cod, nrow = 2) +
  fig_spec$MyGGplotTheme(panel_border = TRUE) +
  labs(y = 'Feto-infant deaths per 100K person-weeks',
       x = 'Week of gestation',
       fill = 'Cause of death')

fig_spec$ExportPDF(
  fig$birthhump_cod_separate,
  '50-birthhump_cod_separate',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.6
)

# Probability of Fetoinfant death ---------------------------------

ProbFetoInfantDeath <-
  function (filt_fit, x) {
    
    filt_fit %>%
      unnest_legacy(pred_summary) %>%
      filter(x == 77) %>%
      select(contains(c(
        'stratum',
        'total_Fx', 'total_iFx',
        'birth_Fx', 'birth_iFx'
      )))
    
  }

# by sex
ProbFetoInfantDeath(fit$sex09)
# by cohort
ProbFetoInfantDeath(fit$cohort)
# by origin
ProbFetoInfantDeath(fit$origin)

# by cohort
# by cause of death

bind_rows(
  maternal = ProbFetoInfantDeath(fit$maternal),
  TreatableNeoplasms = ProbFetoInfantDeath(fit$TreatableNeoplasms),
  UntreatableNeoplasms = ProbFetoInfantDeath(fit$UntreatableNeoplasms),
  ViolenceAccidents = ProbFetoInfantDeath(fit$ViolenceAccidents),
  UnspecificStillbirth = ProbFetoInfantDeath(fit$UnspecificStillbirth),
  Other = ProbFetoInfantDeath(fit$Other),
  SuddenInfantDeath = ProbFetoInfantDeath(fit$SuddenInfantDeath),
  PCML = ProbFetoInfantDeath(fit$PCML),
  Prematurity = ProbFetoInfantDeath(fit$Prematurity),
  InfectionsParacitesOperations = ProbFetoInfantDeath(fit$InfectionsParacitesOperations),
  .id = 'cod'
) %>%
  mutate(p_birth = avg_birth_Fx / sum(avg_birth_Fx)) %>%
  arrange(p_birth)

# Parameter tables --------------------------------------------------------

PrintParameterTable <- function (filt_fit) {

    tab_of_pars <-
      filt_fit %>%
      unnest_legacy(par_summary) %>%
      transmute(
        name,
        stratum,
        avg = exp(avg),
        ci025 = exp(ci025),
        ci975 = exp(ci975)
      ) %>%
      mutate_at(
        c('avg', 'ci025', 'ci975'),
        ~ formatC(., format = 'e', digits = 1)
      )
      

    return(tab_of_pars)

}

# by sex
PrintParameterTable(fit$sex09)
# by cohort
PrintParameterTable(fit$cohort)
# by age
PrintParameterTable(fit$origin)

# Parametric decomposition of mortality differences -------------

DecomposeFetoInfantDeaths <-
  function (filt_fit, pop1, pop2) {

    require(DemoDecomp)
    require(rlang)

    FetoinfantDeaths <-
      function (pars) { (1-FetoinfantSurv(76-24, pars))*1e5 }

    total_difference <-
      filt_fit %>%
      unnest_legacy(par_draws) %>%
      select(draw, stratum, name, value) %>%
      filter(
        stratum %in%
          c(as_name(enquo(pop1)), as_name(enquo(pop2)))
      ) %>%
      spread(stratum, value) %>%
      group_by(draw) %>%
      summarise(
        {{pop1}} := FetoinfantDeaths({{pop1}}),
        {{pop2}} := FetoinfantDeaths({{pop2}}),
        diff = {{pop2}}-{{pop1}},
        reldiff = diff/{{pop1}}
      ) %>%
      ungroup() %>%
      summarise(
        {{pop1}} := mean({{pop1}}),
        {{pop2}} := mean({{pop2}}),
        diff_avg = mean(diff),
        diff_se = sd(diff),
        diff_q025 = quantile(diff, 0.025),
        diff_q975 = quantile(diff, 0.975),
        reldiff_avg = mean(reldiff),
        reldiff_se = sd(reldiff),
        reldiff_q025 = quantile(reldiff, 0.025),
        reldiff_q975 = quantile(reldiff, 0.975)
      ) %>%
      pivot_longer(everything())

    parameter_decomp_draw <-
      filt_fit %>%
      unnest_legacy(par_draws) %>%
      select(draw, stratum, name, value) %>%
      filter(
        stratum %in% c(as_name(enquo(pop1)), as_name(enquo(pop2)))
      ) %>%
      spread(stratum, value) %>%
      group_by(draw) %>%
      mutate(
        contribution =
          horiuchi(
            FetoinfantDeaths,
            pars1 = {{pop1}},
            pars2 = {{pop2}},
            N = 1e2
          )
      )
    
    parameter_decomp_summary <-
      parameter_decomp_draw %>%
      group_by(name) %>%
      summarise(
        '{{pop1}}_avg' :=
          mean({{pop1}}),
        '{{pop2}}_avg' :=
          mean({{pop2}}),
        '{{pop1}}_se' := sd({{pop1}}),
        '{{pop2}}_se' := sd({{pop2}}),
        '{{pop1}}_ci025' := quantile({{pop1}}, 0.025),
        '{{pop2}}_ci975' := quantile({{pop1}}, 0.975),
        contribution_avg = mean(contribution),
        contribution_se = sd(contribution),
        contribution_ci025 = quantile(contribution, 0.025),
        contribution_ci975 = quantile(contribution, 0.975),
      )

    component_decomp <-
      parameter_decomp_draw %>%
      group_by(draw) %>%
      group_modify(~{
        tibble(
          level = .x %>% slice(1) %>% pull('contribution'),
          ontog = .x %>% slice(2) %>% pull('contribution'),
          trans = .x %>% slice(3:5) %>% pull('contribution') %>% sum()
        )  
      }) %>%
      ungroup() %>%
      summarise(
        level_avg = mean(level),
        level_se = sd(level),
        level_ci025 = quantile(level, 0.025),
        level_ci975 = quantile(level, 0.975),
        ontog_avg = mean(ontog),
        ontog_se = sd(ontog),
        ontog_ci025 = quantile(ontog, 0.025),
        ontog_ci975 = quantile(ontog, 0.975),
        trans_avg = mean(trans),
        trans_se = sd(trans),
        trans_ci025 = quantile(trans, 0.025),
        trans_ci975 = quantile(trans, 0.975),
      ) %>%
      pivot_longer(everything())
      
    list(
      diff = total_difference,
      para = parameter_decomp_summary,
      comp = component_decomp
    )

  }

# by sex
DecomposeFetoInfantDeaths(fit$sex09, Female, Male)

# by cohort
DecomposeFetoInfantDeaths(fit$cohort, `1989`, `1999`)
DecomposeFetoInfantDeaths(fit$cohort, `1999`, `2009`)

# by origin
DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic White`, `Non-Hispanic Black`)
DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic White`, `Hispanic`)
DecomposeFetoInfantDeaths(fit$origin, `Non-Hispanic Black`, `Hispanic`)
