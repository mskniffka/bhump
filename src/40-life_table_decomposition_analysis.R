# Decomposition analysis of combined feto-infant mortality
# over gestational age

# Init ------------------------------------------------------------

library(tidyverse)

# figure specs
source('src/00-figure_specifications.R')
# fetoinfant lifetable functions
source('src/00-fnct-feto_infant_lt.R')

# figures
fig <- list()
# tables
tab <- list()

# constants
cnst <- list(
  ar = 0.3,
  scaler = 1e5,
  x_breaks = seq(24, 50, 4),
  x_name = 'Week of gestation',
  x_limits = c(24, 50)
)

# Data ------------------------------------------------------------

# feto-infant lifetables
load('out/30-fetoinfant_lifetables.RData')

# Plot fetoinfant mortality trajectory ----------------------------

tab$total_filt <-
  filt$total09 %>%
  FILTMortalityRates()

tab$birth_distribution <-
  filt$total09 %>%
  FILTAgeDistributionOfBirths()
  
fig$gestational_age_pattern <-
  rbind(
  tab$total_filt,
  tab$total_filt %>% mutate(x = x+n-1e-9)
) %>%
  ggplot(aes(x = x)) +
  geom_step(
    aes(y = m*cnst$scaler)
  ) +
  geom_ribbon(
    aes(
      ymin = (m-2*mstder)*cnst$scaler,
      ymax = (m+2*mstder)*cnst$scaler
    ),
    alpha = 0.2
  ) +
  geom_area(
    aes(x = x, y = pB*100),
    data =
      bind_rows(
        tab$birth_distribution,
        tab$birth_distribution %>% mutate(x = x+n-1e-9)
      )
  ) +
  scale_x_continuous(
    name = cnst$x_name,
    breaks = seq(24, 77, 4)
  ) +
  scale_y_continuous(
    paste0('Feto-infant deaths per ',
           formatC(cnst$scaler, format = 'd', big.mark = ','),
           '\n person-weeks at risk')
  ) +
  fig_spec$MyGGplotTheme(ar = 0.7) +
  coord_cartesian(expand = FALSE)
fig$gestational_age_pattern

fig_spec$ExportPDF(
  fig$gestational_age_pattern,
  '40-gestational_age_pattern',
  'out',
  width = fig_spec$width,
  height = fig_spec$width*0.7
)

# Plot perinatal population dynamics ------------------------------

# feto-infant life-table during the perinatal period
tab$perinatal_filt <-
  filt$total09 %>%
  filter(x >= 24, x <= 50)
  
# relative exposures by weeks of gestation and
# state (fetus vs. neonate vs. postneonate)
tab$relative_exposures <-
  FILTRelativeExposures(tab$perinatal_filt) %>%
  select(x, n, pE_F, pE_N, pE_P) %>%
  pivot_longer(-c('x', 'n'))

# mortality rates by week of gestation and state
tab$mortality_rates <-
  tab$perinatal_filt %>%
  FILTMortalityRates() %>%
  rename(total_m = m, total_stder = mstder) %>%
  pivot_longer(
    -c('x', 'n', 'stratum', 'total_m', 'total_stder'),
    names_to = c('measure', 'population'),
    names_sep = '_'
  ) %>%
  pivot_wider(
    names_from = 'measure',
    values_from = 'value',
  ) %>%
  filter(population != 'I', m != 0) %>%
  mutate(
    name =
      factor(
        population,
        levels = c('F', 'N', 'P'),
        labels = c('Fetal mortality', 'Neonatal mortality',
                   'Postneonatal mortality')
      )
  )

# ratio of neonatal to postneonatal mortality
# by week of gestation
tab$nfmr <-
  tab$perinatal_filt %>%
  FILTNeonatalFetalMortalityRatio()

# plot of fetoinfant mortality and state distribution
fig$fetoinfant_mortality <-
  tab$mortality_rates %>%
  ggplot() +
  geom_area(
    aes(x = x, y = value*50, fill = name),
    data = tab$relative_exposures %>% bind_rows(mutate(., x=x+n-1e-9))
  ) +
  geom_ribbon(
    aes(x = x,
        ymin = (total_m-2*total_stder)*cnst$scaler,
        ymax = (total_m+2*total_stder)*cnst$scaler
    ),
    alpha = 0.2,
    data = tab$mortality_rates %>% bind_rows(mutate(., x=x+n-1e-9))
  ) +
  geom_step(
    aes(x = x,
        y = total_m*cnst$scaler
    ),
    size = fig_spec$line_size_m,
    data = tab$mortality_rates %>% bind_rows(mutate(., x=x+n-1e-9))
  ) +
  scale_fill_manual(
    values = fig_spec$discrete_colors_light
  ) +
  scale_x_continuous(
    name = cnst$x_name,
    breaks = cnst$x_breaks,
    limits = cnst$x_limits
  ) +
  scale_y_continuous(
    paste0('Feto-infant deaths per ',
           formatC(cnst$scaler, format = 'd', big.mark = ','),
           '\n person-weeks at risk')
  ) +
  fig_spec$MyGGplotTheme(ar = cnst$ar) +
  guides(fill = 'none') +
  coord_cartesian(expand = FALSE, clip = 'off')
fig$fetoinfant_mortality

# plot of state specific mortality rates
# by weeks of gestation
fig$separate_rates <-
  tab$mortality_rates %>%
  ggplot() +
  geom_ribbon(
    aes(x = x,
      ymin = (m-2*mstder)*cnst$scaler,
        ymax = (m+2*mstder)*cnst$scaler,
        fill = name),
    data = tab$mortality_rates %>% bind_rows(mutate(., x=x+n-1e-9))
  ) +
  geom_step(
    aes(x = x,
        y = m*cnst$scaler, color = name),
    size = fig_spec$line_size_m,
    data = tab$mortality_rates %>% bind_rows(mutate(., x=x+n-1e-9))
  ) +
  geom_step(
    aes(x = x, y = total_m*cnst$scaler),
    color = 'grey',
    size = fig_spec$line_size_m,
    data = tab$mortality_rates %>% bind_rows(mutate(., x=x+n-1e-9))
  ) +
  facet_wrap(~name) +
  scale_x_continuous(
    name = cnst$x_name,
    breaks = cnst$x_breaks,
    limits = cnst$x_limits
  ) +
  scale_y_continuous(
    paste0('Deaths per ',
           formatC(cnst$scaler, format = 'd', big.mark = ','),
           '\n person-weeks at risk'),
    trans = 'log10'
  ) +
  scale_color_manual(values = fig_spec$discrete_colors) +
  scale_fill_manual(values = fig_spec$discrete_colors_light) +
  fig_spec$MyGGplotTheme(ar = 0.9) +
  guides(color = 'none', fill = 'none') +
  coord_cartesian(expand = FALSE, clip = 'off')
fig$separate_rates

# plot neonatal to fetal mortality
# by week of gestation
fig$feto_neonate_ratio <-
  tab$nfmr %>%
  ggplot(
    aes(x = x+0.5*n)
  ) +
  geom_ribbon(
    aes(
      x = x,
      ymin = exp(logNFMR-2*logNFMR_se),
      ymax = exp(logNFMR+2*logNFMR_se)
    ),
    alpha = 0.2,
    data = bind_rows(tab$nfmr, tab$nfmr %>% mutate(x=x+n-1e-9))
  ) +
  geom_step(
    aes(x = x, y = NFMR),
    size = fig_spec$line_size_m,
    data = bind_rows(tab$nfmr, tab$nfmr %>% mutate(x=x+n-1e-9))
  ) +
  geom_hline(yintercept = 1, size = fig_spec$line_size_s) +
  scale_x_continuous(
    name = cnst$x_name,
    breaks = cnst$x_breaks,
    limits = cnst$x_limits
  ) +
  scale_y_log10(
    'Neonatal-fetal\nmortality ratio',
    labels = scales::label_comma(accuracy = 0.1)
  ) +
  fig_spec$MyGGplotTheme(ar = cnst$ar) +
  coord_cartesian(expand = FALSE, clip = 'off')
fig$feto_neonate_ratio

fig$perinatal_popdynamics <-
  cowplot::plot_grid(
  fig$fetoinfant_mortality, fig$separate_rates, fig$feto_neonate_ratio,
  nrow = 3, align = 'hv', labels = 'AUTO', axis = 'l'
)
fig$perinatal_popdynamics

fig_spec$ExportPDF(
  fig$perinatal_popdynamics,
  filename = '40-perinatal_popdynamics',
  path = 'out',
  width = fig_spec$width,
  height = fig_spec$width*1.1
)

# Decompose feto-infant mortality ---------------------------------

# Perform a Kitagawa decomposition of the difference in
# combined feto-infant mortality over weeks of gestation x1 and x2.
# This requires a FILT object with a single stratum.
DecomposeFetoInfantMortalityDynamics <-
  function(filt, x1, x2, scaler = 1) {
    
    rates_and_relative_exposures <-
      filt %>%
      {left_join(
        FILTMortalityRates(.),
        FILTRelativeExposures(.)
      )}
    
    # proportions and rates of population 1
    pop1 <-
      rates_and_relative_exposures %>%
      filter(x == x1) %>%
      {list(
        m = .$m*scaler,
        p_k =
          c(.$pE_F,
            .$pE_N,
            .$pE_P),
        m_k =
          c(.$m_F*scaler,
            .$m_N*scaler,
            .$m_P*scaler)
      )}
    
    # proportions and rates of population 2
    pop2 <-
      rates_and_relative_exposures %>%
      filter(x == x2) %>%
      {list(
        m = .$m*scaler,
        p_k =
          c(.$pE_F,
            .$pE_N,
            .$pE_P),
        m_k =
          c(.$m_F*scaler,
            .$m_N*scaler,
            .$m_P*scaler)
      )}
    
    # mortality contributions to overall difference
    contribution_m_k <- (pop1$p_k+pop2$p_k)/2*(pop2$m_k-pop1$m_k)
    # in case both m are 0 let the contribution be 0
    contribution_m_k <-
      ifelse(is.nan(contribution_m_k), 0, contribution_m_k)
    names(contribution_m_k) <- c('fetus', 'neonate', 'postneonate')
    # total mortality contribution
    contribution_m <- sum(contribution_m_k)
    
    # exposure contributions to overall difference
    contribution_p_k <- (pop1$m_k+pop2$m_k)/2*(pop2$p_k-pop1$p_k)
    # in case both p are 0 let the contribution be 0
    contribution_p_k <-
      ifelse(is.nan(contribution_p_k), 0, contribution_p_k)
    names(contribution_p_k) <- c('fetus', 'neonate', 'postneonate')
    # total exposure contribution
    contribution_p <- sum(contribution_p_k)
    
    results <-
      list(
        # total diff
        m1 = pop1$m,
        m2 = pop2$m,
        diff_m = pop2$m-pop1$m,
        reldiff_m = (pop2$m-pop1$m)/pop1$m,
        # mortality
        m1_k = pop1$m_k,
        m2_k = pop2$m_k,
        diff_m_k = pop2$m_k - pop1$m_k,
        reldiff_m_k = (pop2$m_k - pop1$m_k)/pop1$m_k,
        contribution_m_k = contribution_m_k,
        contribution_m = contribution_m,
        # population structure
        p1_k = pop1$p_k,
        p2_k = pop2$p_k,
        diff_p_k = pop2$p_k - pop1$p_k,
        reldiff_p_k = (pop2$p_k - pop1$p_k)/pop1$p_k,
        contribution_p_k = contribution_p_k,
        contribution_p = contribution_p
      )
    
    return(results)
    
  }

DecomposeFetoInfantMortalityDynamics(filt$total09, 25, 33, scaler = 1e5)
DecomposeFetoInfantMortalityDynamics(filt$total09, 33, 39, scaler = 1e5)
DecomposeFetoInfantMortalityDynamics(filt$total09, 39, 45, scaler = 1e5)
DecomposeFetoInfantMortalityDynamics(filt$total09, 45, 72, scaler = 1e5)
