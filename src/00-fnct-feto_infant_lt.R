# Event-Exposure aggregator ---------------------------------------

# general function to aggregate multistate transitions
# and exposures within age groups
AggregateStateTransitions <- function (
  df,
  t_in, d_in, t_out, d_out,
  breaks,
  wide = TRUE, drop0exp = TRUE,
  closed_left = TRUE,
  disable_input_checks = FALSE
) {

  require(dplyr)

  t_in = enquo(t_in); d_in = enquo(d_in);
  t_out= enquo(t_out); d_out = enquo(d_out)

  # input checks

  if (identical(disable_input_checks, FALSE)) {
    # check if all transition times are contained in
    # range of breaks
    t_range = c(min(pull(df, !!t_in)), max(pull(df, !!t_out)))
    breaks_range = range(breaks)
    if ( identical(closed_left, TRUE) ) {
      if (any(
        t_range[1] < breaks_range[1] |
        t_range[2] >= breaks_range[2]
      )) {
        stop(paste0(
          'Transition time outside range of breaks. Ensure that all t_ >='),
          breaks_range[1], ' and <', breaks_range[2]
        )
      }
    }
    if ( identical(closed_left, FALSE) ) {
      if (any(
        t_range[1] <= breaks_range[1] |
        t_range[2] > breaks_range[2]
      )) {
        stop(paste0(
          'Transition time outside range of breaks. Ensure that all t_ >'),
          breaks_range[1], ' and <=', breaks_range[2]
        )
      }
    }
  }

  # total number of age intervals
  J_ = length(breaks)-1
  # index of age intervals
  j_ = 1:J_
  # width of age intervals
  n_j_ = diff(breaks)
  # unique origin states
  k_in_ = unique(pull(df, !!d_in))
  # unique destination states
  k_out_ = unique(pull(df, !!d_out))

  # find the index of an interval defined by
  # <breaks> each element in <x> is contained in
  # returns NA if x outside breaks
  FindIntervalJ <-
    function (x, breaks, cl = closed_left) {
      if (identical(cl, TRUE)) {
        # [a, b)
        right = FALSE; lowest = FALSE
      } else {
        # (a, b] with [a0, b0]
        right = TRUE; lowest = TRUE
      }
      .bincode(
        x = x, breaks = breaks,
        right = right, include.lowest = lowest
      )
    }

  # 1. Aggregation

  # tabulate exits by age, origin and destination state
  W_k_tab <-
    df %>%
    select(t_out = !!t_out, d_in = !!d_in, d_out = !!d_out) %>%
    mutate(
      # add age interval index to each exit
      j = FindIntervalJ(pull(., t_out), breaks, closed_left),
    ) %>%
    # for each observed combination of
    # age and
    # origin state and
    # destination state...
    group_by(d_in, d_out, j) %>%
    summarise(
      # ...total number of exits
      W_k = n(),
      # total time lost in age due to exit
      Lw_k = sum(breaks[j+1]-t_out)
    ) %>%
    ungroup()

  # tabulate exits by age and origin state
  # based on prior tabulation on destination specific exits
  W_tab <-
    W_k_tab %>%
    # for each observed combination of
    # age and
    # origin state...
    group_by(j, d_in) %>%
    summarise(
      # ...total exits
      W = sum(W_k),
      # ...total time lost in interval due to exit
      Lw = sum(Lw_k)
    ) %>%
    ungroup() %>%
    # add rows for missing combinations
    # of age interval and origin state
    complete(
      d_in = k_in_, j = j_,
      fill = list(W = 0, Lw = 0)
    )

  # tabulate entries by age and state entered into
  Z_tab <-
    df %>%
    select(d_in = !!d_in, t_in = !!t_in) %>%
    mutate(
      j = FindIntervalJ(pull(., t_in), breaks, closed_left),
    ) %>%
    group_by(j, d_in) %>%
    summarise(
      # ...total entries
      Z = n(),
      # ...total entries right at start of interval
      Z0 = sum(t_in==breaks[j]),
      # ...total time lost in interval due to late-entry
      Lz = sum(t_in-breaks[j])
    ) %>%
    ungroup() %>%
    complete(
      d_in = k_in_, j = j_,
      fill = list(Z = 0, Z0 = 0, Lz = 0)
    )

  # tabulate concurrent entries and exits by interval
  ZW_tab <-
    df %>%
    select(t_in = !!t_in, t_out = !!t_out, d_in = !!d_in) %>%
    # aggregate individual level entry
    # and exit times into predefined age groups
    mutate(
      # add interval index to each entry
      j = FindIntervalJ(pull(., t_in), breaks, closed_left),
      # are entries and exits in same interval?
      zw = j == FindIntervalJ(pull(., t_out), breaks)
    ) %>%
    # for each combination of
    # state and
    # interval
    group_by(d_in, j) %>%
    summarise(
      # ...total concurrent entries and exits
      # there may be NAs in logic vector <zw> when
      # and entry or exit falls outside the range
      # of all intervals. as those cases don't have to
      # be counted na.rm=TRUE is applied
      ZW = sum(zw, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    complete(
      d_in = k_in_, j = j_,
      fill = list(ZW = 0)
    )

  # 2. Determine risk-sets and exposure times

  # exit counts for all possible combinations
  # of origin state, destination state and
  # age interval
  # intrastate transitions are 0 now
  # but are added later
  W_k_tab_complete <-
    W_k_tab %>%
    select(-Lw_k) %>%
    complete(
      d_in = k_in_, d_out = k_out_, j = j_,
      fill = list(W_k = 0)
    )

  # occurence-exposure table
  oe_tab <-
    bind_cols(W_tab, Z_tab[,-(1:2)], ZW_tab[,-(1:2)]) %>%
    mutate(
      x = breaks[j],
      n = n_j_[j]
    ) %>%
    # for each entry state...
    group_by(d_in) %>%
    mutate(
      # number of observations entering j via j-1
      # R_(j+1) = R_j + Z_j - W_j
      R = c(0, head(cumsum(Z) - cumsum(W), -1)),
      # population at risk at x_j
      P = R + Z0,
      # number of observations in j that did neither start
      # nor end during j
      Q = R - W + ZW,
      # number of observations entering j
      # that do not end during j
      U = Z - ZW,
      # total observation time during j
      O = Q*n + (Z + W - ZW)*n - Lz - Lw,
      # number of intrastate transitions
      I = Q + U,
    ) %>%
    ungroup() %>%
    left_join(W_k_tab_complete, by = c('d_in', 'j')) %>%
    # intrastate transitions
    mutate(
      W_k = ifelse(d_in == d_out, I, W_k)
    ) %>%
    select(orig = d_in, dest = d_out, j, x, n, Z, W, P, O, W_k)

  # drop intervals with 0 exposure
  if (identical(drop0exp, TRUE)) {
    oe_tab <-
      oe_tab %>%
      filter(O > 0)
  }

  # convert to wide format
  if (identical(wide, TRUE)) {
    oe_tab <-
      oe_tab %>%
      mutate(dest = paste0('to_', dest)) %>%
      spread(key = dest, value = W_k)
  }

  return(oe_tab)

}


# Rates and standard errors ---------------------------------------

# occurence exposure rate that returns 0 when
# exposure is 0
RobustMortalityRate <-
  function (deaths, exposures) {
    m <- deaths/exposures
    ifelse(is.nan(m), 0, m)
  }

StdErrRate <- function (events, exposures) {
  se <- sqrt(events)/exposures
  ifelse(is.nan(se), 0, se)
}

# FILT (Feto-infant lifetable) object -----------------------------

# Calculate a fetoinfant lifetable (FILT) across strata
# Constructs a FILT object
FILT <-
  function (
    df,
    t_in = 'entry_time', d_in = 'origin',
    t_out = 'exit_time', d_out = 'destination',
    fetal_state_name = 'fetus',
    neonate_state_name = 'neonate',
    postneonate_state_name = 'postneonate',
    death_state_name = 'death',
    censored_state_name = 'censored',
    breaks, stratum
  ) {
    
    require(dplyr)
    require(rlang)
    
    fetoinfant_life_table <-
      df %>%
      group_by({{stratum}}) %>%
      group_modify(~{
        
        fetoinfant_event_table <-
          AggregateStateTransitions(
            .x,
            {{t_in}}, {{d_in}}, {{t_out}}, {{d_out}},
            breaks = breaks,
            wide = TRUE, drop0exp = FALSE
          )
        
        from_fetus <-
          filter(fetoinfant_event_table,
                 orig == {{fetal_state_name}})
        from_neonate <-
          filter(fetoinfant_event_table,
                 orig == {{neonate_state_name}})
        from_postneonate <-
          filter(fetoinfant_event_table,
                 orig == {{postneonate_state_name}})
        
        fetoinfant_life_table <-
          tibble(
            # start of age interval
            x = from_fetus[['x']],
            # width of age interval
            n = from_fetus[['n']],
            # transition counts and exposure times
            # foetus
            N_F =
              from_fetus[['P']],
            E_F =
              from_fetus[['O']],
            D_F =
              from_fetus[[paste0('to_', {{death_state_name}})]],
            # neonates
            B =
              from_fetus[[paste0('to_', {{neonate_state_name}})]],
            N_N =
              from_neonate[['P']],
            E_N =
              from_neonate[['O']],
            D_N =
              from_neonate[[paste0('to_', {{death_state_name}})]],
            # postneonates
            N_P =
              from_postneonate[['P']],
            E_P =
              from_postneonate[['O']],
            D_P =
              from_postneonate[[paste0('to_', {{death_state_name}})]],
            # infants
            N_I =
              N_N + N_P,
            E_I =
              E_N + E_P,
            D_I =
              D_N + D_P,
            # survivors
            C =
              from_postneonate[[paste0('to_', {{censored_state_name}})]],
            # total
            N =
              N_F + N_I,
            E =
              E_F + E_I,
            D =
              D_F + D_I
          )
        
      })
    
    if (missing(stratum)) {
      filt <-
        fetoinfant_life_table %>%
        mutate(stratum = 'Total')
      attr(filt, 'stratum_name') <- NA
      attr(filt, 'stratum_levels') <- 'Total'
    } else {
      filt <-
        fetoinfant_life_table %>%
        rename(stratum = {{stratum}})
      attr(filt, 'stratum_name') <- as_name(enquo(stratum))
      attr(filt, 'stratum_levels') <- as.character(unique(filt$stratum))
    }
    
    class(filt) <-
      c('FILT', 'tbl_df', 'tbl', 'data.frame')
    attr(filt, 'breaks') <- breaks
    
    return(filt)
    
  }

# FILT methods ----------------------------------------------------


FILTCohortSize <- function (filt) {
  filt %>%
    group_by(stratum) %>%
    summarise(N = first(N))
}

FILTCohortSizeTotal <- function (filt) {
  FILTCohortSize(filt) %>%
    pull(N) %>% sum()
}

FILTTransitions <- function (filt) {
  filt %>%
    group_by(x) %>%
    summarise(
      D_F = sum(D_F),
      D_N = sum(D_N),
      D_P = sum(D_P),
      D_I = D_N + D_P,
      D = sum(D),
      B = sum(B),
      C = sum(C)
    )
}

FILTTransitionsTotal <- function (filt) {
  filt %>%
    summarise(
      D_F = sum(D_F),
      D_N = sum(D_N),
      D_P = sum(D_P),
      D_I = D_N + D_P,
      D = sum(D),
      B = sum(B),
      C = sum(C)
    )
}

FILTExposures <- function (filt) {
  filt %>%
    group_by(stratum) %>%
    summarise(
      E_F = sum(E_F),
      E_N = sum(E_N),
      E_P = sum(E_P),
      E_I = E_N + E_P,
      E = sum(E)
    )
}

FILTExposuresTotal <- function (filt) {
  filt %>%
    summarise(
      E_F = sum(E_F),
      E_N = sum(E_N),
      E_P = sum(E_P),
      E_I = E_N + E_P,
      E = sum(E)
    )
}

FILTRelativeExposures <- function (filt) {
  filt %>%
    group_by(stratum) %>%
    transmute(
      x = x,
      n = n,
      pE_F = E_F/E,
      pE_N = E_N/E,
      pE_P = E_P/E,
      pE_I = E_I/E
    ) %>%
    ungroup()
}

FILTAgeDistributionOfBirths <- function (filt) {
  filt %>%
    group_by(stratum) %>%
    transmute(x = x, n = n, pB = B/sum(B)) %>%
    ungroup()
}

FILTNeonatalFetalMortalityRatio <- function (filt) {
  filt %>%
    group_by(stratum) %>%
    transmute(
      x = x,
      n = n,
      m_N = D_N/E_N,
      m_F = D_F/E_F,
      NFMR = m_N/m_F,
      logNFMR = log(NFMR),
      logNFMR_se = sqrt(1/D_F+1/D_N)
    ) %>%
    ungroup()
}

FILTMortalityRates <- function (filt) {
  filt %>%
    group_by(stratum) %>%
    transmute(
      x = x,
      n = n,
      m_F = RobustMortalityRate(D_F, E_F),
      mstder_F = StdErrRate(D_F, E_F),
      m_N = RobustMortalityRate(D_N, E_N),
      mstder_N = StdErrRate(D_N, E_N),
      m_P = RobustMortalityRate(D_P, E_P),
      mstder_P = StdErrRate(D_P, E_P),
      m_I = RobustMortalityRate(D_I, E_I),
      mstder_I = StdErrRate(D_I, E_I),
      m = RobustMortalityRate(D, E),
      mstder = StdErrRate(D, E),
    ) %>%
    ungroup()
}

FILTSurvival <- function (filt) {
  filt %>%
    group_by(stratum) %>%
    transmute(
      x = x,
      n = n,
      p =
        exp(-RobustMortalityRate(D,E)*n),
      S_cnst_hzrd =
        head(cumprod(c(1, p)), -1),
      S_empirical =
        N/N[1],
      F_cnst_hzrd =
        1-S_cnst_hzrd,
      F_empirical =
        c(0, head(cumsum(D),-1))/N[1]
    ) %>%
    ungroup()
}

print.FILT <- function (FILT) {
  IntFormat <- function (x) formatC(x, format = 'd', big.mark = ',')

  N <- IntFormat(FILTCohortSizeTotal(FILT))
  
  exps <- FILTExposuresTotal(FILT)
  E <- IntFormat(exps[['E']])
  E_F <- IntFormat(exps[['E_F']])
  E_N <- IntFormat(exps[['E_N']])
  E_P <- IntFormat(exps[['E_P']])

  trns <- FILTTransitionsTotal(FILT)
  D <- IntFormat(trns[['D']])
  D_F <- IntFormat(trns[['D_F']])
  D_N <- IntFormat(trns[['D_N']])
  D_P <- IntFormat(trns[['D_P']])
  
  cat(
    'Strata\t: ', paste0(attr(FILT, 'stratum_levels'), collapse = ', '), '\n',
    'Initial cohort size\t: ', N, '\n',
    
    'Total exposure\t\t: ', E, '\n',
    '\t Fetus\t\t: ', E_F, '\n',
    '\t Neonatal\t: ', E_N, '\n',
    '\t Postneonatal\t: ', E_P, '\n',
    
    'Total deaths\t\t: ', D, '\n',
    '\t Fetus\t\t: ', D_F, '\n',
    '\t Neonatal\t: ', D_N, '\n',
    '\t Postneonatal\t: ', D_P, '\n',
    sep = ''
  )
  invisible(FILT)
}