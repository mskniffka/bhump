# Aggregate Transitions Counts and Occupancy Times
#
# Episode-split-free Risk-set and Exposure Time Calculation
# from Event History Data
#
# @param df
#   A data frame.
# @param t_in
#   Entry time into state.
# @param d_in
#   State being entered.
# @param t_out
#   Exit time from state.
# @param d_out
#   State being exited into.
# @param breaks
#   A numeric vector of break points for time-scale.
# @param wide
#   Output table in wide format (default=TRUE)?
# @param closed_left
#   Time intervals closed to the left and open to the right (default=TRUE)?
# @param disable_input_checks
#   Should input checks be disabled (default=FALSE)?
#
# @return
#   A data frame with columns
#     orig: origin state
#     j:    age group index
#     x:    starting age of j
#     n:    width of j
#     Z:    number of entries into origin state during j
#     W:    number of exits from origin state during j
#     P:    population number in origin state at beginning of j
#     O:    total observation time of population visiting origin state in j
#     (if wide = FALSE)
#     dest: destination state
#     W_k:  number of exits from origin state to destination state during j
#     (if wide = TRUE)
#     to_*: number of exits from origin state to state * during j
AggregateStateTransitions <- function (
  df,
  t_in, d_in, t_out, d_out,
  breaks,
  wide = TRUE, drop0exp = TRUE,
  closed_left = TRUE,
  disable_input_checks = FALSE
) {

  require(tidyverse)

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
