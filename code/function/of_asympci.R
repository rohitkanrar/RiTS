# of_asympci.R
#
# O'Brien-Fleming Asymptotic Repeated Confidence Intervals
#
# Theory: At each pre-specified look time n_k, the cross-fit AIPW estimator
# satisfies the linear representation
#
#   Delta_hat_n_k - Delta = (1/n_k) * sum_{i=1}^{n_k} W_i + o(1/sqrt(n_k))
#
# where W_i = f_bar(Z_i) - Delta is an MDS. By the strong invariance principle
# (Strassen 1967, applied in Theorem 1 proof eq. A.10), the K-dimensional
# vector of standardized statistics
#
#   B_k = sqrt(n_k) * (Delta_hat_n_k - Delta) / sigma_hat_n_k
#
# converges jointly to MVN(0, Sigma) with Sigma_{jk} = sqrt(tau_j / tau_k),
# j <= k, which is the Brownian motion covariance. OBF critical values c_k are
# calibrated against this distribution, giving:
#
#   C_k^OBF = Delta_hat_n_k +/- c_k * sigma_hat_n_k / sqrt(n_k)
#
# with asymptotic familywise coverage 1 - alpha over all K looks.
#
# Notation follows asymp_cs.R throughout. Return structures are identical to
# get_asympcs() for drop-in compatibility with downstream code.


# --------------------------------------------------------------------------- #
#  get_of_critical_values
#    Compute OBF critical values via gsDesign for a given look schedule.
#
#  Arguments:
#    n_looks    : integer, total number of pre-specified looks K
#    alpha      : two-sided familywise error rate (already Bonferroni-adjusted
#                 for multiplicity across arms if called from get_of_asympci)
#    info_frac  : numeric vector of length n_looks, information fractions
#                 tau_1 < ... < tau_K = 1. Defaults to equally spaced fractions
#                 (1:n_looks)/n_looks, matching Appendix C.8 of the paper.
#
#  Returns:
#    Numeric vector of length n_looks: OBF critical values c_1 >= ... >= c_K.
#    For unknown variance (the asymptotic case here), these are normal-theory
#    bounds; no t-distribution conversion is applied because the method is
#    asymptotic (cf. how lyapunov_asympcs does not apply a t-correction either).
# --------------------------------------------------------------------------- #

get_of_critical_values <- function(n_looks,
                                   alpha,
                                   info_frac = NULL) {
  if (!requireNamespace("gsDesign", quietly = TRUE)) {
    stop("Package 'gsDesign' is required. Install it with install.packages('gsDesign').")
  }
  if (is.null(info_frac)) {
    info_frac <- (1:n_looks) / n_looks
  }
  if (length(info_frac) != n_looks) {
    stop("'info_frac' must have length equal to 'n_looks'.")
  }
  if (any(diff(info_frac) <= 0) || info_frac[n_looks] != 1) {
    stop("'info_frac' must be strictly increasing with last element equal to 1.")
  }

  # gsDesign uses one-sided alpha; divide by 2 for two-sided OBF CI.
  design <- gsDesign::gsDesign(
    k      = n_looks,
    n.I    = info_frac,
    test.type = 1,          # one-sided (we halve alpha to get two-sided bounds)
    alpha  = alpha / 2,
    sfu    = gsDesign::sfLDOF
  )
  return(design$upper$bound)
}


# --------------------------------------------------------------------------- #
#  get_of_asympci
#    Compute OBF asymptotic repeated CIs for all arms and contrasts.
#    Mirrors get_asympcs() in argument names, loop structure, and return value.
#
#  Arguments:
#    trt_hist     : integer vector of length N, treatment assignments in {1,...,K}
#    rwd_hist     : numeric vector of length N, observed rewards
#    prpns_mat    : N x K matrix of propensity scores (clipped, rows sum to 1)
#    context_hist : N x d matrix of covariates
#    placebo_arm  : integer, index of the reference/placebo arm (default 1)
#    times        : integer vector of look times n_1 < ... < n_K (= N at final)
#    aipw_master  : optional pre-computed list from get_aipw_seq(); if NULL it
#                   is computed internally (same default as get_asympcs)
#    alpha        : two-sided FWER before Bonferroni correction (default 0.05)
#    info_frac    : optional numeric vector of information fractions tau_1,...,tau_K.
#                   If NULL, defaults to times/max(times) (calendar fractions).
#                   For the theoretically correct information-based monitoring,
#                   pass hat_tau_k = n_k * sigma_hat_k^2 / (n_K * sigma_hat_K^2),
#                   estimated after a pilot run.
#    c_ks_contr   : optional pre-computed OBF critical values for contrasts
#                   (length = length(times)). If NULL, computed via gsDesign with
#                   Bonferroni correction alpha/(K-1). Overrides info_frac for
#                   contrasts when provided.
#    c_ks_ate     : optional pre-computed OBF critical values for per-arm ATEs
#                   (length = length(times)). If NULL, computed via gsDesign with
#                   Bonferroni correction alpha/K.
#    n_cores      : passed to get_aipw_seq (not used in CI computation itself)
#    learner      : regression learner passed to get_aipw_seq, one of
#                   c("full_ridge", "main_ridge", NULL)
#
#  Returns:
#    list(ate_seq, contr_seq) — identical structure to get_asympcs():
#      ate_seq   : array [length(times) x K x 3], dimnames Times/Arms/CI
#      contr_seq : array [length(times) x (K-1) x 3], dimnames Times/Arms/CI
#    CI dimension: c("Center", "Lower", "Upper")
# --------------------------------------------------------------------------- #

get_of_asympci <- function(trt_hist, rwd_hist, prpns_mat, context_hist,
                           placebo_arm, times,
                           aipw_master  = NULL,
                           alpha        = 0.05,
                           info_frac    = NULL,
                           c_ks_contr   = NULL,
                           c_ks_ate     = NULL,
                           n_cores      = 1,
                           learner      = "full_ridge") {

  N       <- length(trt_hist)
  K       <- length(unique(trt_hist))
  n_looks <- length(times)

  # ---- information fractions (default: calendar time) ---- #
  if (is.null(info_frac)) {
    info_frac <- times / max(times)
  }
  if (info_frac[n_looks] != 1) {
    # Normalise so final fraction is exactly 1, as gsDesign requires
    info_frac <- info_frac / info_frac[n_looks]
  }

  # ---- OBF critical values ---- #
  # Bonferroni-correct alpha for K-1 contrasts and K per-arm ATEs separately,
  # matching the alpha/(K-1) and alpha/K used in get_asympcs().
  if (is.null(c_ks_contr)) {
    c_ks_contr <- get_of_critical_values(
      n_looks   = n_looks,
      alpha     = alpha / (K - 1),    # Bonferroni over K-1 contrasts
      info_frac = info_frac
    )
  }
  if (is.null(c_ks_ate)) {
    c_ks_ate <- get_of_critical_values(
      n_looks   = n_looks,
      alpha     = alpha / K,          # Bonferroni over K arms
      info_frac = info_frac
    )
  }

  # ---- AIPW pseudo-outcomes at each look ---- #
  # Reuse pre-computed master if supplied, otherwise compute.
  # Note: get_asympcs always recomputes; we respect the supplied value here.
  if (is.null(aipw_master)) {
    aipw_master <- get_aipw_seq(
      y          = rwd_hist,
      X          = context_hist,
      treatment  = trt_hist,
      propensity = prpns_mat,
      times      = times,
      learner    = learner,
      n_cores    = n_cores
    )
  }

  # ---- output arrays (same structure as get_asympcs) ---- #
  ate_seq   <- array(NA, dim = c(n_looks, K,     3))
  contr_seq <- array(NA, dim = c(n_looks, K - 1, 3))

  # ---- helper: OBF CI at one look ---- #
  # Constructs est +/- c_k * hat_sigma / sqrt(n_k).
  # hat_sigma is the sample SD of the pseudo-outcomes, matching the cross-fit
  # variance estimator sigma_hat^2_n = var_hat_n(f_hat; a) in eq. (A.8).
  of_ci_from_pseudo <- function(pseudo_vec, c_k) {
    n_k     <- length(pseudo_vec)
    est     <- mean(pseudo_vec)
    se      <- sd(pseudo_vec) / sqrt(n_k)     # sigma_hat_n / sqrt(n_k)
    half_w  <- c_k * se
    c(est, est - half_w, est + half_w)
  }

  # ---- main loop over arms (mirrors get_asympcs loop structure) ---- #
  k <- 1   # index into contr_seq columns (K-1 active arms)

  for (trt_arm in 1:K) {

    # -- per-arm ATE CIs at every look -- #
    for (i in seq_along(times)) {
      time      <- times[i]
      aipw_ate  <- aipw_master[[paste(time)]][1:as.numeric(time), trt_arm]
      ate_seq[i, trt_arm, ] <- of_ci_from_pseudo(aipw_ate, c_ks_ate[i])
    }

    # skip contrast computation for placebo arm itself
    if (trt_arm == placebo_arm) {
      k <- k + 1
      next
    }

    # -- contrast (arm trt_arm - placebo) CIs at every look -- #
    for (i in seq_along(times)) {
      time        <- times[i]
      aipw_all    <- aipw_master[[paste(time)]][1:as.numeric(time), ]
      aipw_contr  <- aipw_all[, trt_arm] - aipw_all[, placebo_arm]
      contr_seq[i, k - 1, ] <- of_ci_from_pseudo(aipw_contr, c_ks_contr[i])
    }

    k <- k + 1
  }

  # ---- dimnames (identical to get_asympcs) ---- #
  dimnames(ate_seq) <- list(
    Times = times,
    Arms  = 1:K,
    CI    = c("Center", "Lower", "Upper")
  )
  dimnames(contr_seq) <- list(
    Times = times,
    Arms  = sort(setdiff(1:K, placebo_arm)),
    CI    = c("Center", "Lower", "Upper")
  )

  return(list(ate_seq, contr_seq))
}


# --------------------------------------------------------------------------- #
#  add_of_asympci
#    Wrapper mirroring add_asympcs(): attaches OBF CIs to a single trial output
#    object and stores them as out$ate_of and out$contr_of.
# --------------------------------------------------------------------------- #

add_of_asympci <- function(out, ate_start, n_looks = 30, placebo_arm = 1,
                           alpha = 0.05, info_frac = NULL,
                           learner = "full_ridge", force_compute = FALSE) {
  # n_looks = 30 matches the default in add_standard_ci() so that both OBF
  # methods operate on identical look times and are directly comparable.
  if (is.null(out$reward_benf)) {
    reward <- out$reward
  } else {
    reward <- out$reward_benf
  }
  trt <- out$trt
  N   <- length(trt)

  # Mirrors the times construction in add_standard_ci() exactly.
  times_ <- floor(seq(ate_start, N, length.out = n_looks))

  if (is.null(out$ate_of) || is.null(out$contr_of) || force_compute) {
    of_ci <- get_of_asympci(
      trt_hist     = trt,
      rwd_hist     = reward,
      prpns_mat    = out$log_dat$prpns_mat,
      context_hist = out$log_dat$context,
      placebo_arm  = placebo_arm,
      times        = times_,
      alpha        = alpha,
      info_frac    = info_frac,
      learner      = learner
    )
    out[["ate_of"]]   <- of_ci[[1]]
    out[["contr_of"]] <- of_ci[[2]]
  }

  out[["alpha"]]     <- alpha
  out[["ate_start"]] <- ate_start
  return(out)
}


# --------------------------------------------------------------------------- #
#  add_of_asympci_sim
#    Wrapper mirroring add_asympcs_sim(): applies add_of_asympci() to a list of
#    simulation replicates, optionally in parallel.
# --------------------------------------------------------------------------- #

add_of_asympci_sim <- function(out_list, ate_start, n_looks = 30, placebo_arm = 1,
                               alpha = 0.05, info_frac = NULL,
                               learner = "full_ridge",
                               n_cores = 1, force_compute = FALSE) {
  out_list <- parallel::mclapply(out_list, function(out) {
    if (is.null(out$ate_of) || is.null(out$contr_of) || force_compute) {
      out <- add_of_asympci(
        out         = out,
        ate_start   = ate_start,
        n_looks     = n_looks,
        placebo_arm = placebo_arm,
        alpha       = alpha,
        info_frac   = info_frac,
        learner     = learner
      )
    }
    return(out)
  }, mc.cores = n_cores)

  return(out_list)
}
