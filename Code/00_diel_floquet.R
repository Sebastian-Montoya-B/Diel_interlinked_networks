# =============================================================================
# 00_diel_floquet.R  —  Floquet stability & perturbation spread, diel networks
#
# Notation follows
# Klausmeier CA (2008) Theor Ecol 1:153-161:
#
#   J_k    community matrix of diel period k, interactions only (zero diagonal)
#   A_k    = J_k - d*I, the full matrix with intraspecific self-regulation d
#   T_k    duration of period k (hours);  T = sum(T_k) = 24
#   E_k    = expm(A_k T_k)
#   X      = E_K ... E_2 E_1, the matrix carrying x through one diel cycle
#          (first period RIGHTMOST; Klausmeier Eq. 10)
#   rho_i  Floquet multipliers, eigenvalues of X
#   mu_i   Floquet exponents, log(rho_i)/T
#   s      = max Re(mu_i) at d = 0. Because every J_k has a zero diagonal,
#          det(X) = 1 and s >= 0 always: s is the self-regulation threshold, the
#          smallest d at which the community is locally stable. Higher s = less
#          stable.
#   C      push accumulated over one cycle; C_k = A_k^-1 (E_k - I)
#   R      press-response matrix (I - X)^-1 C, averaged over the K starting
#          boundaries of the cycle
#
# Input : long data frame, ONE ROW PER OBSERVATION:
#           network_ID, consumer, resource, period, int_type2
#         period    %in% {"morning","afternoon","night"}
#         int_type2 %in% {"mutualistic","antagonistic","dual"}
#         plant abundance table: network_ID, resource, abundance
#
# Main entry point: null_periods_par()  (parallel wrapper of null_periods()).
# =============================================================================

suppressMessages({
  has_expm <- requireNamespace("expm", quietly = TRUE)
  if (!has_expm) library(Matrix)
})
.expm <- function(M) if (has_expm) expm::expm(M) else as.matrix(Matrix::expm(Matrix::Matrix(M)))

# accept Int_type2, INT_TYPE2, network_id, ... without editing the data
.fixnames <- function(df) {
  canon  <- c("network_id","consumer","resource","period","int_type2")
  target <- c("network_ID","consumer","resource","period","int_type2")
  i <- match(tolower(names(df)), canon)
  names(df)[!is.na(i)] <- target[i[!is.na(i)]]
  df
}

# =============================================================================
# 1.  PARAMETERS
# =============================================================================

PARS <- list(
  ## cycle structure. A network lacking one day period (morning or afternoon)
  ## applies the sampled day period to the whole 16 h of daylight.
  phase_hours = c(morning = 8, afternoon = 8, night = 8),
  cycle_order = c("morning", "afternoon", "night"),
  
  ## unmeasured constants
  ant_scale     = 0.15,          # alpha_ant, fixed (only alpha_ant * e is identified)
  mut_ant_ratio = c(0.1, 10),    # r = alpha_mut / alpha_ant, log-uniform
  efficiency    = c(0.05, 0.5),  # e, log-uniform
  
  ## w, the coefficient of a dual visit in the plant rows. The null functions use
  ## min(dual_w_grid); set it to a single value per run (e.g. -0.15 or -0.01).
  dual_w_grid = -0.15,
  
  ## R is computed with d = s + d_margin, separately for each community and for
  ## each randomisation, so every community sits at the same distance above its
  ## own stability threshold.
  d_margin = 0.10
)

draw_pars <- function(P = PARS, dual_w = -P$ant_scale) {
  ratio <- exp(runif(1, log(P$mut_ant_ratio[1]), log(P$mut_ant_ratio[2])))
  list(ant = P$ant_scale, mut = P$ant_scale * ratio, ratio = ratio,
       eff = exp(runif(1, log(P$efficiency[1]), log(P$efficiency[2]))),
       dual_w = dual_w)
}

# =============================================================================
# 2.  DATA -> per-period count matrices
# =============================================================================

prepare_network <- function(df, net, P = PARS) {
  df <- .fixnames(df)
  d <- as.data.frame(df[df$network_ID == net, , drop = FALSE])
  for (v in c("consumer","resource","period","int_type2")) d[[v]] <- as.character(d[[v]])
  
  plants    <- sort(unique(d$resource))
  consumers <- sort(unique(d$consumer))
  periods   <- intersect(P$cycle_order, unique(d$period))
  
  grid <- function(sub) {
    m <- matrix(0, length(plants), length(consumers), dimnames = list(plants, consumers))
    if (nrow(sub) == 0) return(m)
    tb <- stats::xtabs(~ resource + consumer, data = sub)
    m[rownames(tb), colnames(tb)] <- as.matrix(tb); m
  }
  out <- lapply(periods, function(per) {
    s <- d[d$period == per, , drop = FALSE]
    list(mut  = grid(s[s$int_type2 == "mutualistic",  , drop = FALSE]),
         ant  = grid(s[s$int_type2 == "antagonistic", , drop = FALSE]),
         dual = grid(s[s$int_type2 == "dual",         , drop = FALSE]))
  })
  names(out) <- periods
  list(net = net, plants = plants, consumers = consumers, periods = periods, counts = out)
}

# Durations (hours) of the periods a network has. An unsampled day period is
# absorbed by the sampled one, so the day/night split of the cycle is preserved.
phase_T <- function(periods, P = PARS) {
  h <- P$phase_hours[periods]
  day_all <- c("morning", "afternoon")
  day_obs <- intersect(day_all, periods)
  if (length(day_obs) > 0 && length(day_obs) < length(day_all))
    h[day_obs] <- sum(P$phase_hours[day_all]) / length(day_obs)
  h * 24 / sum(h)
}

# Guild: a visitor with any antagonistic or dual record belongs to Ant.
consumer_guild <- function(prep) {
  g <- vapply(prep$consumers, function(cs) {
    a  <- sum(vapply(prep$counts, function(x) sum(x$ant[, cs]),  numeric(1)))
    dd <- sum(vapply(prep$counts, function(x) sum(x$dual[, cs]), numeric(1)))
    if (a > 0 || dd > 0) "antagonistic" else "mutualistic"
  }, character(1))
  names(g) <- prep$consumers; g
}

# =============================================================================
# 3.  Community matrix for one period (Appendix S6, eqns 13-14)
# =============================================================================
# Plant row   : A_ij = (alpha_mut M_ij - alpha_ant H_ij + w D_ij) / c_j,
#               c_j = total records of visitor j in the period
# Visitor row : A_ji = e V_ij / b_i,  b_i = plant abundance (always positive)

build_J <- function(cnt, plants, consumers, plant_ab, pd) {
  nodes <- c(plants, consumers)
  J <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  Ctot <- cnt$mut + cnt$ant + cnt$dual
  cons_ab <- colSums(Ctot)
  
  for (ip in seq_along(plants)) for (ic in seq_along(consumers)) {
    p <- plants[ip]; cs <- consumers[ic]
    if (Ctot[ip, ic] == 0) next
    if (is.null(plant_ab[[p]])) stop("no abundance for plant '", p, "'")
    
    if (cons_ab[ic] > 0)
      J[p, cs] <- (cnt$mut[ip, ic] * pd$mut - cnt$ant[ip, ic] * pd$ant +
                     cnt$dual[ip, ic] * pd$dual_w) / cons_ab[ic]
    if (plant_ab[[p]] > 0)
      J[cs, p] <- pd$eff * Ctot[ip, ic] / plant_ab[[p]]
  }
  J
}

# =============================================================================
# 4.  Perturbation-spread metrics (eqns 9-12)
# =============================================================================

# leak: mean |R| over all (mutualist, antagonist) visitor pairs
leakage <- function(R, consumers, guild) {
  poll <- consumers[guild[consumers] == "mutualistic"]
  herb <- consumers[guild[consumers] == "antagonistic"]
  if (!length(poll) || !length(herb)) return(NA_real_)
  mean(abs(R[poll, herb, drop = FALSE]))
}

# which visitors are active in which period
activity <- function(prep) {
  a <- vapply(prep$counts, function(x)
    colSums(x$mut + x$ant + x$dual) > 0, logical(length(prep$consumers)))
  matrix(a, nrow = length(prep$consumers),
         dimnames = list(prep$consumers, prep$periods))
}

# leak_cross: as leak, restricted to pairs never active in the same period
leakage_crossperiod <- function(R, prep, guild) {
  act  <- activity(prep)
  poll <- prep$consumers[guild[prep$consumers] == "mutualistic"]
  herb <- prep$consumers[guild[prep$consumers] == "antagonistic"]
  if (!length(poll) || !length(herb)) return(c(leak = NA_real_, n = 0))
  co <- act[poll, , drop = FALSE] %*% t(act[herb, , drop = FALSE])  # shared periods
  cross <- co == 0
  if (!any(cross)) return(c(leak = NA_real_, n = 0))
  sub <- R[poll, herb, drop = FALSE]
  c(leak = mean(abs(sub[cross])), n = sum(cross))
}

# =============================================================================
# 5.  Diagnostics on the input
# =============================================================================

inspect_input <- function(df, P = PARS) {
  df <- .fixnames(df)
  need <- c("network_ID","consumer","resource","period","int_type2")
  miss <- setdiff(need, names(df))
  cat("--- columns ---\n  found  :", paste(intersect(need, names(df)), collapse=", "), "\n")
  if (length(miss)) { cat("  MISSING:", paste(miss, collapse=", "), "<-- fix first\n"); return(invisible(FALSE)) }
  cat("\n--- period ---\n"); print(table(df$period, useNA="ifany"))
  bad_p <- setdiff(unique(df$period), names(P$phase_hours))
  if (length(bad_p)) cat("  UNRECOGNISED:", paste(bad_p, collapse=", "),
                         "\n  -> add to PARS$phase_hours AND PARS$cycle_order\n")
  cat("\n--- int_type2 ---\n"); print(table(df$int_type2, useNA="ifany"))
  bad_t <- setdiff(unique(df$int_type2), c("mutualistic","antagonistic","dual"))
  if (length(bad_t)) cat("  UNRECOGNISED:", paste(bad_t, collapse=", "), "\n")
  cat(sprintf("\n--- size ---\n  %d obs | %d networks | %d consumers | %d plants\n",
              nrow(df), length(unique(df$network_ID)),
              length(unique(df$consumer)), length(unique(df$resource))))
  ok <- !length(bad_p) && !length(bad_t)
  cat("\n", if (ok) "READY." else "NOT READY — fix the above.", "\n", sep="")
  invisible(ok)
}

check_data <- function(df, plant_abundance, P = PARS) {
  df <- .fixnames(df); plant_abundance <- .fixnames(as.data.frame(plant_abundance))
  cat("--- periods present ---\n"); print(table(df$network_ID, df$period))
  cat("\n--- guild composition (dual counts as antagonist) ---\n")
  tb <- do.call(rbind, lapply(unique(df$network_ID), function(n) {
    d <- df[df$network_ID == n, ]
    g <- tapply(d$int_type2, d$consumer, function(x)
      if (any(x %in% c("antagonistic","dual"))) "ant" else "mut")
    data.frame(network = n, obs = nrow(d),
               plants = length(unique(d$resource)), consumers = length(g),
               n_mut_spp = sum(g == "mut"), n_ant_spp = sum(g == "ant"),
               pct_ant_obs = round(100*mean(d$int_type2 %in% c("antagonistic","dual")), 1),
               pct_dual_obs = round(100*mean(d$int_type2 == "dual"), 1))
  }))
  print(tb, row.names = FALSE)
  if (any(tb$n_ant_spp == 0)) cat("\n  WARNING: networks with NO antagonist guild cannot yield leakage.\n")
  cat("\n--- plant abundance coverage ---\n")
  for (n in unique(df$network_ID)) {
    need <- unique(df$resource[df$network_ID == n])
    miss <- setdiff(need, plant_abundance$resource[plant_abundance$network_ID == n])
    cat(sprintf("  %-10s %d/%d%s\n", n, length(need)-length(miss), length(need),
                if (length(miss)) paste0("  MISSING: ", paste(miss, collapse=", ")) else ""))
  }
  invisible(tb)
}

# =============================================================================
# 6.  Sauve-style null model on selected diel period(s)
# =============================================================================
# For each randomisation, two independent permutations of the plant labels are
# applied to the rows of the selected period(s): one to the mutualistic matrix,
# one to the antagonistic matrix (dual counts travel with the antagonistic one).
# Other periods, plant abundances and every column are left unchanged.
#
# periods = "night", "afternoon", "morning" (or several). joint = TRUE applies
# the same pair of permutations to every selected period; it is irrelevant when
# a single period is selected. Networks lacking the selected period are skipped.

.permute_counts <- function(x, pM, pA, plants) {
  M <- x$mut [pM, , drop = FALSE]; rownames(M) <- plants
  A <- x$ant [pA, , drop = FALSE]; rownames(A) <- plants
  D <- x$dual[pA, , drop = FALSE]; rownames(D) <- plants   # dual travels with ant
  list(mut = M, ant = A, dual = D)
}

null_periods <- function(df, plant_abundance, n_null = 199, n_par_draws = 3,
                         periods = NULL, joint = TRUE,
                         P = PARS, seed = 1, networks = NULL, verbose = TRUE) {
  df <- .fixnames(df); plant_abundance <- .fixnames(as.data.frame(plant_abundance))
  if (identical(periods, "all")) periods <- NULL
  nets <- unique(df$network_ID)
  if (!is.null(networks)) nets <- intersect(nets, networks)
  preps <- lapply(nets, function(n) prepare_network(df, n, P)); names(preps) <- nets
  guild <- lapply(preps, consumer_guild)
  ab <- lapply(nets, function(n) {
    a <- plant_abundance[plant_abundance$network_ID == n, ]
    setNames(as.list(a$abundance), a$resource) }); names(ab) <- nets
  
  METRICS <- c("s_eff","leak","leak_cross")
  
  one_net <- function(n, pd) {
    pr <- preps[[n]]; g <- guild[[n]]
    sel <- if (is.null(periods)) pr$periods else intersect(pr$periods, periods)
    if (!length(sel)) return(NULL)              # network lacks the asked period(s)
    kperm <- which(pr$periods %in% sel)
    Tv <- phase_T(pr$periods, P); Tc <- sum(Tv); np <- length(pr$plants)
    nn <- length(pr$plants) + length(pr$consumers); I <- diag(nn)
    
    prep_period <- function(cn_k) {
      J <- build_J(cn_k, pr$plants, pr$consumers, ab[[n]], pd)
      list(J = J, E0 = NULL)
    }
    # periods that are never permuted: J and expm computed once per draw
    fixed <- lapply(seq_along(pr$periods), function(k)
      if (k %in% kperm) NULL else {
        z <- prep_period(pr$counts[[k]]); z$E0 <- .expm(z$J * Tv[k]); z })
    
    calc <- function(cnts) {
      pl <- fixed
      for (k in kperm) {
        z <- prep_period(cnts[[k]]); z$E0 <- .expm(z$J * Tv[k]); pl[[k]] <- z
      }
      J  <- lapply(pl, `[[`, "J")
      E0 <- lapply(pl, `[[`, "E0")
      
      # s: dominant Floquet exponent of the undamped cycle (d = 0)
      X <- I; for (k in seq_along(E0)) X <- E0[[k]] %*% X
      s <- log(max(Mod(eigen(X, only.values = TRUE)$values))) / Tc
      d <- s + P$d_margin                       # same distance above own threshold
      
      # R = (I - X)^-1 C, averaged over the K starting boundaries
      nph <- length(E0)
      Ek <- lapply(seq_len(nph), function(k) exp(-d * Tv[k]) * E0[[k]])
      Ck <- lapply(seq_len(nph), function(k) solve(J[[k]] - d*I, Ek[[k]] - I))
      Rsum <- matrix(0, nn, nn)
      for (sft in 0:(nph - 1)) {
        idx <- ((seq_len(nph) + sft - 1) %% nph) + 1
        Xd <- I; Kc <- matrix(0, nn, nn)
        for (k in idx) { Kc <- Ek[[k]] %*% Kc + Ck[[k]]; Xd <- Ek[[k]] %*% Xd }
        Rsum <- Rsum + solve(I - Xd, Kc)
      }
      R <- Rsum / nph
      dimnames(R) <- dimnames(J[[1]])
      
      lc <- leakage_crossperiod(R, pr, g)
      c(s_eff = s, leak = leakage(R, pr$consumers, g),
        leak_cross = unname(lc["leak"]))
    }
    
    obs_cnts <- vector("list", length(pr$periods))
    for (k in kperm) obs_cnts[[k]] <- pr$counts[[k]]
    obs <- calc(obs_cnts)
    
    nul <- vapply(seq_len(n_null), function(b) {
      cn <- vector("list", length(pr$periods))
      if (isTRUE(joint)) {
        pM <- sample(np); pA <- sample(np)       # same pair for every selected period
        for (k in kperm) cn[[k]] <- .permute_counts(pr$counts[[k]], pM, pA, pr$plants)
      } else {
        for (k in kperm) {
          pM <- sample(np); pA <- sample(np)     # independent pair per period
          cn[[k]] <- .permute_counts(pr$counts[[k]], pM, pA, pr$plants)
        }
      }
      calc(cn)
    }, numeric(length(METRICS)))
    
    do.call(rbind, lapply(METRICS, function(m) {
      o <- obs[[m]]; v <- nul[m, ]; v <- v[is.finite(v)]
      if (!is.finite(o) || length(v) < 10)
        return(data.frame(network_ID = n, metric = m, obs = o, null_mean = NA_real_,
                          null_sd = NA_real_, z = NA_real_, p_two = NA_real_, n_ok = length(v)))
      data.frame(network_ID = n, metric = m, obs = o, null_mean = mean(v),
                 null_sd = sd(v), z = (o - mean(v))/sd(v),
                 p_two = 2*min(mean(v <= o), mean(v >= o)), n_ok = length(v))
    }))
  }
  
  all_draws <- do.call(rbind, lapply(seq_len(n_par_draws), function(i) {
    set.seed(seed + i); pd <- draw_pars(P)
    pd$dual_w <- if (is.null(P$dual_w_grid)) -P$ant_scale else min(P$dual_w_grid)
    r <- do.call(rbind, lapply(nets, function(n) {
      if (verbose) cat(sprintf("  draw %d / network %s\n", i, n))
      one_net(n, pd) }))
    if (!is.null(r)) r$par_draw <- i
    r
  }))
  
  out <- .summarise_draws(all_draws)
  lab <- if (is.null(periods)) "all periods" else paste(periods, collapse = "+")
  out$permuted <- lab
  if (verbose) .print_summary(out, lab, n_null, n_par_draws)
  invisible(list(summary = out, draws = all_draws, permuted = lab, joint = joint))
}

# median z and range across parameter draws
.summarise_draws <- function(all_draws) {
  out <- do.call(rbind, lapply(
    split(all_draws, list(all_draws$network_ID, all_draws$metric), drop = TRUE),
    function(k) data.frame(network_ID = k$network_ID[1], metric = k$metric[1],
                           obs = median(k$obs, na.rm = TRUE),
                           z_median = median(k$z, na.rm = TRUE),
                           z_lo = min(k$z, na.rm = TRUE), z_hi = max(k$z, na.rm = TRUE),
                           p_median = median(k$p_two, na.rm = TRUE),
                           stringsAsFactors = FALSE)))
  rownames(out) <- NULL
  out[order(out$metric, out$network_ID), ]
}

.print_summary <- function(out, lab, n_null, n_par_draws) {
  cat(sprintf("\n--- Sauve null on %s: %d randomisations x %d parameter draws ---\n",
              lab, n_null, n_par_draws))
  cat("  z < 0: observed LOWER than chance rewiring of the permuted period(s).\n")
  cat("  For s_eff that means MORE stable; for leak/leak_cross, LESS spread.\n\n")
  for (m in unique(out$metric)) {
    k <- out[out$metric == m, ]
    cat(sprintf("  %s\n", m))
    print(format(k[, c("network_ID","z_median","z_lo","z_hi","p_median")], digits = 3),
          row.names = FALSE); cat("\n")
  }
}

# =============================================================================
# 7.  Parallel wrapper: one parameter draw per worker. Results are identical to
#     null_periods(..., n_par_draws = D, seed = s) run serially.
# =============================================================================
null_periods_par <- function(df, plant_abundance, n_null = 199, n_par_draws = 3,
                             periods = NULL, joint = TRUE, P = PARS, seed = 1,
                             networks = NULL, verbose = TRUE,
                             ncores = max(1L, parallel::detectCores() - 1L),
                             source_file = "Code/00_diel_floquet.R") {
  
  runner <- function(i) {
    r <- null_periods(df, plant_abundance, n_null = n_null, n_par_draws = 1L,
                      periods = periods, joint = joint, P = P,
                      seed = seed + i - 1L, networks = networks, verbose = FALSE)
    d <- r$draws
    if (!is.null(d)) d$par_draw <- i
    d
  }
  
  jobs <- seq_len(n_par_draws)
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, source, normalizePath(source_file))
    parallel::clusterExport(cl,
                            varlist = c("df","plant_abundance","n_null","periods","joint","P","seed","networks"),
                            envir = environment())
    out <- parallel::parLapply(cl, jobs, runner)
  } else {
    out <- parallel::mclapply(jobs, runner, mc.cores = ncores)
  }
  
  bad <- vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) warning(sum(bad), " draw(s) failed and were dropped")
  all_draws <- do.call(rbind, out[!bad])
  
  res <- .summarise_draws(all_draws)
  lab <- if (is.null(periods)) "all periods" else paste(periods, collapse = "+")
  res$permuted <- lab
  if (verbose) .print_summary(res, lab, n_null, n_par_draws)
  invisible(list(summary = res, draws = all_draws, permuted = lab, joint = joint))
}