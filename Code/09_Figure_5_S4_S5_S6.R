packs<-c("tidyverse")

for (i in 1:length(packs)){
  
  if(!require(packs[i],character.only=TRUE)){
    install.packages(packs[i])
  }
  library(packs[i],character.only=TRUE)
}


read.csv("Data/Processed objects/plant_abundance_F.csv") -> plant_ab
read.csv("Data/Processed objects/polli_red_F.csv") -> polli_red

polli_red$period<- gsub(".*\\.", "", polli_red$record_period)  

do.call(rbind, lapply(unique(polli_red$network_ID), function(n) {
  d <- polli_red[polli_red$network_ID == n, ]
  data.frame(network   = n,
             obs       = nrow(d),
             plants    = length(unique(d$resource)),
             consumers = length(unique(d$consumer)),
             matrix_dim= length(unique(d$resource)) + length(unique(d$consumer)),
             n_ant     = sum(d$Int_type2 == "antagonistic"),
             n_dual    = sum(d$Int_type2 == "dual"),
             periods   = paste(sort(unique(d$period)), collapse = "/"))
}))

source("Code/00_diel_floquet.R")
PARS$phase_hours <- c(morning = 8, afternoon = 8, night = 8) # Duration (hours) of diel periods

inspect_input(polli_red)
check_data(polli_red, plant_ab)

# Below we estimate stability and perturbation spread. This process may take several hours. 
#Reduce the number of nulls (currently 1000) if you want to do a quick test.

PARS$dual_w_grid <- -0.15
r_ni15   <-  null_periods_par(polli_red, plant_ab, 1000, 50, periods = "night", P = PARS, verbose = TRUE)
r_af15   <-  null_periods_par(polli_red, plant_ab, 1000, 50, periods = "afternoon", P = PARS, verbose = TRUE)
r_mo15   <-  null_periods_par(polli_red, plant_ab, 1000, 50, periods = "morning", P = PARS, verbose = TRUE)

#saveRDS(r_ni15, "Data/Processed objects/r_ni15.rds")
#saveRDS(r_af15, "Data/Processed objects/r_af15.rds")
#saveRDS(r_mo15, "Data/Processed objects/r_mo15.rds")


PARS$dual_w_grid <- -0.01
r_ni01   <-  null_periods_par(polli_red, plant_ab, 1000, 50, periods = "night", P = PARS, verbose = TRUE)
r_af01   <-  null_periods_par(polli_red, plant_ab, 1000, 50, periods = "afternoon", P = PARS, verbose = TRUE)
r_mo01   <-  null_periods_par(polli_red, plant_ab, 1000, 50, periods = "morning", P = PARS, verbose = TRUE)

#saveRDS(r_ni01, "Data/Processed objects/r_ni01.rds")
#saveRDS(r_af01, "Data/Processed objects/r_af01.rds")
#saveRDS(r_mo01, "Data/Processed objects/r_mo01.rds")
## ===========================================================================
## 1. Z scores by period, and the per-period diel shift
## ===========================================================================
full_data <- readRDS("Data/Processed objects/degree_unweighted.rds")
df_cor <- full_data[-which(full_data$year == 2015 | full_data$year == 2017), ] %>%
  select(year, site, period, z_score, n_species)
df_cor$network_ID <- paste(df_cor$year, df_cor$site, sep = ".")

w <- reshape(df_cor[, c("network_ID", "period", "z_score", "n_species")],
             idvar = "network_ID", timevar = "period", direction = "wide")
names(w) <- sub("z_score\\.[0-9]\\.",  "z_", names(w))
names(w) <- sub("n_species\\.[0-9]\\.", "n_", names(w))

## shift = previous period - current period, cycling morning -> afternoon -> night
w$shift_morning   <- w$z_night     - w$z_morning      # previous = night
w$shift_afternoon <- w$z_morning   - w$z_afternoon    # previous = morning
w$shift_night     <- w$z_afternoon - w$z_night        # previous = afternoon

w$diel_shift <- w$shift_night                         # kept for back-compatibility

## ===========================================================================
## 2. Plotting machinery 
## ===========================================================================
dielcol <- alpha(c("gold", "#FEA501", "#35336b"), 0.4)

metrics <- c(
  s_eff      = "Stability (Z score)",
  leak       = "Perturbation spread (Z score)",
  leak_cross = "Cross-period perturbation spread (Z score)"
)
flip_metric <- c(s_eff = TRUE, leak = FALSE, leak_cross = FALSE)
col_lab     <- c(s_eff = "A)", leak = "B)", leak_cross = "C)")

run_spec <- function(mo, af, ni) list(
  "Morning pattern contribution" =
    list(res = mo, col = dielcol[1], size = "n_morning",
         shift = "shift_morning",   shift_lab = "night - morning"),
  "Afternoon pattern contribution" =
    list(res = af, col = dielcol[2], size = "n_afternoon",
         shift = "shift_afternoon", shift_lab = "morning - afternoon"),
  "Night pattern contribution" =
    list(res = ni, col = dielcol[3], size = "n_night",
         shift = "shift_night",     shift_lab = "afternoon - night")
)

pfmt <- function(p) if (p < 0.001) "P < 0.001" else
  paste0("P = ", formatC(p, format = "f", digits = 3))

panel <- function(d, ylab, col, size_col, shift_col, cexf,
                  trend = "lowess", alpha_cut = 0.05, flip = FALSE) {
  d <- d[is.finite(d[[shift_col]]) & is.finite(d$z_median), ]
  x <- d[[shift_col]]
  y <- d$z_median
  if (isTRUE(flip)) y <- -y
  cx <- cexf(d[[size_col]])
  ct <- suppressWarnings(cor.test(y, x, method = "spearman"))
  
  plot(x, y, type = "n", bty = "n", ylab = ylab, xlab = "")
  abline(h = 0, lty = 2, col = "gray")
  
  if (trend != "none" && is.finite(ct$p.value) && ct$p.value < alpha_cut) {
    if (trend == "lm") {
      xs <- range(x); fit <- lm(y ~ x)
      lines(xs, coef(fit)[1] + coef(fit)[2] * xs, col = "grey60", lwd = 2)
    } else {
      lines(lowess(x, y, f = 1), col = "grey60", lwd = 2)
    }
  }
  
  ord <- order(-cx)
  points(x[ord], y[ord], col = col, pch = 19, cex = cx[ord])
  
  mtext(bquote(rho*" = "*.(formatC(ct$estimate, format = "f", digits = 2))*
                 ", "*.(pfmt(ct$p.value))), side = 3, cex = 0.8)
  
  invisible(data.frame(n = nrow(d),
                       median_z  = round(median(y), 3),
                       n_below_0 = sum(y < 0),
                       rho = round(unname(ct$estimate), 3),
                       p   = signif(ct$p.value, 3)))
}


make_fig <- function(runs, w, file, nets = NULL,
                     width = 6000, height = 5000, res = 600) {
  
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(dev.list())) graphics.off()
  
  if (is.null(nets)) nets <- unique(w$network_ID)
  w3 <- w[w$network_ID %in% nets, ]
  
  sz_rng <- range(sqrt(c(w3$n_morning, w3$n_afternoon, w3$n_night)), na.rm = TRUE)
  cexf <- function(n) {
    out <- 0.8 + 2.4 * (sqrt(n) - sz_rng[1]) / diff(sz_rng)
    out[!is.finite(out)] <- 0.8
    out
  }
  
  png(file, res = res, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  
  layout(matrix(1:9, ncol = 3, byrow = TRUE), widths = c(1.18, 1, 1))
  par(las = 1, oma = c(3, 0, 0, 0))
  
  tab <- list()
  for (rn in names(runs)) {
    dat <- merge(runs[[rn]]$res$summary, w, by = "network_ID")
    dat <- dat[dat$network_ID %in% nets, ]
    for (mn in names(metrics)) {
      first   <- mn == names(metrics)[1]
      top_row <- rn == names(runs)[1]
      par(mar = if (first) c(4, 9, 3, 1) else c(4, 5, 3, 1))
      
      out <- panel(dat[dat$metric == mn, ], metrics[[mn]],
                   runs[[rn]]$col, runs[[rn]]$size, runs[[rn]]$shift,
                   cexf, trend = "lowess", flip = flip_metric[[mn]])
      tab[[length(tab) + 1]] <- cbind(randomised = rn, metric = mn,
                                      shift = runs[[rn]]$shift_lab, out)
      
      if (first)
        mtext(rn, side = 2, line = 7, las = 0, font = 2, cex = 0.85)
      if (top_row)
        mtext(col_lab[[mn]], side = 3, line = 1.3, adj = 0, font = 2, cex = 0.85)
      if (top_row && mn == "leak_cross")
        legend("topright", legend = c(6, 15, 30), pch = 19,
               col = alpha("grey", 0.4), pt.cex = cexf(c(6, 15, 30)),
               bty = "n", cex = 0.9, y.intersp = 1.2, title = "Plant richness")
    }
  }
  mtext("Diel interlinked-pattern shift",
        side = 1, outer = TRUE, line = 1)
  
  tab <- do.call(rbind, tab)
  rownames(tab) <- NULL
  invisible(tab)
}

## ===========================================================================
## 3. The four figures
## ===========================================================================
r_mo15 <- readRDS("Data/Processed objects/r_mo15.rds")
r_af15 <- readRDS("Data/Processed objects/r_af15.rds")
r_ni15 <- readRDS("Data/Processed objects/r_ni15.rds")
r_mo01 <- readRDS("Data/Processed objects/r_mo01.rds")
r_af01 <- readRDS("Data/Processed objects/r_af01.rds")
r_ni01 <- readRDS("Data/Processed objects/r_ni01.rds")

runs15 <- run_spec(r_mo15, r_af15, r_ni15)
runs01 <- run_spec(r_mo01, r_af01, r_ni01)

## networks present in all three randomisations (sampled in all three periods)
nets3 <- Reduce(intersect, lapply(runs15, function(r) unique(r$res$summary$network_ID)))
cat("retained (", length(nets3), "):", paste(sort(nets3), collapse = ", "), "\n")
cat("dropped  :", paste(setdiff(unique(w$network_ID), nets3), collapse = ", "), "\n")

tab_15    <- make_fig(runs15, w, "Figures/Figure 5.png")
tab_01    <- make_fig(runs01, w, "Figures/Figure S4.png")
tab_15_3p <- make_fig(runs15, w, "Figures/Figure S5.png",   nets = nets3)
tab_01_3p <- make_fig(runs01, w, "Figures/Figure S6.png", nets = nets3)

print(tab_15,    row.names = FALSE)
print(tab_01,    row.names = FALSE)
print(tab_15_3p, row.names = FALSE)
print(tab_01_3p, row.names = FALSE)
