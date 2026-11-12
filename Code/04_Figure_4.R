#------------------------------------------------------------------------------#
#                                                                              #
#                               1. Loading data                                #
#                                                                              #
#------------------------------------------------------------------------------#
packs<-c("tidyverse", "lmPerm", "infomapecology")

for (i in 1:length(packs)){
  
  if(!require(packs[i],character.only=TRUE)){
    install.packages(packs[i])
  }
  library(packs[i],character.only=TRUE)
}



plant_list<-readRDS("Data/Processed objects/plant_list.rds")
map_all<-readRDS("Data/Processed objects/map_all.rds")
samp_intensity3<-readRDS("Data/Processed objects/samp_intensity3.rds")
map_null<-readRDS("Data/Processed objects/map_null.rds")



#------------------------------------------------------------------------------#
#                                                                              #
#                        2. Estimating observed NMI                            #
#                                                                              #
#------------------------------------------------------------------------------#

#Function to analyze modularity agreement between two layers
analyze_modules <- function(x, plant_list) {
  # Filter data
  test <- x$modules[x$modules$node_name %in% plant_list, ]
  test <- test[test$layer_id == 1 | test$layer_id == 2, ]
  
  # Check if data exists
  if (nrow(test) == 0) {
    return(list(
      confusion_matrix = NULL,
      agreement_rate = NA,
      changes_count = NA,
      total_species = 0,
      agreement_percentage = NA,
      NMI = NA
    ))
  }
  
  # Create wide format and calculate metrics
  test_wide <- test %>%
    dplyr::select(node_name, module, layer_id) %>%
    pivot_wider(names_from = layer_id, values_from = module, names_prefix = "layer_")
  
  confusion_matrix <- table(Layer_1 = test_wide$layer_1, Layer_2 = test_wide$layer_2)
  agreement_rate <- mean(test_wide$layer_1 == test_wide$layer_2)
  changes <- sum(test_wide$layer_1 != test_wide$layer_2)
  
  list(
    confusion_matrix = confusion_matrix,
    agreement_rate = agreement_rate,
    changes_count = changes,
    total_species = nrow(test_wide),
    agreement_percentage = round(agreement_rate * 100, 2),
    NMI = NMI(confusion_matrix)
  )
}
create_agreement_summary <- function(results_list, samp_intensity3) {
  # Create summary table
  agreement_summary <- map_dfr(results_list, ~{
    data.frame(
      agreement_percentage = .x$agreement_percentage,
      changes_count = .x$changes_count,
      total_species = .x$total_species,
      NMI = .x$NMI
    )
  }, .id = "element_name")
  
  # Separate identifiers
  SI_mod <- samp_intensity3 %>%
    separate(period, into = c("p_ID", "period"), sep = "\\.", remove = TRUE)
  
  agreement_summary <- agreement_summary %>%
    separate(element_name, into = c("year", "site", "p_ID", "period"), sep = "\\.", remove = FALSE) %>%
    left_join(SI_mod, by = c("year", "site", "p_ID", "period"))
  
  # Add network ID
  agreement_summary$net_ID <- paste0(agreement_summary$year, ".", agreement_summary$site)
  
  # Join treatment info
  #treatment_info$year <- as.character(treatment_info$year)
  #agreement_summary <- agreement_summary %>%
  #  left_join(treatment_info, by = c("year", "site"))
  
  # Fix floating point precision
  agreement_summary$NMI_fixed <- ifelse(agreement_summary$NMI > 1 - 1e-10, 1, agreement_summary$NMI)
  agreement_summary$NMI_fixed <- round(agreement_summary$NMI_fixed, 10)
  
  return(agreement_summary)
}

#Observed modularity agreement analysis
results_list <- map(map_all, ~analyze_modules(.x, plant_list))



agreement_summary <- create_agreement_summary(results_list, samp_intensity3)


#------------------------------------------------------------------------------#
#                                                                              #
#                           3. Estimating null NMI                             #
#                                                                              #
#------------------------------------------------------------------------------#
#Null modularity agreement


null_results_list <- map_depth(map_null,2, ~analyze_modules(.x, plant_list)) # Takes a few minutes
null_results_list$`2014.A.1.morning`[1]

create_nested_agreement_summary <- function(nested_results_list, samp_intensity3) {
  # Flatten and tag each sublist with its parent name
  flat_summary <- imap_dfr(nested_results_list, function(sublist, parent_name) {
    map_dfr(sublist, ~{
      data.frame(
        agreement_percentage = .x$agreement_percentage,
        changes_count = .x$changes_count,
        total_species = .x$total_species,
        NMI = .x$NMI
      )
    }) %>% mutate(element_name = parent_name)
  })
  
  # Separate identifiers
  SI_mod <- samp_intensity3 %>%
    separate(period, into = c("p_ID", "period"), sep = "\\.", remove = TRUE)
  
  agreement_summary <- flat_summary %>%
    separate(element_name, into = c("year", "site", "p_ID", "period"), sep = "\\.", remove = FALSE) %>%
    left_join(SI_mod, by = c("year", "site", "p_ID", "period"))
  
  # Add network ID
  agreement_summary$net_ID <- paste0(agreement_summary$year, ".", agreement_summary$site)
  

  # Fix floating point precision
  agreement_summary$NMI_fixed <- ifelse(agreement_summary$NMI > 1 - 1e-10, 1, agreement_summary$NMI)
  agreement_summary$NMI_fixed <- round(agreement_summary$NMI_fixed, 10)
  
  return(agreement_summary)
}

null_results_list$`2014.A.1.morning`[1]

agreement_summary_null <- create_nested_agreement_summary(null_results_list, samp_intensity3)

agreement_summary_null<-
agreement_summary_null %>% split(.$element_name)

# Summarize null results

pmap(list(split(agreement_summary, agreement_summary$element_name), agreement_summary_null), function(obs, null) {
  obs$mean_null_NMI <- mean(null$NMI_fixed , na.rm=TRUE)
  obs$sd_null_NMI <- sd(null$NMI_fixed , na.rm=TRUE)
  obs$Z_NMI <- (obs$NMI_fixed - mean(null$NMI_fixed , na.rm=TRUE)) / sd(null$NMI_fixed , na.rm=TRUE)
  obs$P_NMI <- mean(abs(null$NMI_fixed - mean(null$NMI_fixed , na.rm=TRUE)) >= abs(obs$NMI_fixed - mean(null$NMI_fixed , na.rm=TRUE)), na.rm=TRUE)
  
  obs$mean_null_agreement_percentage <- mean(null$agreement_percentage , na.rm=TRUE)
  obs$sd_null_agreement_percentage <- sd(null$agreement_percentage , na.rm=TRUE)
  obs$Z_agreement_percentage <- (obs$agreement_percentage - mean(null$agreement_percentage) )/ sd(null$agreement_percentage , na.rm=TRUE)
  obs$significant_agreement <- mean(abs(null$agreement_percentage - mean(null$agreement_percentage, na.rm=TRUE)) >= abs(obs$agreement_percentage - mean(null$agreement_percentage, na.rm=TRUE)), na.rm=TRUE)
  
  obs$mean_null_changes_count <- mean(null$changes_count , na.rm=TRUE)
  obs$sd_null_changes_count <- sd(null$changes_count , na.rm=TRUE)
  obs$Z_changes_count <- (obs$changes_count - mean(null$changes_count , na.rm=TRUE)) / sd(null$changes_count , na.rm=TRUE)
  obs$P_changes_count <- mean(abs(null$changes_count  - mean(null$changes_count , na.rm=TRUE)) >= abs(obs$changes_count- mean(null$changes_count , na.rm=TRUE)), na.rm=TRUE)
  
  
  
  return(obs)
}) -> agreement_summary_with_significance
agreement_summary_with_significance <- bind_rows(agreement_summary_with_significance)
agreement_summary_with_significance$Z_NMI[which(is.nan(agreement_summary_with_significance$Z_NMI))]<-0

# Checking sinificance
which(agreement_summary_with_significance $P_NMI <= 0.05)
which(agreement_summary_with_significance $P_agreement_percentage <= 0.05)
which(agreement_summary_with_significance $P_changes_count <= 0.05)




#------------------------------------------------------------------------------#
#                                                                              #
#                            4. Permutation model                              #
#                                                                              #
#------------------------------------------------------------------------------#

mdata<-agreement_summary_with_significance

mdata$source<-paste0(mdata$year, ".", mdata$site) %>% as.factor()

#### Visualize observed vs null NMI

period_num <- as.numeric(relevel(as.factor(mdata$period), ref="morning"))


head(mdata)



mdata<-mdata %>% mutate(
  period = factor(period, levels = c("morning", "afternoon", "night")),
  year = as.numeric(year)
)

set.seed(22)
mod_perm <- aovp(Z_NMI ~ period+log(intensity) + total_species +year + Error(net_ID),
                 data = mdata, perm = "Prob", maxIter = 10000, Ca = 0.001)


summary(mod_perm)

mdata_sub1.1<-
  mdata %>%
  group_by(source) %>%
  filter(all(c("afternoon", "night") %in% period)) %>%
  ungroup()

mdata_sub1.1<-mdata_sub1.1[-which(mdata_sub1.1$period=="morning"),]
set.seed(22)
mod_perm1.1 <- aovp(Z_NMI ~ period+log(intensity) + total_species +year + Error(net_ID),
                    data = mdata_sub1.1, perm = "Prob", maxIter = 10000, Ca = 0.001)

summary(mod_perm1.1)


mdata_sub1.2<-
  mdata[which( mdata$source %in% names(which(( mdata$source %>% table)==3 ))),] 

set.seed(22)
mod_perm1.2 <- aovp(Z_NMI ~ period+log(intensity) + total_species +year + Error(net_ID),
                    data = mdata_sub1.2, perm = "Prob", maxIter = 10000, Ca = 0.001)

summary(mod_perm1.2)




#------------------------------------------------------------------------------#
#                                                                              #
#                                5. Plotting                                   #
#                                                                              #
#------------------------------------------------------------------------------#

png("Figures/Figure 4.png", width = 6000, height = 3000, res = 600)
par(mar = c(4, 4, 3, 1))
layout(matrix(1:2, ncol = 2), widths = c(3, 1))

X_nets <- case_when(
  mdata$period == "morning"   ~ 1,
  mdata$period == "afternoon" ~ 2,
  mdata$period == "night"     ~ 3
)

dielcol  <- c("gold", "#FEA501", "#35336b")
col_nets <- ifelse(mdata$period == "morning", dielcol[1],
                   ifelse(mdata$period == "afternoon", dielcol[2], dielcol[3]))

par(las = 1)
plot(x = 1:3,
     y = NULL,
     xlim = c(0.5, 3.5),
     ylim = c(-3, 3),
     xaxt = "n",
     xlab = "",
     ylab = "Similarity in module affiliation (Z score of NMI)",
     main = "",
     cex = 0,
     bty = "n")
axis(1, at = 1:3, labels = c("Morning", "Afternoon", "Night"))
abline(h = 0, lty = 2, col = "gray")

add_period_segments <- function(df,
                                year_col   = "year",
                                site_col   = "site",
                                period_col = "period",
                                z_col      = "Z_NMI",
                                col        = alpha("lightgray", 0.3),
                                lwd        = 1,
                                lty        = 2) {
  
  # Map period names to numeric order
  period_map <- c("morning" = 1, "afternoon" = 2, "night" = 3)
  df$period_num <- period_map[as.character(df[[period_col]])]
  
  # Split by year x site
  groups <- split(df, interaction(df[[year_col]], df[[site_col]]))
  
  # Draw segments
  for (g in groups) {
    if (nrow(g) > 1) {
      g <- g[order(g$period_num), ]
      segments(
        x0 = g$period_num[-nrow(g)],
        y0 = g[[z_col]][-nrow(g)],
        x1 = g$period_num[-1],
        y1 = g[[z_col]][-1],
        col = col,
        lwd = lwd,
        lty = lty
      )
    }
  }
}

add_period_segments(mdata)
points(X_nets, mdata$Z_NMI, pch = 19, col = alpha(col_nets, 0.2),
       cex = log(mdata$total_species))


#model predictions -------------------------------------------------------------
# Within-network model: network ID as a factor, i.e. the same comparison as the
# within-network block of the permutation ANOVA (year is constant within a
# network, so it is absorbed by net_ID, as in that block).

mod_data <- mdata[!is.na(mdata$Z_NMI), ]
mod_data$period <- factor(mod_data$period, levels = c("morning", "afternoon", "night"))
mod_data$net_ID <- factor(mod_data$net_ID)

lm_fe <- lm(Z_NMI ~ net_ID + period + log(intensity) + total_species, data = mod_data)
print(drop1(lm_fe, test = "F"))   # period SS should equal the within-network R Sum Sq from aovp

# Covariate values held fixed for the predictions
int_ref  <- exp(mean(log(mod_data$intensity)))   # mean on the log scale used in the model
rich_ref <- mean(mod_data$total_species)

# Predicted mean Z_NMI per period for the average network, at fixed covariates
pred_fe <- function(d) {
  d$net_ID <- factor(d$net_ID)
  m <- lm(Z_NMI ~ net_ID + period + log(intensity) + total_species, data = d)
  g <- expand.grid(net_ID = levels(d$net_ID),
                   period = factor(levels(mod_data$period),
                                   levels = levels(mod_data$period)))
  g$intensity     <- int_ref
  g$total_species <- rich_ref
  p <- predict(m, newdata = g)
  as.numeric(tapply(p, g$period, mean))
}

pred_grid <- data.frame(period = levels(mod_data$period),
                        x      = 1:3,
                        fit    = pred_fe(mod_data))

# Cluster bootstrap: resample whole networks, not rows
nets <- levels(mod_data$net_ID)
B    <- 1000
set.seed(123)
boot_mat <- replicate(B, {
  pick <- sample(nets, length(nets), replace = TRUE)
  d_b  <- do.call(rbind, lapply(seq_along(pick), function(k) {
    d <- mod_data[mod_data$net_ID == pick[k], ]
    d$net_ID <- paste0("b", k)        # a network drawn twice counts as two networks
    d
  }))
  tryCatch(suppressWarnings(pred_fe(d_b)), error = function(e) rep(NA_real_, 3))
})

cat("Failed resamples:", sum(is.na(boot_mat[1, ])), "of", B, "\n")

pred_grid$lower <- apply(boot_mat, 1, quantile, 0.025, na.rm = TRUE)
pred_grid$upper <- apply(boot_mat, 1, quantile, 0.975, na.rm = TRUE)
print(pred_grid)

segments(pred_grid$x, pred_grid$lower,
         pred_grid$x, pred_grid$upper,
         lwd = 2, col = "black")
points(pred_grid$x, pred_grid$fit, pch = 15, cex = 1.4)


#legend panel ------------------------------------------------------------------
par(mar = c(0, 0, 0, 0))
plot(x = NULL, y = NULL, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n",
     bty = "n")
legend(x = 0, y = 0.7,
       legend    = seq(min(mdata$total_species), max(mdata$total_species), length = 5) %>% round,
       pch       = 19,
       bty       = "n",
       pt.cex    = seq(min(mdata$total_species), max(mdata$total_species), length = 5) %>% round %>% log,
       col       = alpha("gray", 0.4),
       title     = "Plant richness",
       title.adj = 4)

dev.off()
#------------------------------------------------------------------------------#
#                                                                              #
#                6. Checking modularity significance                           #
#                                                                              #
#------------------------------------------------------------------------------#

map_null$`2014.A.1.morning`[2]

map_depth(map_all, 1, ~.x$L) %>% unlist() -> observed_L_values
map_depth(map_all, 1, ~.x$m) %>% unlist() -> observed_m_values
null_L_values <- map_depth(map_null, 2, ~.x$L)
null_L_values <-map_depth(null_L_values,1,~unlist(.x))
null_m_values <- map_depth(map_null, 2, ~.x$m)
null_m_values <-map_depth(null_m_values,1,~unlist(.x))

p_values_L <- mapply(function(obs, null) {
  mean(null <= obs, na.rm = TRUE)
}, observed_L_values, null_L_values) %>% as.data.frame()

mean_null_L <- mapply(function(null) {
  mean(null, na.rm = TRUE)
}, null_L_values) %>% as.data.frame()

p_values_m <- mapply(function(obs, null) {
  mean(null >= obs, na.rm = TRUE)
}, observed_m_values, null_m_values)%>% as.data.frame()

mean_null_m <- mapply(function(null) {
  mean(null, na.rm = TRUE)
}, null_m_values) %>% as.data.frame()



p_values <- data.frame(
  element_name = rownames(p_values_L),
  observed_L_values = observed_L_values,
  mean_null_L = mean_null_L$.,
  p_value_L = p_values_L$.,
  observed_m_values = observed_m_values,
  mean_null_m = mean_null_m$.,
  p_value_m = p_values_m$.
)

saveRDS(p_values, "Data/Processed objects/modularity_significance_pvalues.rds")



readRDS("Data/Processed objects/modularity_significance_pvalues.rds")->p_values
p_values

p_values <- p_values %>%
  separate(
    col = element_name,
    into = c("year", "site",  "period_code", "period"),
    sep = "\\.",
    remove = FALSE
  )

p_values$year<-as.numeric(p_values$year)
p_values

p_values$observed_L_values %>% mean
p_values$observed_L_values %>% sd

p_values$mean_null_L %>% mean
p_values$mean_null_L %>% sd

p_values$observed_m_values %>% mean
p_values$observed_m_values %>% sd

p_values$mean_null_m %>% mean
p_values$mean_null_m %>% sd

