#------------------------------------------------------------------------------#
#                                                                              #
#                 1. Loading packages and sourcing functions                   #
#                                                                              #
#------------------------------------------------------------------------------#


# Load the packages and helping functions

packs<-c("igraph", 
         "performance",
         "magrittr",
         "scales", 
         "tidyverse",
         "tnet",
         "dendextend",
         "glmmTMB",
         "DHARMa",
         "emmeans",
         "performance",
         "ppcor")

for (i in 1:length(packs)){
  
  if(!require(packs[i],character.only=TRUE)){
    install.packages(packs[i])
  }
  library(packs[i],character.only=TRUE)
}

source("Code/00_helping_functions.R")


#------------------------------------------------------------------------------#
#                                                                              #
#                              2. Fitting GLMMS                                #
#                                                                              #
#------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------
full_data<-readRDS("Data/Processed objects/degree_unweighted_dual.rds")
full_data$size_centered <- scale((full_data$size%>% log), center = TRUE, scale = TRUE)
full_data$log_intensity_centered <- scale(log(full_data$intensity), 
                                          center = TRUE, scale = TRUE)
full_data$log_size_centered <- scale(log(full_data$size), center = TRUE, scale = TRUE)
full_data$source<-paste0(full_data$year, full_data$site)
full_data$yn<-scale(as.numeric(full_data$year))

#fd_subset1.1<-full_data[which(full_data$source %in% names(which((full_data$source %>% table)>=2 ))),] 

fd_subset1.1<-
  full_data %>%
  group_by(source) %>%
  filter(all(c("2.afternoon", "3.night") %in% period)) %>%
  ungroup()
fd_subset1.1<-fd_subset1.1[-which(fd_subset1.1$period=="1.morning"),]

fd_subset1.1$size_centered <- scale((fd_subset1.1$size%>% log), center = TRUE, scale = TRUE)
fd_subset1.1$log_intensity_centered <- scale(log(fd_subset1.1$intensity), 
                                             center = TRUE, scale = TRUE)
fd_subset1.1$log_size_centered <- scale(log(fd_subset1.1$size), center = TRUE, scale = TRUE)
fd_subset1.1$source<-paste0(fd_subset1.1$year, fd_subset1.1$site)
fd_subset1.1$yn<-scale(as.numeric(fd_subset1.1$year))


fd_subset1.2<-full_data[which(full_data$source %in% names(which((full_data$source %>% table)==3 ))),] 
fd_subset1.2$size_centered <- scale((fd_subset1.2$size%>% log), center = TRUE, scale = TRUE)
fd_subset1.2$log_intensity_centered <- scale(log(fd_subset1.2$intensity), 
                                             center = TRUE, scale = TRUE)
fd_subset1.2$log_size_centered <- scale(log(fd_subset1.2$size), center = TRUE, scale = TRUE)
fd_subset1.2$source<-paste0(fd_subset1.2$year, fd_subset1.2$site)
fd_subset1.2$yn<-scale(as.numeric(fd_subset1.2$year))


#--------------------------------------------------------------------------------------
full_data2<-readRDS("Data/Processed objects/degree_weighted_dual.rds")
full_data2$size_centered <- scale((full_data2$size %>% log), center = TRUE, scale = TRUE)
full_data2$log_intensity_centered <- scale(log(full_data2$intensity), 
                                           center = TRUE, scale = TRUE)

full_data2$log_size_centered <- scale(log(full_data2$size), center = TRUE, scale = TRUE)
full_data2$source<-paste0(full_data2$year, full_data2$site)
full_data2$yn<-scale(as.numeric(full_data2$year))


fd_subset2.1<-
  full_data2 %>%
  group_by(source) %>%
  filter(all(c("2.afternoon", "3.night") %in% period)) %>%
  ungroup()
fd_subset2.1<-fd_subset2.1[-which(fd_subset2.1$period=="1.morning"),]
fd_subset2.1$size_centered <- scale((fd_subset2.1$size%>% log), center = TRUE, scale = TRUE)
fd_subset2.1$log_intensity_centered <- scale(log(fd_subset2.1$intensity), 
                                             center = TRUE, scale = TRUE)
fd_subset2.1$log_size_centered <- scale(log(fd_subset2.1$size), center = TRUE, scale = TRUE)
fd_subset2.1$source<-paste0(fd_subset2.1$year, fd_subset2.1$site)
fd_subset2.1$yn<-scale(as.numeric(fd_subset2.1$year))

fd_subset2.2<-full_data2[which(full_data2$source %in% names(which((full_data2$source %>% table)==3 ))),] 
fd_subset2.2$size_centered <- scale((fd_subset2.2$size%>% log), center = TRUE, scale = TRUE)
fd_subset2.2$log_intensity_centered <- scale(log(fd_subset2.2$intensity), 
                                             center = TRUE, scale = TRUE)
fd_subset2.2$log_size_centered <- scale(log(fd_subset2.2$size), center = TRUE, scale = TRUE)
fd_subset2.2$source<-paste0(fd_subset2.2$year, fd_subset2.2$site)
fd_subset2.2$yn<-scale(as.numeric(fd_subset2.2$year))

###----------------------------------------------------------------------------------

## Fitting a model for unweighted degree

fd<-full_data


#night_dummy <- fd$night_dummy




m_centered <- glmmTMB(z_score ~ # z score of correlation in consumer diversity based on unweighted degree
                        period + # diel period
                        size_centered + # plant richness
                        log_intensity_centered + # sampling intensity
                        (1 | site/source)+ # Site/Community
                        yn, # Year
                      data = fd,
                      family="gaussian")

summary(m_centered)
simulateResiduals(m_centered, plot=TRUE, n=3000)
check_collinearity(m_centered)


sot<-emmeans(m_centered, pairwise~period, type="response") %>% 
  summary(., infer = c(TRUE, TRUE))



m_centered_sub1.1 <- glmmTMB(z_score ~ # z score of correlation in consumer diversity based on unweighted degree
                               period + # diel period
                               size_centered + # plant richness
                               log_intensity_centered + # sampling intensity
                               (1 | site/source)+ # Site/Community
                               yn, # Year
                             data = fd_subset1.1,
                             family="gaussian")

summary(m_centered_sub1.1)
simulateResiduals(m_centered_sub1.1, plot=TRUE)
check_collinearity(m_centered_sub1.1) # adj.VIF for categorical predictor looks fine.


sot_sub1.1<-emmeans(m_centered_sub1.1, pairwise~period, type="response") %>% 
  summary(., infer = c(TRUE, TRUE))


m_centered_sub1.2 <- glmmTMB(z_score ~ # z score of correlation in consumer diversity based on unweighted degree
                               period + # diel period
                               size_centered + # plant richness
                               log_intensity_centered + # sampling intensity
                               (1 | site/source)+ # Site/Community
                               yn, # Year
                             data = fd_subset1.2,
                             family="gaussian")

summary(m_centered_sub1.2)
simulateResiduals(m_centered_sub1.2, plot=TRUE, n=1000)
check_collinearity(m_centered_sub1.2) # adj.VIF for categorical predictor looks fine.


sot_sub1.2<-emmeans(m_centered_sub1.2, pairwise~period, type="response") %>% 
  summary(., infer = c(TRUE, TRUE))



## Fitting a model for weighted degree

fd2<-full_data2



m_centered2 <- glmmTMB(z_score ~ 
                         period +
                         size_centered + 
                         log_intensity_centered + 
                         (1 | site/source)+
                         yn,
                       data = fd2,
                       family="gaussian")

summary(m_centered2)
simulateResiduals(m_centered2, plot=TRUE, n=1000)
check_collinearity(m_centered2)

sot2<-emmeans(m_centered2, pairwise~period, type="response") %>%
  summary(., infer = c(TRUE, TRUE))

m_centered_sub2.1 <- glmmTMB(z_score ~ 
                               period +
                               size_centered + 
                               log_intensity_centered + 
                               (1 | site/source)+
                               yn,
                             data = fd_subset2.1,
                             family="gaussian")

summary(m_centered_sub2.1 )
simulateResiduals(m_centered_sub2.1 , plot=TRUE, n=1000)
check_collinearity(m_centered_sub2.1 ) # adj.VIF for categorical predictor looks fine.

sot_subset2.1<-emmeans(m_centered_sub2.1 , pairwise~period, type="response") %>%
  summary(., infer = c(TRUE, TRUE))


m_centered_sub2.2 <- glmmTMB(z_score ~ 
                               period +
                               size_centered + 
                               log_intensity_centered + 
                               (1 | site/source)+
                               yn,
                             data = fd_subset2.2,
                             family="gaussian")

summary(m_centered_sub2.2 )
simulateResiduals(m_centered_sub2.2 , plot=TRUE, n=1000)
check_collinearity(m_centered_sub2.2 ) # adj.VIF for categorical predictor looks fine.

sot_subset2.2<-emmeans(m_centered_sub2.2 , pairwise~period, type="response") %>%
  summary(., infer = c(TRUE, TRUE))


#------------------------------------------------------------------------------#
#                                                                              #
#                               3. Plotting                                    #
#                                                                              #
#------------------------------------------------------------------------------#

# Creating a couple of functions to facilitate the plotting process
plot_assortativ<-function(fd, sot, main, ylab="Correlation in consumer diversity (Z score of ρ)"){
  par(las=1, bty="n")
  plot(x=c(1:3), 
       y=sot$emmean,
       xlim=c(0.5,3.5),
       ylim= c(-2.2,3.5),
       xaxt="n",
       xlab="",
       ylab=ylab,
       main=main,
       cex=0)
  axis(1, at=1:3 ,labels= c("Morning", "Afternoon", "Night"))
  abline(h=0, lty=2, col="gray")
  x_nets<-c(1:3)
  
  jitter_x_by_density <- function(x, y, categories, factor = 0.1, bandwidth = NULL, max_spread = 0.4) {
    # Check for non-finite values in y
    if (any(!is.finite(y))) {
      stop("Y values contain non-finite values (NA, NaN, or Inf). Please clean the data.")
    }
    
    # Split data by categories
    unique_cats <- unique(categories)
    jittered_x <- numeric(length(x))
    
    for (cat in unique_cats) {
      # Get indices for the current category
      idx <- which(categories == cat)
      y_sub <- y[idx]
      
      # Compute density for this category's Y values
      if (is.null(bandwidth)) {
        bandwidth <- bw.nrd0(y_sub)  # Use Silverman's rule of thumb
      }
      dens <- density(y_sub, bw = bandwidth)
      
      # Map Y values to their density
      y_density <- approx(dens$x, dens$y, y_sub)$y
      y_density_scaled <- y_density / max(y_density)
      
      # Generate symmetric jitter
      n_points <- length(idx)
      jitter <- factor * (1 - y_density_scaled) * seq(-1, 1, length.out = n_points)
      
      # Constrain the jitter to a maximum spread
      jitter <- pmin(pmax(jitter, -max_spread), max_spread)
      
      # Randomize order to avoid bias
      jitter <- jitter[order(runif(n_points))]
      
      # Apply jitter to the central X value
      central_x <- x[idx][1]  # Assuming all X values in this category are the same
      jittered_x[idx] <- central_x + jitter
    }
    
    return(jittered_x)
  }
  
  
  X_nets<-
    case_when(
      fd$period=="1.morning"  ~ 1,
      fd$period=="2.afternoon"  ~ 2,
      fd$period=="3.night" ~ 3
    )
  
  
  col_nets<-ifelse(fd$period=="1.morning", dielcol[1],
                   ifelse(fd$period=="2.afternoon", dielcol[2],dielcol[3]))
  
  add_period_segments <- function(df, 
                                  year_col = "year", 
                                  site_col = "site",
                                  period_col = "period",
                                  z_col = "z_score",
                                  col = alpha("lightgray",0.3),
                                  lwd = 1,
                                  lty = 2) {
    
    # Extract numeric period index (1,2,3)
    df$period_num <- as.numeric(sub("\\..*", "", df[[period_col]]))
    
    # Split by year × site
    groups <- split(df, interaction(df[[year_col]], df[[site_col]]))
    
    # Draw segments for each group
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
  
  add_period_segments(fd)
  points_nets<-19 
  
  set.seed(11)
  jittered_x <- jitter_x_by_density(
    x = X_nets, 
    y = fd$z_score,
    categories = X_nets,  
    factor = 0.1, 
    bandwidth = 0.1,
    max_spread = 0.25
  )
  
  points(
    x = X_nets,#jittered_x,
    y = fd$z_score,
    col = alpha(col_nets, 0.2),
    pch = points_nets,
    cex=log(fd$size)
  )
  
  segments(x0=x_nets,
           x1=x_nets,
           y0=sot$asymp.LCL ,
           y1=sot$asymp.UCL , 
           col="black",
           lwd=1.3
           
  )
  
  points(x=x_nets, 
         y=sot$emmean ,
         pch=22,
         cex=1.2,
         bg=rep(c("black"), 3)
         
  )
}


add_outputs<-function(sot2){
  col_ann<-"gray50"
  segments(x0=1:3,
           x1=c(2,3,1),
           y0=sot2$emmeans$asymp.UCL[c(2,2,2)]+c(0.3,0.5,0.8) ,
           y1=sot2$emmeans$asymp.UCL[c(2,2,2)]+c(0.3,0.5,0.8), 
           col=col_ann,
           lwd=1
  )
  
  segments(x0=c(1,2,2,3,3,1),
           x1=c(1,2,2,3,3,1),
           y0=sot2$emmeans$asymp.UCL[c(2,2,2,2,2,2)]+ c(0.3, 0.3, 0.5,0.5,0.8,0.8),
           y1=sot2$emmeans$asymp.UCL[c(2,2,2,2,2,2)]+ c(0.25, 0.25,0.45,0.45,0.75,0.75),
           col=col_ann)
  
  text(x=c(1.5,2, 2.5), 
       y=sot2$emmeans$asymp.UCL[c(2,2,2)] +c(0.4,0.9,0.6),
       labels=ifelse(sot2$contrasts$p.value< 0.001,"P < 0.001",
                     paste0("P = ",sot2$contrasts$p.value %>% round(.,digits=3))
       ) ,
       cex=0.7,
       col=ifelse(sot2$contrasts$p.value <=0.05, "darkred", "gray10"))
  
  points(x=1:3, 
         y=sot2$emmeans$asymp.UCL +0.1, 
         col="darkred",
         pch= ifelse(
           sot2$emmeans$p.value<=0.05, "*", ""
         ))
  
}


## Exporting the plot
dielcol<-c("gold", "#FEA501", "#35336b")
png("Figures/Figure S1.png", res=600, width=5000, height=2000)
layout(matrix(c(1,2,3), ncol=3, byrow=TRUE), widths=c(45,45,20))

par(mar=c(3,4,2,0))

plot_assortativ(fd, sot$emmeans, "A) Unweighted")
add_outputs(sot)

plot_assortativ(fd2, sot2$emmeans, "B) Weighted" , ylab="")


add_outputs(sot2)
par(mar=c(0,0,0,0))

plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n",
     bty="n", xlab="", ylab="" )
legend(x=0, y=0.7, 
       legend=seq(min(full_data$size),max(full_data$size), length=5)%>% round,
       pch=19,
       bty="n",
       pt.cex=seq(min(full_data$size),max(full_data$size), length=5)%>% round %>% log , 
       col=alpha("gray", 0.4),
       title="Plant richness",
       title.adj=4)

dev.off()

