
#------------------------------------------------------------------------------#
#                                                                              #
#                               1. Loading data                                #
#                                                                              #
#------------------------------------------------------------------------------#

packs<-c("tidyverse")

for (i in 1:length(packs)){
  
  if(!require(packs[i],character.only=TRUE)){
    install.packages(packs[i])
  }
  library(packs[i],character.only=TRUE)
}


#------------------------------------------------------------------------------#
#                                                                              #
#                        2. Preparing data for plotting                        #
#                                                                              #
#------------------------------------------------------------------------------# 


smdf<-readRDS("Data/Processed objects/Fig_2.rds")


all_nets<-data.frame(source=unique(smdf$source))

fd_night<-
  smdf[smdf$record_period=="3.night", c("source", "antagonistic", "mutualistic")]

fd_morning<-
  smdf[smdf$record_period=="1.morning", c("source",  "antagonistic", "mutualistic")]

fd_afternoon<-
  smdf[smdf$record_period=="2.afternoon", c("source",  "antagonistic", "mutualistic")]

all_nets %>% left_join(fd_night, by="source") %>%
  replace_na(list(antagonistic=0, mutualistic=0)) -> fd_night



all_nets %>% left_join(fd_morning, by="source") %>%
  replace_na(list(antagonistic=0, mutualistic=0)) -> fd_morning



all_nets %>% left_join(fd_afternoon, by="source") %>%
  replace_na(list(antagonistic=0, mutualistic=0)) -> fd_afternoon

all_pie<-list(morning=fd_morning, afternoon=fd_afternoon, night=fd_night)
all_pie

eee<-all_pie$night[,c("antagonistic", "mutualistic")]
eee<-rowSums(eee)
eee[which(eee==0)]<-NA
eee %>% mean(na.rm=TRUE)
eee %>% sd(na.rm=TRUE)

presence_tbl <- tibble( #order the dataset for plotting
  orig_row = seq_len(nrow(all_pie$night)),
  morning_total   = all_pie$morning$antagonistic   + all_pie$morning$mutualistic,
  afternoon_total = all_pie$afternoon$antagonistic + all_pie$afternoon$mutualistic,
  night_total     = all_pie$night$antagonistic     + all_pie$night$mutualistic
) %>%
  mutate(
    # Count how many periods have >0 observations
    n_periods = (morning_total > 0) +
      (afternoon_total > 0) +
      (night_total > 0)
  ) %>%
  arrange(
    desc(n_periods),
    desc(afternoon_total)
  )

# Extract the new ordering
order_idx <- presence_tbl$orig_row

# Apply ordering to all three data frames
all_pie$morning   <- all_pie$morning[order_idx, ]
all_pie$afternoon <- all_pie$afternoon[order_idx, ]
all_pie$night     <- all_pie$night[order_idx, ]


#------------------------------------------------------------------------------#
#                                                                              #
#                               3. Plotting                                    #
#                                                                              #
#------------------------------------------------------------------------------#

png("Figures/Figure 2.png", width=4000, height=4000,  res=600)
layout(matrix(c(4,4,4,1,2,3), ncol=3, byrow=TRUE), heights=c(1,10))
par( mar=c(4,6,3,0))
### morning
plot(y=c(1:nrow(all_pie$afternoon), 1:nrow(all_pie$afternoon)),
     x=c(all_pie$morning$antagonistic,
         all_pie$morning$antagonistic+
           all_pie$morning$mutualistic),
     xlim=c(0,115),
     xlab = "",
     ylab="Communities",
     xaxt="n",
     yaxt="n",
     bty="n",
     main="A) Morning",
     col="white")

abline(v=0, lty=2, col="gray50")
rect(0, 1:nrow(all_pie$afternoon)-0.4, 
     all_pie$morning$antagonistic+
       all_pie$morning$mutualistic, 1:nrow(all_pie$afternoon)+0.4, col="#D9D9D9", border=NA)
rect(all_pie$morning$antagonistic, 1:nrow(all_pie$afternoon)-0.4,
     0, 1:nrow(all_pie$afternoon)+0.4, col="#7F7F7F", border=NA) 

poin<-c(
  all_pie$morning$antagonistic+
    all_pie$morning$mutualistic)

points(y=c( 1:nrow(all_pie$afternoon)),
       x=poin,
       pch=ifelse(poin>0, 23 ,NA),
       bg="gray20",
       col="black")


axis(1, at=seq(0,max(all_pie$afternoon$antagonistic+
                       all_pie$afternoon$mutualistic), 5 ))
axis(2, at=1:nrow(all_pie$afternoon), labels=all_pie$morning$source, las=1, cex.axis=0.7)


###afternoon
plot(y=c(1:nrow(all_pie$afternoon), 1:nrow(all_pie$afternoon)),
     x=c(all_pie$afternoon$antagonistic,
         all_pie$afternoon$antagonistic+
           all_pie$afternoon$mutualistic),
     xlim=c(0,115),
     xlab = "Visitor richness",
     ylab="",
     xaxt="n",
     yaxt="n",
     bty="n",
     main="B) Afternoon",
     col="white")

abline(v=0, lty=2, col="gray50")

rect(0, 1:nrow(all_pie$afternoon)-0.4, 
     all_pie$afternoon$antagonistic+
       all_pie$afternoon$mutualistic, 1:nrow(all_pie$afternoon)+0.4, col="#D9D9D9", border=NA)
rect(all_pie$afternoon$antagonistic, 1:nrow(all_pie$afternoon)-0.4,
     0, 1:nrow(all_pie$afternoon)+0.4, col="#7F7F7F", border=NA)

poin<-c(
  all_pie$afternoon$antagonistic+
    all_pie$afternoon$mutualistic)

points(y=c( 1:nrow(all_pie$afternoon)),
       x=poin,
       pch=ifelse(poin>0, 23 ,NA),
       bg="gray20",
       col="black")

axis(1, at=seq(0,max(all_pie$afternoon$antagonistic+
                       all_pie$afternoon$mutualistic), 5 ))

###night
plot(y=c(1:nrow(all_pie$afternoon), 1:nrow(all_pie$afternoon)),
     x=c(all_pie$night$antagonistic, 
         all_pie$night$antagonistic+
           all_pie$night$mutualistic),
     xlim=c(0,115),
     xlab = "",
     ylab="",
     xaxt="n",
     yaxt="n",
     bty="n",
     main="C) Night",
     col="white")

abline(v=0, lty=2, col="gray50")

rect(0, 1:nrow(all_pie$afternoon)-0.4, 
     all_pie$night$antagonistic+
       all_pie$night$mutualistic, 1:nrow(all_pie$afternoon)+0.4, col="#D9D9D9", border=NA)
rect(all_pie$night$antagonistic, 1:nrow(all_pie$afternoon)-0.4,
     0, 1:nrow(all_pie$afternoon)+0.4, col="#7F7F7F", border=NA)

poin<-c(
  all_pie$night$antagonistic+
    all_pie$night$mutualistic)
points(y=c( 1:nrow(all_pie$afternoon)),
       x=poin,
       pch=ifelse(poin>0, 23 ,NA),
       bg="gray20",
       col="black")


axis(1, at=seq(0,max(all_pie$afternoon$antagonistic+
                       all_pie$afternoon$mutualistic), 5 ))

par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL,xlim=c(0,1), ylim=c(0,1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
legend(x=0.2, y=0.5,
       legend=c("Herbivores (Antagonists)", "Pollinators (Mutualists)"),
       fill=c("#7F7F7F", "#D9D9D9"),
       bty="n",
       cex=1.2,
       ncol=2)

dev.off()
