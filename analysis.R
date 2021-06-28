library(rethinking)
library(tidyverse)
library(loo)
library(ggridges)
library(viridis)

d <- read.csv("data.csv")

# cortisol needs to be log transformed
# bmi at time of hair cort
# AA = african american (hair type)
# weeks = weeks of preg
# PregRRR = 10 stressful events during pregnancy
# aces10tot = first ten adverse childhood events
# pcl5tot = prenatal ptsd
# xpcl5tot = postnatal ptsd
# language = 1 = english, 2 = spanish
# educ = years edu
# minc = monthly income
# mincnum = num of ppl supported

## Standardizing variables ##
log_cort <- log(d$cortisol)
log_cort_s <- (log_cort - median(log_cort, na.rm=T)) / (sd(log_cort, na.rm=T)*2)

## Indexing missing cort
log_cort_missing <- ifelse(is.na(log_cort_s), 1, 0)
log_cort_missing <- sapply( 1:length(log_cort_missing), 
                            function(n) log_cort_missing[n]*sum(log_cort_missing[1:n]))
log_cort_s <- ifelse(is.na(log_cort_s), -99, log_cort_s)

weeks_s <- (d$weeks - median(d$weeks, na.rm=T)) / (sd(d$weeks, na.rm=T)*2)
age_s <- (d$age - median(d$age, na.rm=T)) / (sd(d$age, na.rm=T)*2)
weeks_s <- (d$weeks - median(d$weeks, na.rm=T)) / (sd(d$weeks, na.rm=T)*2)
spanish <- ifelse(d$language == 1, 1, 0)
educ_s <- (d$educ - median(d$educ, na.rm=T)) / (sd(d$educ, na.rm=T)*2)
secondtrim<-(d$secondtrim)



### Making histogram of outcome ####
par(cex=1.4)
hist(d$pcl5tot, breaks=10,  lty="blank", col=col.alpha(viridis(n=12, option="C")[6], 0.8), main="", xlab="PCL5 Score", xlim=c(0,80))
abline(v=33, lwd=3, lty=2, col="darkred" )
text(x=19.5, y=17, labels="Clinical Cutoff", col="darkred", cex=0.8)
arrows(x0=35, x1=max(d$pcl5tot), y0=16.9, y1=16.9, length=0.1,angle=30,lwd=3, col="darkred")


### Making histogram of outcome ####
par(cex=1.4)
hist(d$pcl5tot, breaks=10,  lty="blank", col=col.alpha(viridis(n=12, option="C")[6], 0.8), main="", xlab="PCL5 Score", xlim=c(0,80))
abline(v=33, lwd=3, lty=2, col="darkred" )
text(x=19.5, y=17, labels="Clinical Cutoff", col="darkred", cex=0.8)
arrows(x0=35, x1=max(d$pcl5tot), y0=16.9, y1=16.9, length=0.1,angle=30,lwd=3, col="darkred")


### Organizing data
data_list <- data.frame(
  pcl5 = d$pcl5tot,
  pcl5_2 = d$pcl5tot,
  log_cort_s = log_cort_s,
  log_cort = log_cort,
  log_cort_missing = log_cort_missing,
  weeks_s = weeks_s,
  age_s = age_s,
  stress = d$PregRRR,
  aces = d$aces10tot,
  AA = d$AA,
  spanish = spanish,
  edu_s = educ_s,
  eth = d$eth
)

data_list$AA <- ifelse(is.na(data_list$AA), -99, data_list$AA)
data_list$stress <- ifelse(is.na(data_list$stress), -99, data_list$stress)
data_list$aces <- ifelse(is.na(data_list$aces), -99, data_list$aces)

data_list <- as.list(data_list)
data_list$N <- nrow(d)
data_list$N_eth <- max(d$eth)
data_list$max_stress <- 10
data_list$max_aces <- 10
data_list$max_pcl5 <- 80
data_list$N_miss_cort <- max(log_cort_missing)

####### Fitting models ####################################################
fit_m1 <- stan( file="stan-models/m1_pre.stan" , data=data_list, chains=3, cores=3, iter=1500, init=0.1, control=list(adapt_delta=0.9) )
fit_m2 <- stan( file="stan-models/m2_pre.stan" , data=data_list, chains=3, cores=3, iter=1500, init=0.1, control=list(adapt_delta=0.9) )

post1 <- extract.samples(fit_m1)
post2 <- extract.samples(fit_m2)

####### LOO CV ###########################################################
loo1 <- loo::loo(fit_m1, save_psis = TRUE, cores = 2, k_threshold = 0.7)
loo2 <- loo::loo(fit_m2, save_psis = TRUE, cores = 2, k_threshold = 0.7)

plot(loo1, label_points = TRUE)
plot(loo2, label_points = TRUE)

m1_points <- loo1$pointwise[,4]
m2_points <- loo2$pointwise[,4]

plot( x=seq(from=1, to=98), y=m2_points-m1_points )

# Weights
model_list <- list(fit_m1, fit_m2)
log_lik_list <- lapply(model_list, extract_log_lik)
weights <- loo::loo_model_weights(log_lik_list, method="stacking")

diff <- loo1$looic - loo2$looic
weight1 <- exp(-0.5 * diff) / ( 1 + exp(-0.5 * diff) )

### Working with posterior samples from m2 ####
post <- extract.samples(fit_m2)

n_samps <- length(post$a[,1])

## Model prediction function, for convienience
mo <- function( scale, x ) {
  if (x == 0)
    return(0)
  
   else 
    return(sum(scale[1:x]))
}

pred_function <- function( aces, log_cort, stress  ) {
  
  preds <- rep(NA, n_samps)
  
  for (i in 1:n_samps)
  preds[i] <- post$a[i] + post$b[i,1]*mo(post$scale_stress[i,], stress) + post$b[i,2]*log_cort + post$b[i,3]*mo(post$scale_aces[i,], aces) + post$b[i,4]*mo(post$scale_stress[i,], stress)*log_cort + post$b[i,5]*mo(post$scale_stress3[i,], stress)*mo(post$scale_aces2[i,], aces) + post$b[i,6]*mo(post$scale_aces3[i,], aces)*log_cort + post$b[i,7]*mo(post$scale_stress4[i,],stress)*mo(post$scale_aces4[i,],aces)*log_cort
  
  # Expected value
  return( 80*logistic(preds) )
}

# predictive sequence for aces
pred_seq <- seq(from=0, to=10)
med_stress <- median(data_list$stress, na.rm=T)
max_stress <- max(data_list$stress, na.rm=T)

### Generating model predictions #######################################
{
  low_stress_low_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    low_stress_low_cort[,j] = pred_function(pred_seq[j], -1, 0)
  }
  
  low_stress_med_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    low_stress_med_cort[,j] = pred_function(pred_seq[j], 0, 0)
  }
  
  low_stress_high_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    low_stress_high_cort[,j] = pred_function(pred_seq[j], 1, 0)
  }
  
  med_stress_low_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    med_stress_low_cort[,j] = pred_function(pred_seq[j], -1, med_stress)
  }
  
  med_stress_med_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    med_stress_med_cort[,j] = pred_function(pred_seq[j], 0, med_stress)
  }
  
  med_stress_high_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    med_stress_high_cort[,j] = pred_function(pred_seq[j], 1, med_stress)
  }
  
  high_stress_low_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    high_stress_low_cort[,j] = pred_function(pred_seq[j], -1, max_stress)
  }
  
  high_stress_med_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    high_stress_med_cort[,j] = pred_function(pred_seq[j], 0, max_stress)
  }
  
  high_stress_high_cort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
  for (j in 1:length(pred_seq)) {
    high_stress_high_cort[,j] = pred_function(pred_seq[j], 1, max_stress)
  }
} ## end prediction loop

{
par(mfrow=c(1,3), cex=1.4, mar=c(5,0,2,1.5), oma=c(0,4,0,1))

plot(NULL, xlim=c(0,10), ylim=c(0,80), xlab="", ylab="", xaxt='n', xaxs='i')
mtext("Prenatal PCL5", side=2, line=2.5, cex=1.4)
axis(1, at=c(0, 5, 10))
## Low ACES
shade(apply(low_stress_low_cort, 2, PI, prob=0.95), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(low_stress_low_cort, 2, PI, prob=0.76), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(low_stress_low_cort, 2, PI, prob=0.57), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(low_stress_low_cort, 2, PI, prob=0.38), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(low_stress_low_cort, 2, PI, prob=0.19), col=col.alpha("skyblue", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(low_stress_low_cort, 2, PI, prob=0.95)[1,], col=col.alpha("skyblue",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(low_stress_low_cort, 2, PI, prob=0.95)[2,], col=col.alpha("skyblue",0.9), lty="dotted")

shade(apply(low_stress_med_cort, 2, PI, prob=0.95), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(low_stress_med_cort, 2, PI, prob=0.76), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(low_stress_med_cort, 2, PI, prob=0.57), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(low_stress_med_cort, 2, PI, prob=0.38), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(low_stress_med_cort, 2, PI, prob=0.19), col=col.alpha("violet", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(low_stress_med_cort, 2, PI, prob=0.95)[1,], col=col.alpha("violet",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(low_stress_med_cort, 2, PI, prob=0.95)[2,], col=col.alpha("violet",0.9), lty="dotted")

shade(apply(low_stress_high_cort, 2, PI, prob=0.95), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(low_stress_high_cort, 2, PI, prob=0.76), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(low_stress_high_cort, 2, PI, prob=0.57), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(low_stress_high_cort, 2, PI, prob=0.38), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(low_stress_high_cort, 2, PI, prob=0.19), col=col.alpha("orange", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(low_stress_high_cort, 2, PI, prob=0.95)[1,], col=col.alpha("orange",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(low_stress_high_cort, 2, PI, prob=0.95)[2,], col=col.alpha("orange",0.9), lty="dotted")

lines(x=pred_seq, y=apply(low_stress_med_cort, 2, median), lwd=2, col="violet")
lines(x=pred_seq, y=apply(low_stress_low_cort, 2, median), lwd=2, col="skyblue")
lines(x=pred_seq, y=apply(low_stress_high_cort, 2, median), lwd=2, col="orange")

text(x=2.5, y=37, labels="Clinical Cutoff", col="darkred", cex=0.7)
abline(h=33, lwd=2, lty="dashed", col="darkred")
mtext("Low Prenatal Stressors",  cex=1.4, line=0.5)

## Med ACES
par(mar=c(5,0,2,1.5))

plot(NULL, xlim=c(0,10), ylim=c(0,80), xlab="", ylab="", xaxt='n', yaxt='n', xaxs='i')
axis(1, at=c(0, 5, 10))
mtext(side=1, text="Early Life Adversity", line=2.5, cex=1.4)
shade(apply(med_stress_low_cort, 2, PI, prob=0.95), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(med_stress_low_cort, 2, PI, prob=0.76), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(med_stress_low_cort, 2, PI, prob=0.57), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(med_stress_low_cort, 2, PI, prob=0.38), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(med_stress_low_cort, 2, PI, prob=0.19), col=col.alpha("skyblue", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(med_stress_low_cort, 2, PI, prob=0.95)[1,], col=col.alpha("skyblue",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(med_stress_low_cort, 2, PI, prob=0.95)[2,], col=col.alpha("skyblue",0.9), lty="dotted")

shade(apply(med_stress_med_cort, 2, PI, prob=0.95), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(med_stress_med_cort, 2, PI, prob=0.76), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(med_stress_med_cort, 2, PI, prob=0.57), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(med_stress_med_cort, 2, PI, prob=0.38), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(med_stress_med_cort, 2, PI, prob=0.19), col=col.alpha("violet", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(med_stress_med_cort, 2, PI, prob=0.95)[1,], col=col.alpha("violet",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(med_stress_med_cort, 2, PI, prob=0.95)[2,], col=col.alpha("violet",0.9), lty="dotted")

shade(apply(med_stress_high_cort, 2, PI, prob=0.95), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(med_stress_high_cort, 2, PI, prob=0.76), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(med_stress_high_cort, 2, PI, prob=0.57), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(med_stress_high_cort, 2, PI, prob=0.38), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(med_stress_high_cort, 2, PI, prob=0.19), col=col.alpha("orange", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(med_stress_high_cort, 2, PI, prob=0.95)[1,], col=col.alpha("orange",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(med_stress_high_cort, 2, PI, prob=0.95)[2,], col=col.alpha("orange",0.9), lty="dotted")

lines(x=pred_seq, y=apply(med_stress_med_cort, 2, median), lwd=2, col="violet")
lines(x=pred_seq, y=apply(med_stress_low_cort, 2, median), lwd=2, col="skyblue")
lines(x=pred_seq, y=apply(med_stress_high_cort, 2, median), lwd=2, col="orange")

legend(x=0.25, y=80, legend=c("Low Cort", "Med Cort", "High Cort"), col=c("skyblue", "violet", "orange"), lty=1, cex=0.65, lwd=2, bty="n")
abline(h=33, lwd=2, lty="dashed", col="darkred")
mtext("Med Prenatal Stressors",  cex=1.4, line=0.5)

## High ACES
par(mar=c(5,0,2,1.5))
plot(NULL, xlim=c(0,10), ylim=c(0,80), xlab="", ylab="", xaxt='n', yaxt='n', xaxs='i')
axis(1, at=c(0, 5, 10))
shade(apply(high_stress_low_cort, 2, PI, prob=0.95), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(high_stress_low_cort, 2, PI, prob=0.76), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(high_stress_low_cort, 2, PI, prob=0.57), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(high_stress_low_cort, 2, PI, prob=0.38), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(high_stress_low_cort, 2, PI, prob=0.19), col=col.alpha("skyblue", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(high_stress_low_cort, 2, PI, prob=0.95)[1,], col=col.alpha("skyblue",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(high_stress_low_cort, 2, PI, prob=0.95)[2,], col=col.alpha("skyblue",0.9), lty="dotted")

shade(apply(high_stress_med_cort, 2, PI, prob=0.95), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(high_stress_med_cort, 2, PI, prob=0.76), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(high_stress_med_cort, 2, PI, prob=0.57), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(high_stress_med_cort, 2, PI, prob=0.38), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(high_stress_med_cort, 2, PI, prob=0.19), col=col.alpha("violet", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(high_stress_med_cort, 2, PI, prob=0.95)[1,], col=col.alpha("violet",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(high_stress_med_cort, 2, PI, prob=0.95)[2,], col=col.alpha("violet",0.9), lty="dotted")

shade(apply(high_stress_high_cort, 2, PI, prob=0.95), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(high_stress_high_cort, 2, PI, prob=0.76), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(high_stress_high_cort, 2, PI, prob=0.57), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(high_stress_high_cort, 2, PI, prob=0.38), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(high_stress_high_cort, 2, PI, prob=0.19), col=col.alpha("orange", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(high_stress_high_cort, 2, PI, prob=0.95)[1,], col=col.alpha("orange",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(high_stress_high_cort, 2, PI, prob=0.95)[2,], col=col.alpha("orange",0.9), lty="dotted")

lines(x=pred_seq, y=apply(high_stress_med_cort, 2, median), lwd=2, col="violet")
lines(x=pred_seq, y=apply(high_stress_low_cort, 2, median), lwd=2, col="skyblue")
lines(x=pred_seq, y=apply(high_stress_high_cort, 2, median), lwd=2, col="orange")

abline(h=33, lwd=2, lty="dashed", col="darkred")
mtext("High Prenatal Stressors",  cex=1.4, line=0.5)
}

dev.off() # clears graphics par

### Interaction of Prenatal Stress x log(Cort)
pred_seq <- seq(from=0,to=7)
stress_lowcort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
stress_medcort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))
stress_highcort <- matrix(NA, nrow=n_samps, ncol=length(pred_seq))

med_aces <- median(data_list$aces, na.rm=T)

for (j in 1:length(pred_seq)) {
  stress_lowcort[,j] <- pred_function(med_aces, -1, pred_seq[j])
  stress_medcort[,j] <- pred_function(med_aces, 0, pred_seq[j])
  stress_highcort[,j] <- pred_function(med_aces, 1, pred_seq[j])
}

par(cex=1.4, mar=c(5,5,3,3))
plot(NULL, xlim=c(0,7), ylim=c(0,80), xlab="Prenatal Stressors", ylab="Prenatal PCL5", xaxs='i')
mtext(side=1, text="", line=2.5, cex=1.4)
shade(apply(stress_lowcort, 2, PI, prob=0.95), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(stress_lowcort, 2, PI, prob=0.76), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(stress_lowcort, 2, PI, prob=0.57), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(stress_lowcort, 2, PI, prob=0.38), col=col.alpha("skyblue", 0.08), pred_seq)
shade(apply(stress_lowcort, 2, PI, prob=0.19), col=col.alpha("skyblue", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(stress_lowcort, 2, PI, prob=0.95)[1,], col=col.alpha("skyblue",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(stress_lowcort, 2, PI, prob=0.95)[2,], col=col.alpha("skyblue",0.9), lty="dotted")

shade(apply(stress_medcort, 2, PI, prob=0.95), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(stress_medcort, 2, PI, prob=0.76), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(stress_medcort, 2, PI, prob=0.57), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(stress_medcort, 2, PI, prob=0.38), col=col.alpha("violet", 0.08), pred_seq)
shade(apply(stress_medcort, 2, PI, prob=0.19), col=col.alpha("violet", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(stress_medcort, 2, PI, prob=0.95)[1,], col=col.alpha("violet",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(stress_medcort, 2, PI, prob=0.95)[2,], col=col.alpha("violet",0.9), lty="dotted")

shade(apply(stress_highcort, 2, PI, prob=0.95), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(stress_highcort, 2, PI, prob=0.76), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(stress_highcort, 2, PI, prob=0.57), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(stress_highcort, 2, PI, prob=0.38), col=col.alpha("orange", 0.08), pred_seq)
shade(apply(stress_highcort, 2, PI, prob=0.19), col=col.alpha("orange", 0.08), pred_seq)

#lines(x=pred_seq, y=apply(stress_highcort, 2, PI, prob=0.95)[1,], col=col.alpha("orange",0.9), lty="dotted")
#lines(x=pred_seq, y=apply(stress_highcort, 2, PI, prob=0.95)[2,], col=col.alpha("orange",0.9), lty="dotted")

lines(x=pred_seq, y=apply(stress_medcort, 2, median), lwd=2, col="violet")
lines(x=pred_seq, y=apply(stress_lowcort, 2, median), lwd=2, col="skyblue")
lines(x=pred_seq, y=apply(stress_highcort, 2, median), lwd=2, col="orange")

#abline(h=33.5, lwd=2, col="red")
legend(x=0.25, y=80, legend=c("Low Cort", "Med Cort", "High Cort"), col=c("skyblue", "violet", "orange"), lty=1, cex=0.7, lwd=2, bty='n')
text(x=1.25, y=36, labels="Clinical Cutoff", col="darkred", cex=0.7)
abline(h=33, lwd=2, lty="dashed", col="darkred")

dev.off()
#### Making density plots #######################

## Estimating observation-level variance
olv_beta_binom <- function( theta ) {
  olv = theta * (pi^2/3)
  return(olv)
}

olv <- olv_beta_binom( post$theta_pcl5 )

forest_df <- data.frame(
  prenatal_stress = post$b[,1],
  log_cortisol = post$b[,2],
  early_life_stress = post$b[,3],
  prenatal_x_cort = post$b[,4],
  prenatal_x_early = post$b[,5],
  early_x_cort = post$b[,6],
  prenatal_x_cort_early = post$b[,7],
  age = post$b[,8],
  spanish = post$b[,9],
  african_american = post$b[,10],
  edu = post$b[,11],
  weeks = post$b[,12]
)

names(forest_df) <- c("Prenatal Stressors", "ln(Cortisol)","Early Life Adversity", "Prenatal Stressors*Cortisol", "Prenatal Stressors*Early Life Adversity", "Early Life Adversity*Cortisol", "Prenatal*Cortisol*Early Life Adversity", "Age", "Spanish Language", "African American", "Years Education", "Weeks of Pregnancy")

# Converting to d
forest_df <- forest_df / sqrt(olv)

# We want to sort by the absolute value of the median d
median_d <- abs( apply(forest_df, 2, median) )

plot_cols <- viridis(n=12, option="C")[6]

cohens_df <- forest_df %>% gather(key="var", value="est")
cohens_df$var <- factor(cohens_df$var, levels=names( sort(median_d) ))

ggplot(cohens_df, aes(x=est, y=var)) + geom_hline(yintercept = c(1:ncol(forest_df)), alpha=0.6) + geom_density_ridges2(aes(),fill=plot_cols, color=plot_cols, alpha=0.5, scale=0.8, rel_min_height=0.01) + ylab("") + xlab("Cohen's d") + theme_bw(base_size=20) + scale_y_discrete(expand = c(0.00, 0)) + theme(legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", axis.text.y=element_text(vjust=0), axis.ticks.y = element_blank()) + geom_vline(xintercept = 0, linetype="dashed", lwd=1) + annotate("blank", x = 0, y=13) + scale_x_continuous(limits=c(-1.2,1.2))

###############################################################