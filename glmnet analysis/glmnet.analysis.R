rm(list=ls())

library(glmnet)
library(dplyr)
library(randomForest)

description <- "Aligned.to.targets.w.primers"
samps.2.exclude <- c("z0132662","z0190864")
min.cov <- 100

load("LASSO.regression.input.and.results.Rdata")
site.names <- names(corrected.pct.meth)
corrected.pct.meth <- corrected.pct.meth[-which(rownames(corrected.pct.meth) %in% samps.2.exclude),]
corrected.pct.meth <- cbind(id=rownames(corrected.pct.meth),corrected.pct.meth)
names(age.data)[which(names(age.data) %in% c("swfsc.labid","age.best","age.confidence"))] <- c("id","age.point","confidence")
age.data$numeric.sex <- 0
age.data$numeric.sex[which(age.data$sex=="Female")] <- 1
age.data$numeric.sex[which(age.data$sex=="Male")] <- 0
age.data$decade <- floor(age.data$age.point/10)
age.data$decade[which(age.data$decade>3)] <- 3
first.meth.col <- dim(age.data)[2] + 1

all.samples <- right_join(age.data,corrected.pct.meth)
calibration.set <- subset(all.samples, subset=(all.samples$confidence >= 3))

################################################################################
# Pearson rank correlations

sex.corr.coeff <- do.call('rbind',lapply(first.meth.col:length(all.samples), function(site){
  temp <- cor.test(all.samples[[site]],all.samples$numeric.sex,method="pearson", use="na.or.complete")
  return(c(corr.coeff=temp$estimate, p.val = temp$p.value))
}))
rownames(sex.corr.coeff) <- site.names

age.corr.coeff <- do.call('rbind',lapply(first.meth.col:length(calibration.set), function(site){
  temp <- cor.test(calibration.set[[site]],calibration.set$age.point,method="pearson", use="na.or.complete")
  return(c(corr.coeff=temp$estimate, p.val = temp$p.value))
}))
rownames(age.corr.coeff) <- site.names

################################################################################
# glmnet analyses
################################################################################

alpha.vals <- seq(from = 0.1, to = 0.9, by = 0.1)

#alternately remove sites and samples to eliminate all NAs

NAs.by.site <- sapply(8:214, function(s){
  length(which(is.na(calibration.set[,s])))
})

calibration.set.complete <- calibration.set[,-(which(NAs.by.site>10)+7)]

NAs.by.sample <- sapply(1:nrow(calibration.set.complete), function(s){
  length(which(is.na(calibration.set.complete[s,])))
})

calibration.set.complete <- calibration.set.complete[-which(NAs.by.sample>6),]

NAs.by.site <- sapply(8:ncol(calibration.set.complete), function(s){
  length(which(is.na(calibration.set.complete[,s])))
})

calibration.set.complete <- calibration.set.complete[,-(which(NAs.by.site>0)+7)]

#################################################
# Cross-validated binomial sex fit over a range of alphas

x.meth <- as.matrix(calibration.set.complete[,first.meth.col:dim(calibration.set.complete)[2]])
y.sex <- as.matrix(calibration.set.complete$numeric.sex)

sex.model <- cv.glmnet(x.meth,y.sex, alpha = 1, family = "binomial", type.measure = "class")
predicted.sex <- predict(sex.model, x.meth, type="class", s = "lambda.min")
sex.errors <- y.sex - as.numeric(predicted.sex)



#################################################
# Cross-validated linear age fit of age over a range of alphas

y.age <- as.matrix(calibration.set.complete$age.point)

test.alpha <- lapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth,y.age, alpha = a)
  corr.coef <- coef(cvfit, s = "lambda.min")
  predicted.age <- predict(cvfit, x.meth, s = "lambda.min")
  age.error <- y.age - predicted.age
  return(list(corr.coef=corr.coef,age.error=age.error, predicted.age=predicted.age,
              mean.age.error=mean(abs(age.error)), median.age.error = median(abs(age.error))))
})

age.error.sum <- t(sapply(test.alpha, function(a){
  c(a$mean.age.error, a$median.age.error)
}))

age.sum <- cbind(calibration.set.complete[,1:13], test.alpha[[9]]$predicted.age)
names(age.sum)[14] <- "predicted.age"
plot(age.sum$age.point, age.sum$predicted.age, pch=(age.sum$confidence+20))
lines(c(0,50),c(0,50))

#################################################
# Cross-validated multinomial fit of decade over a range of alphas

y.decade <- as.matrix(calibration.set.complete$decade)

test.alpha.decade <- lapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth,y.decade, alpha = a, family="multinomial")
  corr.coef <- coef(cvfit, s = "lambda.min")
  predicted.decade <- predict(cvfit, x.meth, s = "lambda.min", type="class")
  decade.error <- y.decade - as.numeric(predicted.decade)
  return(list(corr.coef=corr.coef,decade.error=decade.error, 
              mean.decade.error=mean(abs(decade.error)), median.decade.error = median(abs(decade.error))))
})

decade.error.sum <- t(sapply(test.alpha.decade, function(a){
  c(a$mean.decade.error, a$median.decade.error)
}))

