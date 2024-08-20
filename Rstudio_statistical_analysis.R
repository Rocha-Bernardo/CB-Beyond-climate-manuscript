#libraries
library(data.table)
library(Hmisc)
library(broom)
library(tidyr)
library(corrplot)
library(scales)
library(dplyr)
library(lavaan)
library(semTools)
library(nonnest2)
library(sp)
library(caret)
library(piecewiseSEM)
library(fitdistrplus)
library(writexl)
library(ggplot2)
library(ggpubr)
library(splines)
library(betareg)
#


##Calculate metrics####
#
#libraries.
library(data.table)

#
abs_quad<-read.csv("BD_frequency.csv", header=T, sep=";", row.names = 1)

abs_site<-original%>%
   group_by(SITE)%>%
   summarise_all("mean")%>%
   select(-QUADRAT)

write.csv(abs_site, "abs_site.csv", row.names = FALSE)

#
##Correlations####

#upload files.
BD_Geral<-read.csv("BD_Geral.csv", row.names = 1, header=T, sep=";",)

#calculate correlations values
corstarsl <- function(x){ 
   require(Hmisc) 
   x <- as.matrix(x) 
   R <- rcorr(x,type=c("spearman"))$r 
   p <- rcorr(x,type=c("spearman"))$P 
   
   #define notions for significance levels; spacing is important.
   mystars <- ifelse(p < .001,"***", ifelse(p < .01,"**", ifelse(p < .05,"*","")))
   
   #truncate the matrix that holds the correlations to two decimal
   R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
   
   #build a new matrix that includes the correlations with their stars 
   Rnew <- matrix(paste(R,mystars, sep=""), ncol=ncol(x)) 
   diag(Rnew) <- paste(diag(R), " ", sep="") 
   rownames(Rnew) <- colnames(x) 
   colnames(Rnew) <- paste(colnames(x), " ", sep="") 
   
   #remove upper triangle
   Rnew <- as.matrix(Rnew)
   Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
   Rnew <- as.data.frame(Rnew) 
   
   #remove last column and return the matrix (which is now a data frame)
   Rnew <- cbind(Rnew[1:length(Rnew)-1])
   return(Rnew) 
}

Final_corr <- corstarsl(BD_Geral)
write.csv(Final_corr, "Correlations_spearman.csv")

#
##Linear models ####
#LMs with and without transformation

#without
BD_biocrostas<-read.csv("BD_Geral.csv", sep = ";", row.names = 1)

#Log transform of all variables
BD_biocrostas<- BD_biocrostas+20
log_biocrostas<-log(BD_biocrostas)

#Log transform of response variables
log_biocrostas<-log(BD_biocrostas[,1:3])
log_biocrostas[log_biocrostas == "-Inf"] <- 0

BD_biocrostas[,1:3]<-log_biocrostas

#Square root transform of response variables
sqrt_biocrostas<-sqrt(BD_biocrostas[,1:3])

BD_biocrostas[,1:3]<-sqrt_biocrostas

#Square root transform of response variables and Log of predictors
sqrt_biocrostas<-sqrt(BD_biocrostas[,1:3])

BD_biocrostas[,1:3]<-sqrt_biocrostas

BD_biocrostas[,4:62]<- BD_biocrostas[,4:62]+20
log_biocrostas<-log(BD_biocrostas[,4:62])

BD_biocrostas[,4:62]<- log_biocrostas


#Run Models

BD_BRIO<- log_biocrostas[,c(1,39,40,45,55)]
BD_LICHEN<- log_biocrostas[,c(2,28,39,45,52)]
BD_PLANT<- log_biocrostas[,c(3,28,39,43,52)]

#Loop for all combinations
vars <- colnames(BD_BRIO[,2:5]) #change data here

fit_model <- function(vars) {
   formula <- as.formula(paste("BRIO ~", paste(vars, collapse = "+"))) #change taxa here
   model <- lm(formula, data = BD_BRIO) #change data here
   return(model)
}

combinations <- lapply(1:length(vars), 
                       function(i) combn(vars, i, simplify = FALSE))

A <- lapply(combinations, function(combo) lapply(combo, fit_model))
A <- unlist(A, recursive = FALSE)

#Compare models performance (R2, AdjR2, P, AIC)
summary(A[[8]])
AIC(A[[8]])

#Calculate variance partitioning of all predictors
af <- anova(A[[8]])
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

#
##Structural equation modelling (lavaan) (old version) ####
#https://lavaan.ugent.be/tutorial/index.html

#upload files.
BD_biocrostas<-read.csv("BD_geral.csv", sep = ";", row.names = 1)
str(BD_biocrostas)

#prior to modelling, variables were normalized (Min-Max normalization)
MinMax_normalization <- function(x) {
   (x - min(x)) / (max(x) - min(x))
}

BD_biocrostas[,4:6] <- as.data.frame(lapply(BD_biocrostas[4:6], 
                                            MinMax_normalization))

#model sintax.
SEM1 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2
LICHEN ~ d*Altitude

PLANT ~~ 0*BRIO 
PLANT ~~ 0*LICHEN 
BRIO ~~ 0*LICHEN 

#by manually fixing each value of zero (thus indicating you are not interested 
#in estimating them/are comfortable assuming they take on a value of 0)

#Removing the covariance of End1 and End2 is the same as fixing its covariance 
#to zero through the following syntax End1 ~~ 0*End2
"

SEM2 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2 + d*PLANT
LICHEN ~ e*Altitude + f*PLANT

BRIO ~~ 0*LICHEN 

#indirect effect of Altitude on BRIO via PLANT
BIE1 := a*d
#indirect effect of BIO3 on BRIO via PLANT
BIE2 := b*d
#total effect of Altitude on BRIO
BTE1 := BIE1
#total effect of BIO3 on BRIO
BTE2 := BIE2

#indirect effect of Altitude on LICHEN via PLANT
LIE1 := a*f
#indirect effect of BIO3 on LICHEN via PLANT
LIE2 := b*f
#total effect of Altitude on LICHEN
LTE1 := LIE1+e
#total effect of BIO3 on LICHEN
LTE2 := LIE2
"

SEM3 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2 + d*PLANT
LICHEN ~ e*Altitude

BRIO ~~ 0*LICHEN 

#indirect effect of Altitude on BRIO via PLANT
BIE1 := a*d
#indirect effect of BIO3 on BRIO via PLANT
BIE2 := b*d
#total effect of Altitude on BRIO
BTE1 := BIE1
#total effect of BIO3 on BRIO
BTE2 := BIE2
"

SEM4 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2
LICHEN ~ e*Altitude + f*PLANT

BRIO ~~ 0*LICHEN 

#indirect effect of Altitude on LICHEN via PLANT
LIE1 := a*f
#indirect effect of BIO3 on LICHEN via PLANT
LIE2 := b*f
#total effect of Altitude on LICHEN
LTE1 := LIE1+e
#total effect of BIO3 on LICHEN
LTE2 := LIE2
"

SEM5 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*BRIO
BRIO ~ c*BIO2
LICHEN ~ e*Altitude + f*BRIO

PLANT ~~ 0*LICHEN 

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*d
#total effect of BIO2 on PLANT
PTE1 := PIE1

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*f
#total effect of BIO2 on LICHEN
LTE1 := LIE1
"

SEM6 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*BRIO
BRIO ~ c*BIO2
LICHEN ~ e*Altitude

PLANT ~~ 0*LICHEN 

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*d
#total effect of BIO2 on PLANT
PTE1 := PIE1
"

SEM7 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2
LICHEN ~ e*Altitude + f*BRIO

PLANT ~~ 0*LICHEN 

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*f
#total effect of BIO2 on LICHEN
LTE1 := LIE1
"

SEM8 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*LICHEN
BRIO ~ c*BIO2 + f*LICHEN
LICHEN ~ e*Altitude

PLANT ~~ 0*BRIO 

#indirect effect of Altitude on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude on PLANT
PTE1 := PIE1+a

#indirect effect of Altitude on BRIO via LICHEN
BIE1 := e*f
#total effect of Altitude on BRIO
BTE1 := BIE1
"

SEM9 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*LICHEN
BRIO ~ c*BIO2
LICHEN ~ e*Altitude

PLANT ~~ 0*BRIO 

#indirect effect of Altitude on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude on PLANT
PTE1 := PIE1+a
"

SEM10 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2 + f*LICHEN
LICHEN ~ e*Altitude

PLANT ~~ 0*BRIO 

#indirect effect of Altitude on BRIO via LICHEN
BIE1 := e*f
#total effect of Altitude on BRIO
BTE1 := BIE1
"

SEM11 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2 + d*PLANT
LICHEN ~ e*Altitude + f*PLANT + g*BRIO

#indirect effect of Altitude on BRIO via PLANT
BIE1 := a*d
#indirect effect of BIO3 on BRIO via PLANT
BIE2 := b*d
#total effect of Altitude on BRIO
BTE1 := BIE1
#total effect of BIO3 on BRIO
BTE2 := BIE2

#indirect effect of Altitude on LICHEN via PLANT
LIE1 := a*f
#indirect effect of BIO3 on LICHEN via PLANT
LIE2 := b*f
#total effect of Altitude on LICHEN
LTE1 := LIE1+e
#total effect of BIO3 on LICHEN
LTE2 := LIE2

#indirect effect of BIO2 on LICHEN via BRIO
LIE3 := c*g
#total effect of BIO2 on LICHEN
LTE3 := LIE3
"

SEM12 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2 + d*PLANT + g*LICHEN
LICHEN ~ e*Altitude + f*PLANT 

#indirect effect of Altitude on BRIO via PLANT
BIE1 := a*d
#indirect effect of BIO3 on BRIO via PLANT
BIE2 := b*d
#indirect effect of Altitude on BRIO via LICHEN
BIE3 := e*g
#total effect of Altitude on BRIO
BTE1 := BIE1 + BIE3
#total effect of BIO3 on BRIO
BTE2 := BIE2

#indirect effect of Altitude on LICHEN via PLANT
LIE1 := a*f
#indirect effect of BIO3 on LICHEN via PLANT
LIE2 := b*f
#total effect of Altitude on LICHEN
LTE1 := LIE1+e
#total effect of BIO3 on LICHEN
LTE2 := LIE2
"

SEM13 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2 + d*PLANT + g*LICHEN
LICHEN ~ e*Altitude 

#indirect effect of Altitude on BRIO via PLANT
BIE1 := a*d
#indirect effect of BIO3 on BRIO via PLANT
BIE2 := b*d
#indirect effect of Altitude on BRIO via LICHEN
BIE3 := e*g
#total effect of Altitude on BRIO
BTE1 := BIE1 + BIE3
#total effect of BIO3 on BRIO
BTE2 := BIE2
"

SEM14 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2 + g*LICHEN
LICHEN ~ e*Altitude + f*PLANT 

#indirect effect of Altitude on BRIO via LICHEN
BIE1 := e*g
#total effect of Altitude on BRIO
BTE1 := BIE1


#indirect effect of Altitude on LICHEN via PLANT
LIE1 := a*f
#indirect effect of BIO3 on LICHEN via PLANT
LIE2 := b*f
#total effect of Altitude on LICHEN
LTE1 := LIE1+e
#total effect of BIO3 on LICHEN
LTE2 := LIE2
"

SEM15 <- 
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2 + d*PLANT
LICHEN ~ e*Altitude + g*BRIO

#indirect effect of Altitude on BRIO via PLANT
BIE1 := a*d
#indirect effect of BIO3 on BRIO via PLANT
BIE2 := b*d
#total effect of Altitude on BRIO
BTE1 := BIE1
#total effect of BIO3 on BRIO
BTE2 := BIE2

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*g
#total effect of BIO2 on LICHEN
LTE1 := LIE1
"

SEM16 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*BRIO
BRIO ~ c*BIO2 
LICHEN ~ e*Altitude + f*PLANT + g*BRIO

#indirect effect of BIO2 on PLANT via BRIO
BIE1 := c*d
#total effect of BIO2 on PLANT
BTE1 := BIE1

#indirect effect of Altitude on LICHEN via PLANT
LIE1 := a*f
#indirect effect of BIO3 on LICHEN via PLANT
LIE2 := b*f
#total effect of Altitude on LICHEN
LTE1 := LIE1+e
#total effect of BIO3 on LICHEN
LTE2 := LIE2

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*g
#total effect of BIO2 on LICHEN
LTE1 := LIE1
"

SEM17 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*BRIO + f*LICHEN
BRIO ~ c*BIO2 
LICHEN ~ e*Altitude + g*BRIO

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*d
#indirect effect of Altitude on PLANT via LICHEN
PIE2 := e*f
#total effect of BIO2 on PLANT
PTE1 := PIE1
#total effect of Altitude on PLANT
PTE2 := PIE2


#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*g
#total effect of BIO2 on LICHEN
LTE1 := LIE1
"

SEM18 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*BRIO + f*LICHEN
BRIO ~ c*BIO2 
LICHEN ~ e*Altitude

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*d
#indirect effect of Altitude on PLANT via LICHEN
PIE2 := e*f
#total effect of BIO2 on PLANT
PTE1 := PIE1
#total effect of Altitude on PLANT
PTE2 := PIE2
"

SEM19 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*LICHEN
BRIO ~ c*BIO2 
LICHEN ~ e*Altitude + g*BRIO

#indirect effect of Altitude on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude on PLANT
PTE1 := PIE1

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*g
#total effect of BIO2 on LICHEN
LTE1 := LIE1
"

SEM20 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*BRIO
BRIO ~ c*BIO2 
LICHEN ~ e*Altitude + f*PLANT

#indirect effect of BIO2 on PLANT via BRIO
BIE1 := c*d
#total effect of BIO2 on PLANT
BTE1 := BIE1

#indirect effect of Altitude on LICHEN via PLANT
LIE1 := a*f
#indirect effect of BIO3 on LICHEN via PLANT
LIE2 := b*f
#total effect of Altitude on LICHEN
LTE1 := LIE1+e
#total effect of BIO3 on LICHEN
LTE2 := LIE2
"

SEM21 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*LICHEN
BRIO ~ c*BIO2 + f*PLANT + g*LICHEN
LICHEN ~ e*Altitude

#indirect effect of Altitude on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude on PLANT
PTE1 := PIE1+a

#indirect effect of Altitude on BRIO via PLANT
BIE1 := a*f
#indirect effect of Altitude on BRIO via LICHEN
BIE2 := e*g
#indirect effect of BIO3 on BRIO via PLANT
BIE3 := b*f
#total effect of Altitude on BRIO
BTE1 := BIE1+BIE2
#total effect of BIO3 on BRIO
BTE3 := BIE3
"

SEM22 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*LICHEN + f*BRIO
BRIO ~ c*BIO2 + g*LICHEN
LICHEN ~ e*Altitude

#indirect effect of Altitude on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude on PLANT
PTE1 := PIE1+a

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*f
#total effect of BIO2 on PLANT
PTE1 := PIE1

#indirect effect of Altitude on BRIO via LICHEN
BIE1 := e*g
#total effect of Altitude on BRIO
BTE1 := BIE1
"

SEM23<-
   "
PLANT ~ a*Altitude + b*BIO3
BRIO ~ c*BIO2
LICHEN ~ e*Altitude + f*PLANT + g*BRIO

#indirect effect of Altitude on LICHEN via PLANT
LIE1 := a*f
#indirect effect of BIO3 on LICHEN via PLANT
LIE2 := b*f
#total effect of Altitude on LICHEN
LTE1 := LIE1+e
#total effect of BIO3 on LICHEN
LTE2 := LIE2

#indirect effect of BIO2 on LICHEN via BRIO
LIE3 := c*g
#total effect of BIO2 on LICHEN
LTE3 := LIE3
"

SEM24 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + f*BRIO
BRIO ~ c*BIO2 + g*LICHEN
LICHEN ~ e*Altitude

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*f
#total effect of BIO2 on PLANT
PTE1 := PIE1

#indirect effect of Altitude on BRIO via LICHEN
BIE1 := e*g
#total effect of Altitude on BRIO
BTE1 := BIE1
"

SEM25 <- 
   "
PLANT ~ a*Altitude + b*BIO3 + d*LICHEN
BRIO ~ c*BIO2 + f*PLANT
LICHEN ~ e*Altitude

#indirect effect of Altitude on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude on PLANT
PTE1 := PIE1+a

#indirect effect of Altitude on BRIO via PLANT
BIE1 := a*f
#indirect effect of BIO3 on BRIO via PLANT
BIE3 := b*f
#total effect of Altitude on BRIO
BTE1 := BIE1
#total effect of BIO3 on BRIO
BTE3 := BIE3
"

SEM11 <- sem(SEM11, data = BD_biocrostas, 
             estimator = "MLM", se="robust.sem")   
summary(SEM11, fit.measures=T, standardized=T, rsquare=T)

modindices(SEM11, sort = T, maximum.number = 10)

#
#plot SEM.
semPaths(SEM5, what='std')

#
#compare between models fitness.
anova(SEM11,SEM12,SEM16,SEM17,SEM18,SEM21,SEM22)

#
#smaller difference between matrices, better fit of the SEM to our data
lavInspect(SEM1, what = "sampstat")
lavInspect(SEM1, what = "implied")

#
##Structural equation modelling (lavaan) (new version) ####

BD_biocrostas<- read.csv("BD_geral.csv", sep = ";", row.names = 1)

BD_biocrostas<- BD_biocrostas+20
log_biocrostas<-log(BD_biocrostas)

#model sintax.
SEM1 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ d*Altitude_2m

PLANT ~~ 0*BRIO 
PLANT ~~ 0*LICHEN 
BRIO ~~ 0*LICHEN 

#by manually fixing each value of zero (thus indicating you are not interested 
#in estimating them/are comfortable assuming they take on a value of 0)

#Removing the covariance of End1 and End2 is the same as fixing its covariance 
#to zero through the following syntax End1 ~~ 0*End2
"

SEM2 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + d*PLANT
LICHEN ~ e*Altitude_2m + f*PLANT

BRIO ~~ 0*LICHEN 

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*d
#total effect of Altitude_2m on BRIO
BTE1 := BIE1

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e
"

SEM3 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + d*PLANT
LICHEN ~ e*Altitude_2m

BRIO ~~ 0*LICHEN 

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*d
#total effect of Altitude_2m on BRIO
BTE1 := BIE1
"

SEM4 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m + f*PLANT

BRIO ~~ 0*LICHEN 

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e
"

SEM5 <- 
   "
PLANT ~ a*Altitude_2m + d*BRIO
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m + f*BRIO

PLANT ~~ 0*LICHEN 

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*d
#indirect effect of Slope_2m on PLANT via BRIO
PIE2 := b*d
#total effect of BIO2 on PLANT
PTE1 := PIE1
#total effect of Slope_2m on PLANT
PTE2 := PIE2

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*f
#indirect effect of Slope_2m on LICHEN via BRIO
LIE2 := b*f
#total effect of BIO2 on LICHEN
LTE1 := LIE1
#total effect of Slope_2m on LICHEN
LTE2 := LIE2
"

SEM6 <- 
   "
PLANT ~ a*Altitude_2m + d*BRIO
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m

PLANT ~~ 0*LICHEN 

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*d
#indirect effect of Slope_2m on PLANT via BRIO
PIE2 := b*d
#total effect of BIO2 on PLANT
PTE1 := PIE1
#total effect of Slope_2m on PLANT
PTE2 := PIE2
"

SEM7 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m + f*BRIO

PLANT ~~ 0*LICHEN 

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*f
#indirect effect of Slope_2m on LICHEN via BRIO
LIE2 := b*f
#total effect of BIO2 on LICHEN
LTE1 := LIE1
#total effect of Slope_2m on LICHEN
LTE2 := LIE2
"

SEM8 <- 
   "
PLANT ~ a*Altitude_2m + d*LICHEN
BRIO ~ b*Slope_2m + c*BIO2 + f*LICHEN
LICHEN ~ e*Altitude_2m

PLANT ~~ 0*BRIO 

#indirect effect of Altitude_2m on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude_2m on PLANT
PTE1 := PIE1+a

#indirect effect of Altitude_2m on BRIO via LICHEN
BIE1 := e*f
#total effect of Altitude_2m on BRIO
BTE1 := BIE1
"

SEM9 <- 
   "
PLANT ~ a*Altitude_2m + d*LICHEN
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m

PLANT ~~ 0*BRIO 

#indirect effect of Altitude_2m on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude_2m on PLANT
PTE1 := PIE1+a
"

SEM10 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + f*LICHEN
LICHEN ~ e*Altitude_2m

PLANT ~~ 0*BRIO 

#indirect effect of Altitude_2m on BRIO via LICHEN
BIE1 := e*f
#total effect of Altitude_2m on BRIO
BTE1 := BIE1
"

SEM11 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + d*PLANT
LICHEN ~ e*Altitude_2m + f*PLANT + g*BRIO

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*d
#total effect of Altitude_2m on BRIO
BTE1 := BIE1

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e

#indirect effect of BIO2 on LICHEN via BRIO
LIE2 := c*g
#indirect effect of Slope_2m on LICHEN via BRIO
LIE3 := b*g
#total effect of BIO2 on LICHEN
LTE2 := LIE2
#total effect of Slope_2m on LICHEN
LTE3 := LIE3
"

SEM12 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + d*PLANT + g*LICHEN
LICHEN ~ e*Altitude_2m + f*PLANT 

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*d
#indirect effect of Altitude_2m on BRIO via LICHEN
BIE2 := e*g
#total effect of Altitude_2m on BRIO
BTE1 := BIE1 + BIE2

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e
"

SEM13 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + d*PLANT + g*LICHEN
LICHEN ~ e*Altitude_2m 

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*d
#indirect effect of Altitude_2m on BRIO via LICHEN
BIE2 := e*g
#total effect of Altitude_2m on BRIO
BTE1 := BIE1 + BIE2
"

SEM14 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + g*LICHEN
LICHEN ~ e*Altitude_2m + f*PLANT 

#indirect effect of Altitude_2m on BRIO via LICHEN
BIE1 := e*g
#total effect of Altitude_2m on BRIO
BTE1 := BIE1

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e
"

SEM15 <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + d*PLANT
LICHEN ~ e*Altitude_2m + g*BRIO

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*d
#total effect of Altitude_2m on BRIO
BTE1 := BIE1

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*g
#indirect effect of Slope_2m on LICHEN via BRIO
LIE2 := b*g
#total effect of BIO2 on LICHEN
LTE1 := LIE1
#total effect of Slope_2m on LICHEN
LTE2 := LIE2
"

SEM16 <- 
   "
PLANT ~ a*Altitude_2m + d*BRIO
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m + f*PLANT + g*BRIO

#indirect effect of BIO2 on PLANT via BRIO
BIE1 := c*d
#indirect effect of Slope_2m on PLANT via BRIO
BIE2 := b*d
#total effect of BIO2 on PLANT
BTE1 := BIE1
#total effect of BIO2 on PLANT
BTE2 := BIE2

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e

#indirect effect of BIO2 on LICHEN via BRIO
LIE2 := c*g
#indirect effect of Slope on LICHEN via BRIO
LIE3 := b*g
#total effect of BIO2 on LICHEN
LTE2 := LIE2
#total effect of Slope on LICHEN
LTE3 := LIE3
"

SEM17 <- 
   "
PLANT ~ a*Altitude_2m + d*BRIO + f*LICHEN
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m + g*BRIO

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*d
#indirect effect of Slope_2m on PLANT via BRIO
PIE2 := b*d
#indirect effect of Altitude_2m on PLANT via LICHEN
PIE3 := e*f
#total effect of BIO2 on PLANT
PTE1 := PIE1
#total effect of Slope on PLANT
PTE2 := PIE2
#total effect of Altitude_2m on PLANT
PTE3 := PIE3+a

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*g
#indirect effect of Slope on LICHEN via BRIO
LIE2 := b*g
#total effect of BIO2 on LICHEN
LTE1 := LIE1
#total effect of Slope on LICHEN
LTE2 := LIE2
"

SEM18 <- 
   "
PLANT ~ a*Altitude_2m + d*BRIO + f*LICHEN
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*d
#indirect effect of Slope_2m on PLANT via BRIO
PIE2 := b*d
#indirect effect of Altitude_2m on PLANT via LICHEN
PIE3 := e*f
#total effect of BIO2 on PLANT
PTE1 := PIE1
#total effect of Slope on PLANT
PTE2 := PIE2
#total effect of Altitude_2m on PLANT
PTE3 := PIE3+a
"

SEM19 <- 
   "
PLANT ~ a*Altitude_2m + d*LICHEN
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m + g*BRIO

#indirect effect of Altitude_2m on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude_2m on PLANT
PTE1 := PIE1+a

#indirect effect of BIO2 on LICHEN via BRIO
LIE1 := c*g
#indirect effect of Slope on LICHEN via BRIO
LIE2 := b*g
#total effect of BIO2 on LICHEN
LTE1 := LIE1
#total effect of Slope on LICHEN
LTE2 := LIE2
"

SEM20 <- 
   "
PLANT ~ a*Altitude_2m + d*BRIO
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m + f*PLANT

#indirect effect of BIO2 on PLANT via BRIO
BIE1 := c*d
#indirect effect of Slope on PLANT via BRIO
BIE2 := b*d
#total effect of BIO2 on PLANT
BTE1 := BIE1
#total effect of Slope on PLANT
BTE2 := BIE2

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e
"

SEM21 <- 
   "
PLANT ~ a*Altitude_2m + d*LICHEN
BRIO ~ b*Slope_2m + c*BIO2 + f*PLANT + g*LICHEN
LICHEN ~ e*Altitude_2m

#indirect effect of Altitude_2m on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude_2m on PLANT
PTE1 := PIE1+a

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*f
#indirect effect of Altitude_2m on BRIO via LICHEN
BIE2 := e*g
#total effect of Altitude_2m on BRIO
BTE1 := BIE1+BIE2
"

SEM22 <- 
   "
PLANT ~ a*Altitude_2m + d*LICHEN + f*BRIO
BRIO ~ b*Slope_2m + c*BIO2 + g*LICHEN
LICHEN ~ e*Altitude_2m

#indirect effect of Altitude_2m on PLANT via LICHEN
PIE1 := e*d
#indirect effect of BIO2 on PLANT via BRIO
PIE2 := c*f
#indirect effect of Slope on PLANT via BRIO
PIE3 := b*f
#total effect of Altitude_2m on PLANT
PTE1 := PIE1+a
#total effect of BIO2 on PLANT
PTE2 := PIE2
#total effect of Slope on PLANT
PTE3 := PIE3

#indirect effect of Altitude_2m on BRIO via LICHEN
BIE1 := e*g
#total effect of Altitude_2m on BRIO
BTE1 := BIE1
"

SEM23<-
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2
LICHEN ~ e*Altitude_2m + f*PLANT + g*BRIO

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e

#indirect effect of BIO2 on LICHEN via BRIO
LIE2 := c*g
#indirect effect of Slope on LICHEN via BRIO
LIE3 := b*g
#total effect of BIO2 on LICHEN
LTE2:= LIE2
#total effect of Slope on LICHEN
LTE3 := LIE3
"

SEM24 <- 
   "
PLANT ~ a*Altitude_2m + f*BRIO
BRIO ~ b*Slope_2m + c*BIO2 + g*LICHEN
LICHEN ~ e*Altitude_2m

#indirect effect of BIO2 on PLANT via BRIO
PIE1 := c*f
#indirect effect of Slope on PLANT via BRIO
PIE2 := b*f
#total effect of BIO2 on PLANT
PTE1 := PIE1
#total effect of Slope on PLANT
PTE2 := PIE2

#indirect effect of Altitude_2m on BRIO via LICHEN
BIE1 := e*g
#total effect of Altitude_2m on BRIO
BTE1 := BIE1
"

SEM25 <- 
   "
PLANT ~ a*Altitude_2m + d*LICHEN
BRIO ~ b*Slope_2m + c*BIO2 + f*PLANT
LICHEN ~ e*Altitude_2m

#indirect effect of Altitude_2m on PLANT via LICHEN
PIE1 := e*d
#total effect of Altitude_2m on PLANT
PTE1 := PIE1+a

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*f
#total effect of Altitude_2m on BRIO
BTE1 := BIE1
"

SEM <- sem(SEM1, data = log_biocrostas, 
           estimator = "MLM", se="robust.sem")   
summary(SEM, fit.measures=T, standardized=T, rsquare=T)

modindices(SEM, sort = T, maximum.number = 10)

#
##Calculate net effects from SEM (new version) ####
#(https://jslefche.github.io/sem_book/composite-variables.html)

#upload files
BD_biocrostas<- read.csv("BD_geral.csv", sep = ";", row.names = 1)

BD_biocrostas<- BD_biocrostas+20
log_biocrostas<-log(BD_biocrostas)

#calculate bryophythes net effects
Brios_lm<-lm(BRIO~Slope_2m + BIO2 + PLANT, data=log_biocrostas)

summary(Brios_lm)
summary(Brios_lm)$coefficients[2, 1]
summary(Brios_lm)$coefficients[3, 1]
summary(Brios_lm)$coefficients[4, 1]

Brios_sem <- '
Abiotic <~ -0.83816*Slope_2m + -92.77523*BIO2

BRIO ~ Abiotic+PLANT
'
SEM <- sem(Brios_sem, data = log_biocrostas, 
           estimator = "MLM", se="robust.sem")   
summary(SEM, standardized=T, rsquare=T)

#calculate lichens net effects
Lichens_lm<-lm(LICHEN~Altitude_2m + PLANT + BRIO, data=log_biocrostas)

summary(Lichens_lm)
summary(Lichens_lm)$coefficients[2, 1]
summary(Lichens_lm)$coefficients[3, 1]
summary(Lichens_lm)$coefficients[4, 1]

Lichens_sem <- '
biotic <~ -0.5101518*PLANT + -0.1191359*BRIO

LICHEN ~ biotic+Altitude_2m
'

SEM <- sem(Lichens_sem, data = log_biocrostas, 
           estimator = "MLM", se="robust.sem")   
summary(SEM, standardized=T, rsquare=T)

#
##Mapping vegetation from SEM (lavaan)####

#upload files
load(file="BD_predictVeg.Rdata")
BD_biocrostas<- read.csv("BD_geral.csv", sep = ";", row.names = 1)

#remove -9999 values from slope
BD_predictVeg<- BD_predictVeg %>%filter(Slope_2m!='-9999')

#isolate coordinates and environmental variables
coordinates<-BD_predictVeg[,1:2]
BD_predictVeg<- BD_predictVeg [,c(3:6)]

#Log transform data
log_predictVeg<-log(BD_predictVeg)
log_predictVeg[log_predictVeg == "-Inf"] <- 0

##(Winner SEM)
BD_biocrostas<- BD_biocrostas %>% 
   select(BRIO,LICHEN,PLANT, Altitude_2m, Slope_2m, BIO2, F_BIO2)

log_biocrostas<-log(BD_biocrostas)
log_biocrostas[log_biocrostas == "-Inf"] <- 0

SEM11_Present <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*BIO2 + d*PLANT
LICHEN ~ e*Altitude_2m + f*PLANT + g*BRIO

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*d
#total effect of Altitude_2m on BRIO
BTE1 := BIE1

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e

#indirect effect of P_BIO2 on LICHEN via BRIO
LIE2 := c*g
#indirect effect of Slope_2m on LICHEN via BRIO
LIE3 := b*g
#total effect of P_BIO2 on LICHEN
LTE2 := LIE2
#total effect of Slope_2m on LICHEN
LTE3 := LIE3
"

SEM11_Future <- 
   "
PLANT ~ a*Altitude_2m
BRIO ~ b*Slope_2m + c*F_BIO2 + d*PLANT
LICHEN ~ e*Altitude_2m + f*PLANT + g*BRIO

#indirect effect of Altitude_2m on BRIO via PLANT
BIE1 := a*d
#total effect of Altitude_2m on BRIO
BTE1 := BIE1

#indirect effect of Altitude_2m on LICHEN via PLANT
LIE1 := a*f
#total effect of Altitude_2m on LICHEN
LTE1 := LIE1+e

#indirect effect of F_BIO2 on LICHEN via BRIO
LIE2 := c*g
#indirect effect of Slope_2m on LICHEN via BRIO
LIE3 := b*g
#total effect of F_BIO2 on LICHEN
LTE2 := LIE2
#total effect of Slope_2m on LICHEN
LTE3 := LIE3
"

SEM_F <- sem(SEM11_Future, data = log_biocrostas, 
             estimator = "MLM", se="robust.sem")

#Predict present vegetation values 
#(Intended response variables to be predicted must enter as columns with NAs)
BRIO = NA
LICHEN = NA
PLANT = NA
log_predictVeg$PLANT = PLANT

colnames(log_predictVeg)[1] <- "BIO2"
colnames(log_predictVeg)[3] <- "Altitude_2m"


Log_predict_P <-as.data.frame(lavPredictY(SEM_P, newdata = log_predictVeg,
                                          ynames = c("BRIO", "LICHEN", "PLANT"),
                                          xnames = c("Altitude_2m", "Slope_2m", "BIO2"),
                                          method = "conditional.mean", label = T, assemble = F))

#Predict future vegetation values 
Log_predict_F <-as.data.frame(lavPredictY(SEM_F, newdata = log_predictVeg,
                                          ynames = c("BRIO","LICHEN","PLANT"),
                                          xnames = c("Altitude_2m", "Slope_2m", "F_BIO2"),
                                          method = "conditional.mean", 
                                          label = T, assemble = F))


#undo log transform and remove over 100
Veg_predict_P<- exp(Log_predict_P)
Veg_predict_F<- exp(Log_predict_F)

Veg_predict_P[Veg_predict_P > 100]<- 100
Veg_predict_F[Veg_predict_F > 100]<- 100

#add coordinates to points
Veg_predict_P<-cbind(coordinates,Veg_predict_P)
Veg_predict_F<-cbind(coordinates,Veg_predict_F)

#export
Veg_predict_P1<- Veg_predict_P[1:1030000,]
Veg_predict_P2<- Veg_predict_P[1030001:2060000,]
Veg_predict_P3<- Veg_predict_P[2060001:3090000,]
Veg_predict_P4<- Veg_predict_P[3090001:4120000,]
Veg_predict_P5<- Veg_predict_P[4120001:5144134,]

Veg_predict_F1<- Veg_predict_F[1:1030000,]
Veg_predict_F2<- Veg_predict_F[1030001:2060000,]
Veg_predict_F3<- Veg_predict_F[2060001:3090000,]
Veg_predict_F4<- Veg_predict_F[3090001:4120000,]
Veg_predict_F5<- Veg_predict_F[4120001:5144134,]


write_xlsx(Veg_predict_F5, "Veg_predict_F5.xlsx")

#
##Generalized linear models ####
BD_biocrostas<-read.csv("BD_Geral.csv", sep = ";", row.names = 1)

hist(BD_biocrostas$BRIO)
hist(BD_biocrostas$LICHEN)
hist(BD_biocrostas$PLANT)

descdist(BD_biocrostas$BRIO, discrete = T)
Fit_BRIO<- fitdist(BD_biocrostas$BRIO, "pois")
Fit_BRIO<- fitdist(BD_biocrostas$BRIO, "nbinom") #ganha este
plot(Fit_BRIO)

descdist(BD_biocrostas$LICHEN, discrete = T)
Fit_LICHEN<- fitdist(BD_biocrostas$LICHEN, "pois")
Fit_LICHEN<- fitdist(BD_biocrostas$LICHEN, "nbinom") #ganha este
plot(Fit_LICHEN)

descdist(BD_biocrostas$PLANT, discrete = T)
Fit_PLANT<- fitdist(BD_biocrostas$PLANT, "pois")
Fit_PLANT<- fitdist(BD_biocrostas$PLANT, "nbinom") #ganha este
plot(Fit_PLANT)

#Poisson or Quasipoisson
mean(BD_biocrostas$BRIO)
var(BD_biocrostas$BRIO)

mean(BD_biocrostas$LICHEN)
var(BD_biocrostas$LICHEN)

mean(BD_biocrostas$PLANT)
var(BD_biocrostas$PLANT)

#Quasipoisson is more appropriate as there is overdispersion 
#However Poisson is simpler so we use Poisson

BD_BRIO<- BD_biocrostas[,c(1,39,40,46,56)]
BD_LICHEN<- BD_biocrostas[,c(2,28,39,46,53)]
BD_PLANT<- BD_biocrostas[,c(3,28,39,43,53)]

#Loop for all combinations
vars <- colnames(BD_BRIO[,2:5]) #change data here

fit_model <- function(vars) {
   formula <- as.formula(paste("BRIO ~", paste(vars, collapse = "+"))) #change taxa here
   model <- glm(formula, data = BD_BRIO, family = poisson()) #change data here
   return(model)
}

combinations <- lapply(1:length(vars), 
                       function(i) combn(vars, i, simplify = FALSE))

A <- lapply(combinations, function(combo) lapply(combo, fit_model))
A <- unlist(A, recursive = FALSE)

#Compare models performance (AIC, BIC, Deviance)
compareGLM(A[[1]],A[[2]],A[[3]],A[[4]],A[[5]],
           A[[6]],A[[7]],A[[8]],A[[9]],A[[10]],
           A[[11]],A[[12]],A[[13]],A[[14]],A[[15]])

GLM_compare$Models
summary(A[[1]])

write.csv(GLM_compare$Fit.criteria, "GLMs_Plants.csv")

#Best model for each taxa
Brios<-glm(BRIO ~ BIO2 + Slope_2m, "poisson", 
           data=BD_biocrostas)
Lichens<- glm(LICHEN ~ Altitude_2m, "poisson", 
              data=BD_biocrostas)
Plants<- glm(PLANT ~ TWI + Ice_512 + Altitude_2m, "poisson", 
             data=BD_biocrostas)

summary(Lichens)


#
##Structural equation modelling (piecewiseSEM) ####
#(https://www.youtube.com/watch?v=VT-gw_VVP1E)

#upload files.
BD_biocrostas<-read.csv("BD_Geral.csv", sep = ";", row.names = 1)

#Original model sintax without biotic interactions
SEMglm1 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

#Check for suggested alterations in the regression models
dSep(SEMglm1, conserve = T)

##All remaining model sintax with biotic interactions
SEMglm2 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2 + PLANT, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m + PLANT, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

SEMglm3 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2 + PLANT, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

SEMglm4 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m + PLANT, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

SEMglm5 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + BRIO, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m + BRIO, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

SEMglm6 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + BRIO, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

SEMglm7 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m + BRIO, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

SEMglm8 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + LICHEN, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2 + LICHEN, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

SEMglm9 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + LICHEN, "poisson",
                    data=BD_biocrostas),
                glm(BRIO ~ Slope_2m + BIO2, "poisson",
                    data=BD_biocrostas),
                glm(LICHEN ~ Altitude_2m, "poisson", 
                    data=BD_biocrostas),
                data=BD_biocrostas)

SEMglm10 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm11 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + PLANT, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + PLANT + BRIO, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm12 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + PLANT + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + PLANT, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm13 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + PLANT + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm14 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + PLANT, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm15 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + PLANT, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + BRIO, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm16 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + BRIO, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + PLANT + BRIO, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm17 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + BRIO + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + BRIO, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm18 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + BRIO + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm19 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + BRIO, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm20 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + BRIO, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + PLANT, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm21 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + PLANT + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm22 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + LICHEN + BRIO, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm23 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m + PLANT + BRIO, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm24 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + BRIO, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

SEMglm25 <- psem(glm(PLANT ~ Ice_512 + Altitude_2m + TWI + LICHEN, "poisson",
                     data=BD_biocrostas),
                 glm(BRIO ~ Slope_2m + BIO2 + PLANT, "poisson",
                     data=BD_biocrostas),
                 glm(LICHEN ~ Altitude_2m, "poisson", 
                     data=BD_biocrostas),
                 data=BD_biocrostas)

#Overall model evaluation and (pseudo)R-square
summary(SEMglm1,conserve = T)

rsquared(SEMglm2, method = "mcfadden") #(>40%) model fits the data very well

#Check model residuals of endogenous variables 
SEM_residuals<- residuals(SEMglm1)

#Compare between SEMs based on Chi-Square
anova(SEMglm2, SEMglm3)

#Calculate Standard estimates for all direct pathways
##https://jslefche.github.io/sem_book/coefficients.html#scaling-to-other-non-normal-distributions (CHAPTER 4.5)

Glm_Brio <- glm(BRIO ~ BIO2 + Slope_2m + Altitude_2m + BIO15 + PLANT, "poisson", 
                data=BD_biocrostas)
Glm_Lichen <-glm(LICHEN ~ BIO2 + Ice_512 + Altitude_2m + BIO9 + PLANT + BRIO, "poisson", 
                 data=BD_biocrostas)
Glm_Plant <-glm(PLANT ~ TWI + Ice_512 + Altitude_2m + BIO9, "poisson", 
                data=BD_biocrostas)

coef(Glm_Brio) #extrair estimates de cada regressao

R2 <- cor(BD_biocrostas$BRIO, predict(Glm_Brio, type = "response"))^2 #entra a variavel dependente
sd.yhat <- sqrt(var(predict(Glm_Brio, type = "link"))/R2)

coefVALUE * sd(BD_biocrostas$BIO2)/sd.yhat #entra a variavel independente e respetivo estimate

#Calculate Standard estimates for all indirect pathways




##Boxplots ####
load(file="Boxplots.Rdata")

BRIO<-Boxplots[Boxplots$Vegetation == 'BRIO',]
LICHEN<-Boxplots[Boxplots$Vegetation == 'LICHEN',]
PLANT<-Boxplots[Boxplots$Vegetation == 'PLANT',]

ggplot(BRIO, aes(x=factor(Time, level=c('Sampled', 'Present', 'Future')), 
                 y=Abundance))+ 
   geom_boxplot(fill="cadetblue4")+
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.title = element_blank())

ggplot(LICHEN, aes(x=factor(Time, level=c('Sampled', 'Present', 'Future')), 
                   y=Abundance))+ 
   geom_boxplot(fill="darkgoldenrod2")+
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.title = element_blank())

ggplot(PLANT, aes(x=factor(Time, level=c('Sampled', 'Present', 'Future')), 
                  y=Abundance))+ 
   geom_boxplot(fill="chartreuse4")+
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.title = element_blank())
#
##Abundance shifts across altitudinal gradient####
BD_transect<-read.csv("Transect_lines.csv", sep = ";", row.names = 1)
Altitude_Hurd<-BD_transect[c(1:200),c(4,5)]

ggplot(Altitude_Hurd, aes(x=Distance_coast, y=Altitude)) +
   geom_smooth(fill="white", color="black")+
   ylim(0, 300)

BRIO<-BD_transect[BD_transect$Vegetation == 'BRIO',]
LICHEN<-BD_transect[BD_transect$Vegetation == 'LICHEN',]
PLANT<-BD_transect[BD_transect$Vegetation == 'PLANT',]

ggplot(PLANT, aes(x=Distance_coast, y=Abundance, 
                  group=Time, color=Time))+
   scale_color_manual(values=c("Present"="yellow", "Future"="green"))+
   geom_smooth(show.legend = F, size=1)+
   scale_x_reverse()+
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         axis.title = element_blank())+
   ylim(0,100)
#