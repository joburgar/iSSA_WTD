##Step Selection Function Code for R and GME

### Load packages
library(tidyverse)
library(MASS)
# library(lubridate)

install.packages("suncalc")
citation("suncalc")
##WTD telemetry data from 2012-2014
## For all Deer
Deer_all <- read.csv(file="Deer_used_ordered_ta_sl_cleaned_4.csv", header = T)

glimpse(Deer_all)
Date <- dmy(Deer_all$Date_)
Time <- as.factor(Deer_all$Time_)
DT <- ymd_hms(paste(Date,Time))
Deer_all$Datetime <- DT

min(Date)
max(Date)

# subset data to >360 observations and 30 sampling days in a season
# seasons = snow-free (May through Sep) and snow (Nov through Mar)
# exclude 4 and 10
head(Deer_all)
Deer_all %>% count(Month_)
Deer_all <- Deer_all %>% filter(!Month_ %in% c(4,10))
Deer_all <- Deer_all %>% mutate(Season = ifelse(Month_ < 4,"Snow", ifelse(Month_ > 10, "Snow","Snow-free")))
Deer_all %>% count(Season, Month_)
Deer_all$SeasonYr <- paste(Deer_all$Year_,Deer_all$Season, sep="_")
Deer_all %>% group_by(SeasonYr, CollarID) %>% count(CollarDay)

num_fixes <- Deer_all %>% count(CollarID)
num_fixes <- num_fixes %>% mutate(Use = ifelse(n>360,"Yes","No"))
num_fixes %>% count(Use)
deer_to_keep <- num_fixes %>% filter(Use=="Yes")

Deer_all <- Deer_all %>% filter(CollarID %in% deer_to_keep$CollarID)
Deer_all %>% count(CollarID)

# ##Need to generate random step lengths with Gamma distribution (Length only)
# COLID_32194 <- read.csv(file="COLID_32194_ta_sl.csv", header=T)
# # recalculate shape and scale parameters:
# max(COLID_32194$STEPLENGTH)
# mean(COLID_32194$STEPLENGTH) #129.545
# quantile(COLID_32194$STEPLENGTH, probs=c(0.25, 0.95)) # 95% quantile = 469.06968
# myfun <- function(shape) {
#   scale <- 129.545/shape
#   pgamma(469.06968, shape, scale=scale) - 0.95
# }
# 
# tmp <- uniroot( myfun, lower=0, upper=10 )
# 
# myshape <- tmp$root
# myscale <- 129.545/tmp$root
# myshape #0.5896872
# rate <- 1/0.5896872
# 
# rate #1.695814
# myscale #219.6843
# 
# hist(COLID_32194$STEPLENGTH)
# 
# COLID_32194_samples <- read.csv(file="COLID_32194_ssfsamples4.csv", header=T)
# hist(COLID_32194_samples$STEPLENGTH)
# mean((COLID_32194_samples$STEPLENGTH) & (COLID_32194_samples$OBSERVED<1))
# 
# 
# 
# Available <- COLID_32194_samples[(COLID_32194_samples$STEPLENGTH) & (COLID_32194_samples$OBSERVED<1), ]
# Used <- COLID_32194_samples[(COLID_32194_samples$STEPLENGTH) & (COLID_32194_samples$OBSERVED=1), ]
# 
# 
# min(Available$STEPLENGTH)
# max(Available$STEPLENGTH)
# mean(Available$STEPLENGTH)
# head(Available)
# hist(Available$STEPLENGTH)
# 
# t.test(Used$STEPLENGTH, Available$STEPLENGTH)
# mean(Used$STEPLENGTH)
# max(Used$STEPLENGTH)
# min(Used$STEPLENGTH)
# hist(Used$STEPLENGTH)
# 
# 
# ##Attach gamma dist. to dataframe
# 
# require(MASS)
# fitdistr(Length, "gamma")
# COL_32194 <- transform (COL_32194, StepL.gam = rgamma(2034, shape = 0.5871214, scale=219.6919))
# head(COL_32194)
# plot(Step.length.gamma)
# hist(Step.length.gamma)
# hist(rgamma(2046, shape= 0.5871214, scale=219.6919))
# 
# 
# 
# set.seed(123)
# fitdistr(Length, rgamma, list(shape= 0.5871214, rate =219.6919), lower=0.001)
# hist(rgamma(2046, shape=0.5871214, rate = 219.6919))
# 
# 
# ### gamma might not be a good fit - try a lognormal?
# 
# fitdistr(Angle, "log-normal")
# 
# 
# fitdistr(COLID_32194_ta_sl$STEPLENGTH, "log-normal")
# max(COLID_32194_ta_sl$STEPLENGTH)
# min(COLID_32194_ta_sl$STEPLENGTH)
# 
# ##meanlog       sdlog   
# #4.03115653   1.45637691
# #(0.03226649) (0.02281586)




## Sort by CollarID and then Datetime 
Deer_all <- Deer_all[with(Deer_all, order(Deer_all$CollarID, Deer_all$Datetime)),]
glimpse(Deer_all)
head(Deer_all)

# recalculate shape and scale parameters:
max(Deer_all$STEPLENGTH) #10744
aveSL <- mean(Deer_all$STEPLENGTH) #165.4324
quant95SL <- as.numeric(quantile(Deer_all$STEPLENGTH, probs=c(0.95))) # 95% quantile = 620.38043

myfun <- function(shape) {
  scale <- aveSL/shape
  pgamma(quant95SL, shape, scale=scale) - 0.95
}

tmp <- uniroot(myfun, lower=0, upper=10 )

myshape <- tmp$root
myscale <- aveSL/myshape
myshape # 0.5348888
rate <- 1/myshape
rate # 1.869547
myscale # 309.2838

## Lognormal for all deer
head(Deer_all)
hist(Deer_all$STEPLENGTH)
fitdistr(Deer_all$STEPLENGTH, "log-normal")
##meanlog        sdlog   
#4.115064196   1.573299329 
#(0.005754498) (0.004069045)
head(Deer_all)

## Run in GME
# R code provided by GME for movementssfsamples function - tool to generate available steps
movement.ssfsamples(in="Deer_all_used.gdb!Deer_used_ordered_ta_sl_cleaned", 
                    uidfield="FID_", order="OBJECTID", tad=c("UNIFORM", -3.14, 3.14), 
                    sld=c("GAMMA", myscale, rate), nsamples=2, 
                    out="Deer_all_ssfsamples1.shp", include=TRUE, radians=TRUE)

## Steps between collar IDs connected - must run individuals in GME separately. 
## Subset deer_all with the turning angles and steplengths generated into collar IDs.


##COLID 32194
COLID_32194_ta_sl <- Deer_all[ which(Deer_all$CollarID=='32194'),]
head(COLID_32194_ta_sl)
write.csv(COLID_32194_ta_sl, file = "COLID_32194_ta_sl.csv")
COLID_32194_ta_sl <- read.csv(file="COLID_32194_ta_sl.csv", header = T)


max(COLID_32194_ta_sl$STEPLENGTH)

##with shape and scale parameters from total deer
movement.ssfsamples(in="E:\AITF Computer\Alberta_Boreal_Deer_Project_2012-2014\SSF\Deer_all_used.gdb!Deer_used_ordered_ta_sl_cleaned", uidfield="FID_", order="OBJECTID", tad=c("UNIFORM", -3.14, 3.14), sld=c("GAMMA", 308.7792, 1.771425), nsamples=2, out="E:\AITF Computer\Alberta_Boreal_Deer_Project_2012-2014\SSF\COLID_32194_ssfsamples_deer_all.shp", include=TRUE, radians=TRUE, where="CollarID=32194");
#WHERE clause: CollarID=32194

#6120 sampled steps created.
## Joinged GME line file to DEER_ALL points file in ARC.
COLID_32194_steps <- read.csv(file="COLID_32194_ssfsamples5_log.csv", header=T)
head(COLID_32194_steps)
summary(COLID_32194_steps)
COLID_32194_steps <- read.csv(file="COLID_32194_steps.csv", header=T)


Available <- COLID_32194_steps[(COLID_32194_steps$STEPLENGTH) & (COLID_32194_steps$OBSERVED<1), ]
Used <- COLID_32194_steps[(COLID_32194_steps$STEPLENGTH) & (COLID_32194_steps$OBSERVED>0), ]
Used
plot(used.model)


cor.test(Used$STEPLENGTH, Used$TURNANGLE, method = 'spearman')

plot(Used$TURNANGLE)
used.test <- t(Used$STEPLENGTH ~ Used$TURNANGLE)

min(Available$STEPLENGTH)
max(Available$STEPLENGTH)
mean(Available$STEPLENGTH)
head(Available)
hist(Available$STEPLENGTH)

t.test(Used$STEPLENGTH, Available$STEPLENGTH)
mean(Used$STEPLENGTH)
max(Used$STEPLENGTH)
min(Used$STEPLENGTH)
hist(Used$STEPLENGTH)

##COLID32195
COLID_32195_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32195'),]
head(COLID_32195_ta_sl)
write.csv(COLID_32195_ta_sl, file = "COLID_32195_ta_sl.csv")

fitdistr(COLID_32195_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32195_ta_sl$STEPLENGTH)
min(COLID_32195_ta_sl$STEPLENGTH)

##meanlog       sdlog   
#4.68305084  1.42125210
#(0.03226649) (0.02281586)

##COLID32198
COLID_32198_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32198'),]
head(COLID_32198_ta_sl)
write.csv(COLID_32198_ta_sl, file = "COLID_32198_ta_sl.csv")

fitdistr(COLID_32198_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32198_ta_sl$STEPLENGTH)
min(COLID_32198_ta_sl$STEPLENGTH)

##meanlog       sdlog   
#3.12233423  1.65937000
#(0.03226649) (0.02281586)


##COLID32440
COLID_32440_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32440'),]
head(COLID_32440_ta_sl)
write.csv(COLID_32440_ta_sl, file = "COLID_32440_ta_sl.csv")

fitdistr(COLID_32440_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32440_ta_sl$STEPLENGTH)
min(COLID_32440_ta_sl$STEPLENGTH)

##meanlog       sdlog   
##4.32449883   1.44467315 
##(0.02609482) (0.01845183)

##COLID32787
COLID_32787_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32787'),]
head(COLID_32787_ta_sl)
write.csv(COLID_32787_ta_sl, file = "COLID_32787_ta_sl.csv")

fitdistr(COLID_32787_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32787_ta_sl$STEPLENGTH)
min(COLID_32787_ta_sl$STEPLENGTH)

##meanlog       sdlog   
##4.39051061   1.15986452 
##(0.03459586) (0.02446297)


##COLID32788
COLID_32788_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32788'),]
head(COLID_32788_ta_sl)
write.csv(COLID_32788_ta_sl, file = "COLID_32788_ta_sl.csv")

fitdistr(COLID_32788_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32788_ta_sl$STEPLENGTH)
min(COLID_32788_ta_sl$STEPLENGTH)

#meanlog       sdlog   
##3.65265154   1.65649825 
#(0.12346809) (0.08730512)

##COLID32789
COLID_32789_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32789'),]
head(COLID_32789_ta_sl)
write.csv(COLID_32789_ta_sl, file = "COLID_32789_ta_sl.csv")

fitdistr(COLID_32789_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32789_ta_sl$STEPLENGTH)
min(COLID_32789_ta_sl$STEPLENGTH)

    ## meanlog       sdlog   
   ## 3.63509331   1.55857379 
   ## (0.07754194) (0.05483043) 

    
##COLID32790
COLID_32790_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32790'),]
head(COLID_32790_ta_sl)
write.csv(COLID_32790_ta_sl, file = "COLID_32790_ta_sl.csv")

fitdistr(COLID_32790_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32790_ta_sl$STEPLENGTH)
min(COLID_32790_ta_sl$STEPLENGTH)

##meanlog       sdlog   
#3.94076289   1.45515339 
#(0.04197170) (0.02967847)

##COLID32791
COLID_32791_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32791'),]
head(COLID_32791_ta_sl)
write.csv(COLID_32791_ta_sl, file = "COLID_32791_ta_sl.csv")

fitdistr(COLID_32791_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32791_ta_sl$STEPLENGTH)
min(COLID_32791_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#4.24698965   1.59602908 
#(0.03272228) (0.02313815)


##COLID32792
COLID_32792_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32792'),]
head(COLID_32792_ta_sl)
write.csv(COLID_32792_ta_sl, file = "COLID_32792_ta_sl.csv")

fitdistr(COLID_32792_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32792_ta_sl$STEPLENGTH)
min(COLID_32792_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.68939558   1.38648250 
#(0.04755596) (0.03362714)


##COLID32793
COLID_32793_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32793'),]
head(COLID_32793_ta_sl)
write.csv(COLID_32793_ta_sl, file = "COLID_32793_ta_sl.csv")

fitdistr(COLID_32793_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32793_ta_sl$STEPLENGTH)
min(COLID_32793_ta_sl$STEPLENGTH)

## meanlog       sdlog   
#4.85151279   1.63474539 
#(0.03820375) (0.02701413)

##COLID32794
COLID_32794_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32794'),]
head(COLID_32794_ta_sl)
write.csv(COLID_32794_ta_sl, file = "COLID_32794_ta_sl.csv")

fitdistr(COLID_32794_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32794_ta_sl$STEPLENGTH)
min(COLID_32794_ta_sl$STEPLENGTH)

# meanlog       sdlog   
#4.58900981   1.72400079 
#(0.04461774) (0.03154951)

##COLID32795
COLID_32795_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32795'),]
head(COLID_32795_ta_sl)
write.csv(COLID_32795_ta_sl, file = "COLID_32795_ta_sl.csv")

fitdistr(COLID_32795_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32795_ta_sl$STEPLENGTH)
min(COLID_32795_ta_sl$STEPLENGTH)

#  meanlog       sdlog   
#3.05201085   1.52522863 
#(0.04718205) (0.03336275)

##COLID32796
COLID_32796_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32796'),]
head(COLID_32796_ta_sl)
write.csv(COLID_32796_ta_sl, file = "COLID_32796_ta_sl.csv")

fitdistr(COLID_32796_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32796_ta_sl$STEPLENGTH)
min(COLID_32796_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#4.04133924   1.59495278 
#(0.02034968) (0.01438940)

##COLID32797
COLID_32797_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32797'),]
head(COLID_32797_ta_sl)
write.csv(COLID_32797_ta_sl, file = "COLID_32797_ta_sl.csv")

fitdistr(COLID_32797_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32797_ta_sl$STEPLENGTH)
min(COLID_32797_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.11820405   1.45042933 
#(0.07198360) (0.05090009)


##COLID32798
COLID_32798_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32798'),]
head(COLID_32798_ta_sl)
write.csv(COLID_32798_ta_sl, file = "COLID_32798_ta_sl.csv")

fitdistr(COLID_32798_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32798_ta_sl$STEPLENGTH)
min(COLID_32798_ta_sl$STEPLENGTH)

# meanlog       sdlog   
#4.13347271   1.47352995 
#(0.01527733) (0.01080270)

##COLID32799
COLID_32799_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32799'),]
head(COLID_32799_ta_sl)
write.csv(COLID_32799_ta_sl, file = "COLID_32799_ta_sl.csv")

fitdistr(COLID_32799_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32799_ta_sl$STEPLENGTH)
min(COLID_32799_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.12230465   1.38231431 
#(0.06501860) (0.04597509)

#COLID327800
COLID_32800_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32800'),]
head(COLID_32800_ta_sl)
write.csv(COLID_32800_ta_sl, file = "COLID_32800_ta_sl.csv")

fitdistr(COLID_32800_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32800_ta_sl$STEPLENGTH)
min(COLID_32800_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.15116965   1.44446527 
#(0.04040561) (0.02857108)

##COLID327803
COLID_32803_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32803'),]
head(COLID_32803_ta_sl)
write.csv(COLID_32803_ta_sl, file = "COLID_32803_ta_sl.csv")

fitdistr(COLID_32803_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32803_ta_sl$STEPLENGTH)
min(COLID_32803_ta_sl$STEPLENGTH)

# meanlog       sdlog   
#2.70096271   1.32797998 
#(0.12383477) (0.08756441)

##COLID327805
COLID_32805_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32805'),]
head(COLID_32805_ta_sl)
write.csv(COLID_32805_ta_sl, file = "COLID_32805_ta_sl.csv")

fitdistr(COLID_32805_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32805_ta_sl$STEPLENGTH)
min(COLID_32805_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.92095031   1.60023204 
#(0.03590818) (0.02539091)

##COLID327805
COLID_32805_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32805'),]
head(COLID_32805_ta_sl)
write.csv(COLID_32805_ta_sl, file = "COLID_32805_ta_sl.csv")

fitdistr(COLID_32805_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32805_ta_sl$STEPLENGTH)
min(COLID_32805_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.92095031   1.60023204 
#(0.03590818) (0.02539091)

##COLID32196
COLID_32196_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32196'),]
head(COLID_32196_ta_sl)
write.csv(COLID_32196_ta_sl, file = "COLID_32196_ta_sl.csv")

fitdistr(COLID_32196_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32196_ta_sl$STEPLENGTH)
min(COLID_32196_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#2.83193805   1.18103539 
#(0.07045466) (0.04981897)

##COLID32197
COLID_32197_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32197'),]
head(COLID_32197_ta_sl)
write.csv(COLID_32197_ta_sl, file = "COLID_32197_ta_sl.csv")

fitdistr(COLID_32197_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32197_ta_sl$STEPLENGTH)
min(COLID_32197_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.09763308   1.29025841 
#(0.03477040) (0.02458639)

####COLID32806
COLID_32806_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32806'),]
head(COLID_32806_ta_sl)
write.csv(COLID_32806_ta_sl, file = "COLID_32806_ta_sl.csv")

fitdistr(COLID_32806_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32806_ta_sl$STEPLENGTH)
min(COLID_32806_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#4.29018981   1.58443546 
#(0.01982093) (0.01401552)

####COLID327807
COLID_32807_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32807'),]
head(COLID_32807_ta_sl)
write.csv(COLID_32807_ta_sl, file = "COLID_32807_ta_sl.csv")

fitdistr(COLID_32807_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32807_ta_sl$STEPLENGTH)
min(COLID_32807_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#4.22515830   1.61948279 
#(0.04359501) (0.03082633)

####COLID327809
COLID_32809_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32809'),]
head(COLID_32809_ta_sl)
write.csv(COLID_32809_ta_sl, file = "COLID_32809_ta_sl.csv")

fitdistr(COLID_32809_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32809_ta_sl$STEPLENGTH)
min(COLID_32809_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#4.33404655   1.53452782 
#(0.02056668) (0.01454284)

####COLID327810
COLID_32810_ta_sl <-Deer_all[ which(Deer_all$CollarID=='32810'),]
head(COLID_32810_ta_sl)
write.csv(COLID_32810_ta_sl, file = "COLID_32810_ta_sl.csv")

fitdistr(COLID_32810_ta_sl$STEPLENGTH, "log-normal")
max(COLID_32810_ta_sl$STEPLENGTH)
min(COLID_32810_ta_sl$STEPLENGTH)

## meanlog       sdlog   
#4.60335918   1.38141319 
#(0.03550264) (0.02510416)


###COLID33452
COLID_33452_ta_sl <-Deer_all[ which(Deer_all$CollarID=='33452'),]
head(COLID_33452_ta_sl)
write.csv(COLID_33452_ta_sl, file = "COLID_33452_ta_sl.csv")

fitdistr(COLID_33452_ta_sl$STEPLENGTH, "log-normal")
max(COLID_33452_ta_sl$STEPLENGTH)
min(COLID_33452_ta_sl$STEPLENGTH)

#meanlog      sdlog  
#3.1765333   1.8502685 
#(0.2312836) (0.1635422)

###COLID33453
COLID_33453_ta_sl <-Deer_all[ which(Deer_all$CollarID=='33453'),]
head(COLID_33453_ta_sl)
write.csv(COLID_33453_ta_sl, file = "COLID_33453_ta_sl.csv")

fitdistr(COLID_33453_ta_sl$STEPLENGTH, "log-normal")
max(COLID_33453_ta_sl$STEPLENGTH)
min(COLID_33453_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#2.94690986   1.49010582 
#(0.09386784) (0.06637459)

###COLID33454
COLID_33454_ta_sl <-Deer_all[ which(Deer_all$CollarID=='33454'),]
head(COLID_33454_ta_sl)
write.csv(COLID_33454_ta_sl, file = "COLID_33454_ta_sl.csv")

fitdistr(COLID_33454_ta_sl$STEPLENGTH, "log-normal")
max(COLID_33454_ta_sl$STEPLENGTH)
min(COLID_33454_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#4.25547379   1.59402911 
#(0.02416305) (0.01708586)

##COLID_33456
COLID_33456_ta_sl <-Deer_all[ which(Deer_all$CollarID=='33456'),]
head(COLID_33456_ta_sl)
write.csv(COLID_33456_ta_sl, file = "COLID_33456_ta_sl.csv")

fitdistr(COLID_33456_ta_sl$STEPLENGTH, "log-normal")
max(COLID_33456_ta_sl$STEPLENGTH)
min(COLID_33456_ta_sl$STEPLENGTH)


## meanlog       sdlog   
#4.07224381   1.57351900 
#(0.02579539) (0.01824010)

##COLID_33457
COLID_33457_ta_sl <-Deer_all[ which(Deer_all$CollarID=='33457'),]
head(COLID_33457_ta_sl)
write.csv(COLID_33457_ta_sl, file = "COLID_33457_ta_sl.csv")

fitdistr(COLID_33457_ta_sl$STEPLENGTH, "log-normal")
max(COLID_33457_ta_sl$STEPLENGTH)
min(COLID_33457_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.88656654   1.64162200 
#(0.02796099) (0.01977140)

##COLID_328052
COLID_328052_ta_sl <-Deer_all[ which(Deer_all$CollarID=='328052'),]
head(COLID_328052_ta_sl)
write.csv(COLID_328052_ta_sl, file = "COLID_328052_ta_sl.csv")

fitdistr(COLID_328052_ta_sl$STEPLENGTH, "log-normal")
max(COLID_328052_ta_sl$STEPLENGTH)
min(COLID_328052_ta_sl$STEPLENGTH)

#   meanlog       sdlog   
#3.97537766   1.54161565 
#(0.03071581) (0.02171936)

##COLID_328082
COLID_328082_ta_sl <-Deer_all[ which(Deer_all$CollarID=='328082'),]
head(COLID_328082_ta_sl)
write.csv(COLID_328082_ta_sl, file = "COLID_328082_ta_sl.csv")

fitdistr(COLID_328082_ta_sl$STEPLENGTH, "log-normal")
max(COLID_328082_ta_sl$STEPLENGTH)
min(COLID_328082_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#3.57736171   1.73761875 
#(0.06831284) (0.04830447)

##COLID_328083
COLID_328083_ta_sl <-Deer_all[ which(Deer_all$CollarID=='328083'),]
head(COLID_328083_ta_sl)
write.csv(COLID_328083_ta_sl, file = "COLID_328083_ta_sl.csv")

fitdistr(COLID_328083_ta_sl$STEPLENGTH, "log-normal")
max(COLID_328083_ta_sl$STEPLENGTH)
min(COLID_328083_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#4.00421909   1.50154211 
#(0.02006162) (0.01418571)

##COLID_334552
COLID_334552_ta_sl <-Deer_all[ which(Deer_all$CollarID=='334552'),]
head(COLID_334552_ta_sl)
write.csv(COLID_334552_ta_sl, file = "COLID_334552_ta_sl.csv")

fitdistr(COLID_334552_ta_sl$STEPLENGTH, "log-normal")
max(COLID_334552_ta_sl$STEPLENGTH)
min(COLID_334552_ta_sl$STEPLENGTH)

#meanlog       sdlog   
#4.06388548   1.55121553 
#(0.03003165) (0.02123559)
