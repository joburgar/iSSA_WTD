## Load R packages


install.packages("caret")
install.packages("glme")
install.packages("nlme")

library(caret)
library(glme)
library(nlme)


## Import deer_data_all CSV
all_data <- read.csv(file="deer_all_data.csv", header=T)
deerdata
colSums(is.na(all_data))


## getwd()

setwd("/Users/admin/Desktop/Alberta Boreal Deer Project/Siobhan Winter RSF files")
setwd("/Users/siobhandarlington-moore/Desktop")




##Create Composite Variables

all_data <- transform(all_data, Rail = pmin(RailLine, RailVegeta))
head(all_data)

all_data<- transform (all_data, Settlment = pmin(Cultivatio, MunicipalW, Urban, RuralResid))

all_data <- transform (all_data, IndustrialMine = pmin(BorrowPits, Industrial, MineSite, WellSite, PeatMine))

 all_data <- transform (all_data, OtherLinear = pmin(OneLaneGra, OneLaneUnd, TwoLaneGra, TwoLaneUnd, RailLineSp, SingleTrac, WinterRoad))
 
 all_data <- transform (all_data, Road = pmin(Unimproved, RoadVegeta, RoadHardSu, RoadTrailV))

all_data <- transform (all_data, NearestBlock = pmin(CutBlocks, HighDensit, OtherDistu, IndustrialMine, Settlment, WellSite))

all_data <- transform (all_data, NearestLinear = pmin(Cutline, Cutline3D, Electrical, Pipeline, Trail, TruckTrail, Rail, Road))

# # List of Composite Variables:
# 
# NearestBlock = CutBlocks, HighDensit, OtherDistu, IndustrialMine, Settlement, Wellsite
# Settlment = Cultivatio, MunicipalW, Urban, RuralResid
# IndustrialMine = BorrowPits, Industrial, MineSite, WellSite, PeatMine
# Rail = RailLine, RailVegeta
# OtherLinear = OneLaneGra, OneLaneUnd, TwoLaneGra, TwoLaneUnd, RailLineSp, SingleTrac, WinterRoad
# Road = RoadHardSu, RoadTrailV, RoadVegeta, Unimproved

# # Additional anthropogenic composite variables created:
# 
# NearestBlock = CutBlocks, HighDensit, OtherDistu, IndustrialMine, Settlement, Wellsite
# NearestLinear = Cutline, Cutline3D, Electrical, Pipeline, Trail, TruckTrail, Rail, Road


## Before subsetting - standardize and run collinearity and spearman rank test
##Standardize all global model covariates:

all_data$IndustrialMine.std <- scale(all_data$IndustrialMine) 
all_data$Reservoirs.std <- scale(all_data$Reservoirs)
all_data$CutBlocks.std <- scale(all_data$CutBlocks) 
all_data$OtherDistu.std <- scale(all_data$OtherDistu) 
all_data$WellSite.std <- scale(all_data$WellSite)
all_data$Cutline.std <- scale(all_data$Cutline)
all_data$Cutline3D.std <- scale(all_data$Cutline3D)
all_data$Electrical.std <- scale(all_data$Electrical)
all_data$Pipeline.std <- scale(all_data$Pipeline)
all_data$Trail.std <- scale(all_data$Trail)
all_data$TruckTrail.std <- scale(all_data$TruckTrail)
all_data$Settlment.std <- scale(all_data$Settlment)
all_data$Road.std <- scale(all_data$Road)
all_data$NearestBlock.std <- scale(all_data$NearestBlock)
all_data$NearestLinear.std <- scale(all_data$NearestLinear)
all_data$OtherLinear.std <- scale(all_data$OtherLinear)
all_data$Rail.std <- scale(all_data$Rail)
all_data$OtherDistu.std <- scale(all_data$OtherDistu)
all_data$HighDensit.std <- scale(all_data$HighDensit)

summary(Road)
Road <- as.numeric(as.character(all_data$Road))
#standardizing cover classes

all_data$PCT_AW.std <- scale(all_data$PCT_AW)
all_data$PCT_BW.std <- scale(all_data$PCT_BW)
all_data$PCT_PB.std <- scale(all_data$PCT_PB)
all_data$PCT_PJ.std <- scale(all_data$PCT_PJ)
all_data$PCT_SB.std <- scale(all_data$PCT_SB)
all_data$PCT_SW.std <- scale(all_data$PCT_SW)
all_data$PCT_FB.std <- scale(all_data$PCT_FB)
all_data$PCT_LT.std <- scale(all_data$PCT_LT)
all_data$uPCT_AW.std <- scale(all_data$uPCT_AW)
all_data$uPCT_BW.std <- scale(all_data$uPCT_BW)
all_data$uPCT_PB.std <- scale(all_data$uPCT_PB)    	
all_data$uPCT_PJ.std <- scale(all_data$uPCT_PJ)
all_data$uPCT_SB.std <- scale(all_data$uPCT_SB)
all_data$uPCT_SW.std <- scale(all_data$uPCT_SW)
all_data$uPCT_FB.std <- scale(all_data$uPCT_FB)
all_data$uPCT_LT.std <- scale(all_data$uPCT_LT)
head(all_data)

##Alex's standardizing code
stdize = function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))} 
Cutline.std <- stdize(all_data$Cutline, na.rm = T)
summary(all_data$Cutline.std)
all_data$Cutline.std

stdize= function(x,...){(x-min(x,...))/(max(x,...)-min(x,...))}
Pipeline.std <- stdize(all_data$Pipeline, na.rm = T)
summary(all_data$Pipeline.std)
summary(all_data$Cutline)

## Categorical data
Individual <- as.factor(all_data$CollarID)
Seralstage <- as.factor(all_data$SERAL_STGE)
FC_DOM <- as.factor(all_data$FC_DOM)
Height <- as.factor(all_data$HEIGHT)
Density <- as.factor(all_data$Density)
####



####
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,use="na.or.complete"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = 1.5)#cex.cor * r)
}

##Substrate PCA includes: organics.clay, sand, gravel, cobble, boulder
## pairs plot
pairs(Substrate.tot, upper.panel = panel.cor, cex=1.2, pch=16, labels=c("fines", "sand", "gravel", "cobble", "boulder"))






### New Data Frame

Use <- all_data$Use
Individual <- all_data$CollarID
IndustrialMine <- all_data$IndustrialMine.std 
Reservoirs <- all_data$Reservoirs.std 
CutBlocks <- all_data$CutBlocks.std 
OtherDistu <- all_data$OtherDistu.std  
WellSite <- all_data$WellSite.std 
Cutline <- all_data$Cutline.std 
NearestBlock <- all_data$NearestBlock.std
Road <- all_data$Road.std
Pipeline <- all_data$Pipeline.std
OtherLinear <- all_data$OtherLinear.std
NearestLinear <- all_data$NearestLinear.std
Cutline3D <- all_data$Cutline3D.std
Electrical <- all_data$Electrical.std
Trail <- all_data$Trail.std
Settlement <- all_data$Settlment.std
HighDensit <- all_data$HighDensit.std
PCT_AW <- all_data$PCT_AW.std 
PCT_BW <- all_data$PCT_BW.std 
PCT_PB <- all_data$PCT_PB.std 
PCT_PJ <- all_data$PCT_PJ.std 
PCT_SB <- all_data$PCT_SB.std 
PCT_SW <- all_data$PCT_SW.std 
PCT_FB <- all_data$PCT_FB.std 
PCT_LT <- all_data$PCT_LT.std 
uPCT_AW <- all_data$uPCT_AW.std 
uPCT_BW <- all_data$uPCT_BW.std 
uPCT_PB <- all_data$uPCT_PB.std     	
uPCT_PJ <- all_data$uPCT_PJ.std 
uPCT_SB <- all_data$uPCT_SB.std 
uPCT_SW <- all_data$uPCT_SW.std 
uPCT_FB <- all_data$uPCT_FB.std 
uPCT_LT <- all_data$uPCT_LT.std 

##Create New Data Frame

deerdata <- data.frame(IndustrialMine, Individual, Reservoirs, CutBlocks, OtherDistu, HighDensit, WellSite, Cutline, NearestBlock, NearestLinear, Road, Pipeline, OtherLinear, Cutline3D, Electrical, Trail, Settlement, PCT_AW, PCT_BW, PCT_PB, PCT_PJ, PCT_SB, PCT_SW, PCT_FB, PCT_LT, uPCT_AW, uPCT_BW, uPCT_PB, uPCT_PJ, uPCT_SB, uPCT_SW, uPCT_FB, uPCT_LT)

head(deerdata)

HF <- data.frame(IndustrialMine, Reservoirs, CutBlocks, OtherDistu, HighDensit, WellSite, Cutline, NearestBlock, NearestLinear, Road, Pipeline, OtherLinear, Cutline3D, Electrical, Trail, Settlement)

NF <- data.frame(PCT_AW, PCT_BW, PCT_PB, PCT_PJ, PCT_SB, PCT_SW, PCT_FB, PCT_LT, uPCT_AW, uPCT_BW, uPCT_PB, uPCT_PJ, uPCT_SB, uPCT_SW, uPCT_FB, uPCT_LT)


## Test for normality
shapiroo.test ## doesn't work for large data sets
can look at QQ plot


## Colinearity 

cor.test(Cutline, NearestLinear, method="spearman")
cor.test(CutBlocks, NearestBlock, method="spearman")
cor.test(Road, OtherLinear, method="spearman")
cor.test(Pipeline, OtherLinear, method="spearman")
cor.test(Cutline, OtherLinear, method="spearman")
plot(Use, all_data$CutBlocks)

cor.test(Road, NearestLinear, method="spearman")


Use <- as.factor(all_data$Use)
summary(Individual)
Individual <- as.factor(all_data$CollarID)
cor(all_data$Cutline.std, all_data$CutBlocks.std, method = "spearman")

## pairs plot
pairs(~ Road + Pipeline + Cutline + Trail + NearestLinear + OtherLinear + Electrical, data = HF,
  main = "Linear Features", method = "spearman")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="na.or.complete", method="spearman"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

Use <- all_data$Use
summary(Use)
## Variance Inflation Factor

full.model<-lm(Use~Cutline+CutBlocks+Pipeline+Road+Trail+NearestLinear+Settlement+Reservoirs+WellSite+OtherDistu+Cutline3D, data=HF)
summary(full.model)
vif(full.model)
##drop Industrial Mine, Electrical, Other linear, High Densit



full.model2 <- lm(Use~PCT_AW + PCT_BW + PCT_PB + PCT_PJ + PCT_SB + PCT_SW + PCT_FB + PCT_LT + uPCT_AW + uPCT_BW + uPCT_PB + uPCT_PJ + uPCT_SB + uPCT_SW + uPCT_FB + uPCT_LT, data=NF)
vif(full.model2)

##Global model
install.packages("lme4")
library(lme4)
globalHF <- glmer(Use ~ Cutline + CutBlocks + Pipeline + Road + Trail + NearestLinear + Settlement + Reservoirs + WellSite + OtherDistu + Cutline3D + (1 | Individual), deerdata, family= binomial)
summary(globalHF)   

globalHF.2 <- glmer(Use~Cutline + CutBlocks + Pipeline + Road + Trail + Settlement + WellSite + Cutline3D + (1|Individual), deerdata, family=binomial)
summary(globalHF.2)

globalHF.3 <- glmer(Use~ Cutline+ CutBlocks+ Pipeline+ Road+ Trail+ WellSite+ Cutline3D+ (1|Individual), deerdata, family=binomial)
summary(globalHF.3)

global.model <- glmer(Use~ Cutline+ CutBlocks+ Pipeline+ Road+ Trail+ WellSite+ Cutline3D+ PCT_AW + PCT_BW + PCT_PB + PCT_PJ + PCT_SB + PCT_SW + PCT_FB + PCT_LT + (1|Individual), deerdata, family=binomial)

globalNF <- glmer(Use~ PCT_AW + PCT_BW + PCT_PB + PCT_PJ + PCT_SB + PCT_SW + PCT_FB + PCT_LT + (1|Individual), data= deerdata, family=binomial)
summary(globalNF)


globalNFu <- glmer(Use~ uPCT_AW + uPCT_BW + uPCT_PB + uPCT_PJ + uPCT_SB + uPCT_SW + uPCT_FB + uPCT_LT + (1|Individual), data=deerdata, family = binomial)
summary(globalNFu)

##Drop uPCT_PJ, uPCT_BW, uPCT_AW

globalNFu.2 <- glmer(Use~ uPCT_PB + uPCT_SB + uPCT_SW + uPCT_FB + uPCT_LT + (1|Individual), data=deerdata, family = binomial)
summary(globalNFu.2)

####
set.seed(123)
inTraining <- createDataPartition(deerdata$FID, p= .8, list=FALSE)
training <- deerdata[inTraining,]
testing <- deerdata[-inTraining,]
training

describe()
structure(training)
summary(training)

####

> globalHF.3 <- glmer(Use~ Cutline+ CutBlocks+ Pipeline+ Road+ Trail+ WellSite+ Cutline3D+ (1|Individual), deerdata, family=binomial)
> summary(globalHF.3)
Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: binomial  ( logit )
Formula: Use ~ Cutline + CutBlocks + Pipeline + Road + Trail + WellSite +      Cutline3D + (1 | Individual)
   Data: deerdata

     AIC      BIC   logLik deviance df.resid 
 86404.6  86488.0 -43193.3  86386.6    78762 

Scaled residuals: 
     Min       1Q   Median       3Q      Max 
-12.2115  -0.6931  -0.0212   0.7148  14.1530 

Random effects:
 Groups     Name        Variance Std.Dev.
 Individual (Intercept) 0.4386   0.6623  
Number of obs: 78771, groups:  Individual, 34

Fixed effects:
             Estimate Std. Error z value Pr(>|z|)    
(Intercept) -0.049846   0.115590   -0.43    0.666    
Cutline     -0.188617   0.009801  -19.24   <2e-16 ***
CutBlocks    0.298646   0.011850   25.20   <2e-16 ***
Pipeline     0.585627   0.010919   53.63   <2e-16 ***
Road        -1.087295   0.014932  -72.82   <2e-16 ***
Trail       -0.384733   0.010072  -38.20   <2e-16 ***
WellSite    -0.430954   0.011563  -37.27   <2e-16 ***
Cutline3D    1.043590   0.013370   78.05   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
          (Intr) Cutlin CtBlck Pipeln Road   Trail  WellSt
Cutline    0.002                                          
CutBlocks -0.006 -0.019                                   
Pipeline   0.002  0.048  0.052                            
Road       0.014  0.091 -0.311 -0.244                     
Trail      0.003 -0.017 -0.160  0.138 -0.063              
WellSite  -0.002 -0.162  0.181 -0.291 -0.180 -0.071       
Cutline3D -0.009  0.088  0.412  0.143 -0.219 -0.082 -0.267


