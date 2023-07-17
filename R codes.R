#Open .csv file from the working directory
setwd("")
data=read.csv("Scheuermann.uninst.csv")
data #Displays data

#References for R codes (Baek et al., 2016; Moeyaert et al., 2020; Pustejovsky et al., 2014)
# Step 1: Effect size calculation 
#create phase-by-time interaction
data$phase_time <- with(data, 
                         unlist(tapply((phase=="B") * time, 
                                       list(phase,case), 
                                       function(x) x - min(x))))
data #Displays data

#Time-point constants
A <- 3
B <- 7

#Center at follow-up time
Center <- B
data$time <- data$time - Center

#Two-level model: varying intercepts, no trend, varying treatment effect, varying phase-by-time interaction 
install.packages("nlme")
library(nlme)
ctrl <- lmeControl(opt='optim')
two.level <- lme(fixed = outcome ~ phase + phase_time, 
                 random = ~ phase + phase_time | case,
                 correlation = corAR1(0, ~ time | case), control=ctrl,
                 data = data,
                 method = "REML")
summary(two.level)

VarCorr(two.level) #Extracts both variances and standard deviations from the lme object.
ICC.two.level <- as.numeric(VarCorr(two.level)[1:4])  # vector of variance estimates
ICC <- ICC.two.level[1]/sum(ICC.two.level[1],ICC.two.level[4])
ICC    #Displays value

#Reference for R codes (Pustejovsky et al., 2014) 
#Calculate between-case standardized mean difference (BC-SMD)
install.packages("scdhlm")
library(scdhlm)
g <- g_REML(m_fit = two.level, p_const = c(0,1, B - A),  r_const = c(1,0,1,0,0,0,0,0))
g[c("delta_AB", "g_AB","V_g_AB", "phi", "nu")]
CI_g(g)


#Step 2: Meta-analysis for effect sizes  
#Open .csv file from the working directory
setwd("")
WP=read.csv("WP.csv")
WP

#References for R codes (Hedges et al., 2010; Tipton, 2015; Tipton & Pustejovsky, 2015)
#robu() meta-analysis
install.packages("clubsandwich")
library(clubSandwich)

install.packages("robumeta")
library(robumeta)

#Intercept only (overall ES)
wp_intercept <- robu(formula = effect.size ~ 1, data = WP, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(wp_intercept)

#Create forest plot 
forest.robu(wp_intercept, es.lab = "Crit.Cat", study.lab = "Study", "BC-SMD" = effect.size, "Weight" = Weight) 

#Sensitivity analysis
sensitivity(wp_intercept)

#Publication bias
install.packages("metafor")
library(metafor)
res <- rma(effect.size, var.eff.size, data = WP)
regtest(res, model="lm")
funnel(res, xlab="BC-SMD")
regtest(res, model="lm", ret.fit=TRUE)


#Table 1. Subgroup analysis (CEC quality indicators) 
CEC1<- subset(WP, CEC1==1,)
CEC1_intercept <- robu(formula = effect.size ~ 1, data = CEC1, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(CEC1_intercept)

CEC2<- subset(WP, CEC2==1,)
CEC2_intercept <- robu(formula = effect.size ~ 1, data = CEC2, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(CEC2_intercept)

CEC3<- subset(WP, CEC3==1,)
CEC3_intercept <- robu(formula = effect.size ~ 1, data = CEC3, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(CEC3_intercept)

CEC4<- subset(WP, CEC4==1,)
CEC4_intercept <- robu(formula = effect.size ~ 1, data = CEC4, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(CEC4_intercept)

CEC5<- subset(WP, CEC5==1,)
CEC5_intercept <- robu(formula = effect.size ~ 1, data = CEC5, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(CEC5_intercept)

CEC6<- subset(WP, CEC6==1,)
CEC6_intercept <- robu(formula = effect.size ~ 1, data = CEC6, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(CEC6_intercept)

CEC7<- subset(WP, CEC7==1,)
CEC7_intercept <- robu(formula = effect.size ~ 1, data = CEC7, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(CEC7_intercept)

CEC8<- subset(WP, CEC8==1,)
CEC8_intercept <- robu(formula = effect.size ~ 1, data = CEC8, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(CEC8_intercept)


#Table 1. Subgroup analysis (CCSSM content standard) 
OA<- subset(WP, OA==1,)
OA_intercept <- robu(formula = effect.size ~ 1, data = OA, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(OA_intercept)

NF<- subset(WP, NF==1,)
NF_intercept <- robu(formula = effect.size ~ 1, data = NF, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(NF_intercept)

G<- subset(WP, G==1,)
G_intercept <- robu(formula = effect.size ~ 1, data = G, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(G_intercept)

RP<- subset(WP, RP==1,)
RP_intercept <- robu(formula = effect.size ~ 1, data = RP, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(RP_intercept)

NS<- subset(WP, NS==1,)
NS_intercept <- robu(formula = effect.size ~ 1, data = NS, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(NS_intercept)

EE<- subset(WP, EE==1,)
EE_intercept <- robu(formula = effect.size ~ 1, data = EE, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(EE_intercept)


#Table 1. Subgroup analysis (CCSSM practice standard) 
MP1 <- subset(WP, MP1==1,)
MP1_intercept <- robu(formula = effect.size ~ 1, data = MP1, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(MP1_intercept)

MP2 <- subset(WP, MP2==1,)
MP2_intercept <- robu(formula = effect.size ~ 1, data = MP2, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(MP2_intercept)

MP3 <- subset(WP, MP3==1,)
MP3_intercept <- robu(formula = effect.size ~ 1, data = MP3, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(MP3_intercept)

MP4 <- subset(WP, MP4==1,)
MP4_intercept <- robu(formula = effect.size ~ 1, data = MP4, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(MP4_intercept)

MP5 <- subset(WP, MP5==1,)
MP5_intercept <- robu(formula = effect.size ~ 1, data = MP5, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(MP5_intercept)

MP6 <- subset(WP, MP6==1,)
MP6_intercept <- robu(formula = effect.size ~ 1, data = MP6, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(MP6_intercept)

MP7 <- subset(WP, EE==1,)
MP7_intercept <- robu(formula = effect.size ~ 1, data = MP7, studynum = Study, var.eff.size = var.eff.size, rho = .8, small = TRUE)
print(MP7_intercept)


#Table 2. Meta-regression by CEC quality indicator (centered at the grand mean) 
wp_robu <- robu(effect.size ~ CEC.cnl,
                data = WP, studynum = Study, var.eff.size = var.eff.size, 
                modelweights = "HIER")
print(wp_robu, digits=3)


#Table 2. Meta-regression by CCSSM content standards (reference variable = OA) 
ccssm.c_robu <- robu(effect.size ~ NF + G + RP + NS + EE-1,  
                data = WP, studynum = Study, var.eff.size = var.eff.size, 
                modelweights = "HIER")
print(ccssm.c_robu)


#Table 2. Meta-regression by CCSSM practice standards (reference variable = 2.4.5) 
ccssm.p_robu <- robu(effect.size ~ MP1.2.4.5.6 + MP1.2.4.5.6.7 + MP1.2.3.4.5.6 + MP1.2.3.4.5.6.7-1,  
                data = WP, studynum = Study, var.eff.size = var.eff.size, 
                modelweights = "HIER")
print(ccssm.p_robu)

