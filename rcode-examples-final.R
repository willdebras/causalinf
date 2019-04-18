#Jason Roy's R script



##################
# Toy example
#install.packages("sandwich")
library(sandwich)

n<-200
L<-rbinom(n,1,0.5)
A<-rbinom(n,1,(0.3+0.3*L))
Y<-rnorm(n,(35+10*L+20*A),5)
mydata<-data.frame(cbind(A,L,Y))

#unadjusted
#mean of Y for treated subjects
mean(Y[A==1])
#mean of Y for untreated subjects
mean(Y[A==0])
#difference
mean(Y[A==1])-mean(Y[A==0])

#regression (ok here since continuous outcome, additive linear true model)
reg<-lm(Y~A+L,data=mydata)
summary(reg)

#g-comp / standardization
#can do by hand since one binary covariate
pr.l <- prop.table(table(mydata$L))
pr.l

# Table of Classified Means
tab.out <- aggregate(Y ~ A + L, mydata, mean)
tab.out

((mean(mydata$Y[mydata$A==1 & mydata$L==1]) - 
    mean(mydata$Y[mydata$A==0 & mydata$L==1]))*pr.l[2]) + 
  (mean(mydata$Y[mydata$A==1 & mydata$L==0]) - 
     mean(mydata$Y[mydata$A==0 & mydata$L==0]))*pr.l[1]


#IPTW
# Estimate Propensity Score
p.score <- glm(A ~ as.factor(L), data=mydata, family=binomial)
p.a <- ifelse(mydata$A == 0, 1 - predict(p.score, type = "response"),
                  predict(p.score, type = "response"))

# Table of probability of A
table(p.a)

# Generate IP Weights                     
mydata$w <- 1/p.a
table(mydata$w)
summary(mydata$w)
sd(mydata$w)

# Estimate ATE
mean(mydata$w*as.numeric(mydata$A==1)*mydata$Y) - 
  mean(mydata$w*as.numeric(mydata$A==0)*mydata$Y)

# marginal structural model MSM
msm <- lm(Y  ~ A, data = mydata, weights = w)

# Now with Inference
SE <-sqrt(diag(vcovHC(msm, type="HC0"))) # robust standard errors
beta <- coef(msm)
lcl <- beta-1.96*SE 
ucl <- beta+1.96*SE
cbind(beta, lcl, ucl)[2,]



###################
###################
#RHC Example
###################


#load packages
#install.packages(tableone)
#install.packages(Matching)
#install.packages(geepack)
#install.packages(stats)
library(tableone)
library(Matching)
library(geepack)
library(stats)

#read in data
load(url("https://urldefense.proofpoint.com/v2/url?u=http-3A__biostat.mc.vanderbilt.edu_wiki_pub_Main_DataSets_rhc.sav&d=DwIGAg&c=lEzKI_JJakPtcnbAQ6Q5xQ&r=ADBnQF_uxktNSAkWyomgLeFb-eOkP5LRzvDAXj-4ESQ&m=GJ1JdJTpKzJvLWlb0bGX87lB6Gx5rPtIV9m0KvRo6AI&s=639__Kuc3BOyiS-CFnTwBH3AVQyrBK3GMJStuyusjrk&e="))
#view data
View(rhc)

#treatment variables is swang1
#x variables that we will use
#cat1: primary disease category
#age
#sex
#aps1: APACHE score
#meanbp1: mean blood pressure

#create a data set with just these variables, for simplicity
ARF<-as.numeric(rhc$cat1=='ARF')
CHF<-as.numeric(rhc$cat1=='CHF')
Cirr<-as.numeric(rhc$cat1=='Cirrhosis')
colcan<-as.numeric(rhc$cat1=='Colon Cancer')
Coma<-as.numeric(rhc$cat1=='Coma')
COPD<-as.numeric(rhc$cat1=='COPD')
lungcan<-as.numeric(rhc$cat1=='Lung Cancer')
MOSF<-as.numeric(rhc$cat1=='MOSF w/Malignancy')
sepsis<-as.numeric(rhc$cat1=='MOSF w/Sepsis')
female<-as.numeric(rhc$sex=='Female')
died<-as.numeric(rhc$death=='Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1=='RHC')
meanbp1<-rhc$meanbp1

#new dataset
mydata<-cbind(ARF,CHF,Cirr,colcan,Coma,lungcan,MOSF,sepsis,
              age,female,meanbp1,treatment,died)
mydata<-data.frame(mydata)

#covariates we will use (shorter list than you would use in practice)
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

#look at a table 1
table1<- CreateTableOne(vars=xvars,strata="treatment", data=mydata, test=FALSE)
## include standardized mean difference (SMD)
print(table1,smd=TRUE)

#unadjusted mean difference
mean(died[treatment==1])
mean(died[treatment==0])

#########################################
#g-computation / standardization
#########################################

#logistic regression for E(Y|A,L)

outmodel<-glm(died~treatment+ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1,
             family=binomial(),data=mydata)

#now get predicted values for each person if they had been treated
mydata1<-mydata
mydata1$treatment<-1
pred1<-predict(outmodel,mydata1,type='response')

#get predicted values for each person if they had not been treated
mydata0<-mydata
mydata0$treatment<-0
pred0<-predict(outmodel,mydata0,type='response')

mean(pred1)
mean(pred0)

##########################################
#IPTW
#########################################

#propensity score

psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1,
             family=binomial(),data=mydata)

#show coefficients etc
summary(psmodel)
#create propensity score
pscore<-psmodel$fitted.values


#unstabilized weights
wt<-ifelse(treatment==1,1/pscore,1/(1-pscore))
hist(wt)
summary(wt)

#stabilized weights
#probability of treatment (unconditional)
probt<-mean(treatment)
sw<-ifelse(treatment==1,probt/pscore,(1-probt)/(1-pscore))
hist(sw)
summary(sw)

#IPTW estimator of E(Y^1-Y^0)

#directly computed weighted means
weighted.mean(died[treatment==1],wt[treatment==1])
weighted.mean(died[treatment==0],wt[treatment==0])

#now do weighted regression (get same estimates)
id<-1:length(died)
msm1 <- geeglm(as.numeric(died) ~ treatment, 
               id = id,
               weights = wt,
               family = binomial(link=logit))
summary(msm1)
msm1$fitted.values


#using stabilized weights
msm2 <- geeglm(as.numeric(died) ~ treatment, 
               id = id,
               weights = sw,
               data=mydata,
               family = binomial(link=logit))
summary(msm2)
msm2$fitted.values

############################################
#do greedy matching on Mahalanobis distance
############################################

greedymatch<-Match(Tr=treatment,M=1,X=mydata[xvars],replace=FALSE)
matched<-mydata[unlist(greedymatch[c("index.treated","index.control")]), ]

#get table 1 for matched data with standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)

#outcome analysis
y_trt<-matched$died[matched$treatment==1]
y_con<-matched$died[matched$treatment==0]
mean(y_trt)
mean(y_con)

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)

#McNemar test
table(y_trt,y_con)

mcnemar.test(matrix(c(973,513,395,303),2,2))



##########################
#propensity score matching
#########################

#fit a propensity score model. logistic regression

psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1+aps,
    family=binomial(),data=mydata)

#show coefficients etc
summary(psmodel)
#create propensity score
pscore<-psmodel$fitted.values


#do greedy matching on logit(PS) using Match with a caliper

logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treatment,M=1,X=logit(pscore),replace=FALSE,caliper=.2)
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

#get standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)

#outcome analysis
y_trt<-matched$died[matched$treatment==1]
y_con<-matched$died[matched$treatment==0]
mean(y_trt)
mean(y_con)

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)

