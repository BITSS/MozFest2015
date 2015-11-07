#Replicating Young 2009 to learn how to do logit with clustered se in R
#June 23, 2015

setwd("~/Dropbox/Cambridge_University/Data_Downloads/Replication/Young2009")

#load Young's data
library(foreign)
mydata <- read.dta("young_humanrights.dta") 

summary(mydata)

#   ***********************
#   *Ordered Logit Models TABLE 2                   *
#   *************************************************
#   
#   ************************************
#   * Baseline Model 1a                *
#   * DV PHYSINT                       *
#   * OLOGIT                           *
#   ************************************
#   #delimit ;
#   ologit physintflip polity2 gdpenl jobinsecurity_bdm
# lpop diss_lag cowwarl lmtnest physint_lag, nolog cluster(cowcode); 


# trying with lrm
library(rms)
res1 <- lrm(physintflip ~ polity2 + gdpenl + jobinsecurity_bdm +
                      +lpop + diss_lag + cowwarl + lmtnest + physint_lag,data=mydata,x=TRUE, y=TRUE)
print(res1) 

g <- robcov(res1,mydata$cowcode) #robust errors
print(g) #good match with Young 2009 (estimates, se, observations, countries)

library(texreg,include.pseudors = F, 
        include.lr = F, include.nobs = TRUE,)
t <- extract(g)

# get pseudo Rsquared 
#install.packages("pscl") #http://www.inside-r.org/packages/cran/pscl/docs/pR2
library(pscl)
pR2(g) # 0.303534 suddenly McFadden is right! -- I tried it again and it is 0.62!?!? not reproducible!

# extract for latex, rename variables without Pseudo R
varnames <- c("DEMOC", "GDP", "INSECURE", "POP", "DISSACT","WAR","MOUNTAINS","PAST")
dvnames <- "PHYSINT"
stargazer(g,type="text",omit=c("y>=1","y>=2","y>=3","y>=4","y>=5","y>=6","y>=7","y>=8"),omit.stat=c("rsq","chi2"), covariate.labels=varnames, dep.var.labels=dvnames)

# calculate McFadden myself http://stat.ethz.ch/R-manual/R-patched/library/stats/html/logLik.html
res0 <- lrm(physintflip ~ 1,data=mydata,x=TRUE, y=TRUE)
McF.pR2 <- 1 - (logLik(res1)[1]/logLik(res0)[1] )  #-2064.157 #-5493.97 
# result: 0.6242868 as in pR2 package (they also give the loglikelihood for null model and my model = same numbers)

# try out adjusted McFADDEN (see my code "R-squared logistic.R")
p <- length(res1$coeff) 
adjusted <- 1 - (logLik(res1)[1]-p-1)/logLik(res0)[1] 
#0.6211925

# stata provides McFadden's R-squared http://www.ats.ucla.edu/stat/stata/output/stata_ologit_output.htm although stata manual itself does not say so http://www.stata.com/manuals13/rologit.pdf

# Summary: I could reproduce everything except R-squared




# I replicated some things in STATA as well (he did not provide the code, but I did the commands I think are right)

# Binary Decomposition 
res1 <- lrm(physintflip ~ factor(polity2) + gdpenl + jobinsecurity_bdm +
              +lpop + diss_lag + cowwarl + lmtnest + physint_lag,data=mydata,x=TRUE, y=TRUE)
print(res1)
g <- robcov(res1,mydata$cowcode) #robust errors
print(g) #good match with Young 2009 (estimates, se, observations, countries)

library(texreg,include.pseudors = F, 
        include.lr = F, include.nobs = TRUE,)
t <- extract(g)

# get pseudo Rsquared 
#install.packages("pscl") #http://www.inside-r.org/packages/cran/pscl/docs/pR2
library(pscl)
pR2(g) # 0.303534 suddenly McFadden is right! -- I tried it again and it is 0.62!?!? not reproducible!

# extract for latex, rename variables without Pseudo R
varnames <- c("DEMOC", "GDP", "INSECURE", "POP", "DISSACT","WAR","MOUNTAINS","PAST")
dvnames <- "PHYSINT"
stargazer(g,type="text",omit.stat=c("rsq","chi2"), covariate.labels=varnames, dep.var.labels=dvnames)


# 
# testing if coefficients are significantly different from each other (H1: b !=b2 H0: b=B2)
mydata$polity2 <- factor(mydata$polity2)
mydata<-mydata[!is.na(mydata$polity2),]
res1 <- lrm(physintflip ~ polity2 + gdpenl + jobinsecurity_bdm +
              +lpop + diss_lag + cowwarl + lmtnest + physint_lag,data=mydata,x=TRUE, y=TRUE)
print(res1)
vcov(res1)
#install.packages("multcomp")
library(multcomp)
summary(glht(res1, mcp(polity2="Tukey"))) #error
# Note: In STATA this is the "test" command; in R this is done with package multcomp http://stats.stackexchange.com/questions/59085/how-to-test-for-simultaneous-equality-of-choosen-coefficients-in-logit-or-probit




# # # # # # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # #
# Try rerunning the ologit with polr instead of lrm (ignore the clustered se for now that I can only get with lrm) as Torsten from multcomp package advised in email
# # # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # #
# email to multcomp package author, who said: "rms does not implement proper coef, vcov and model.matrix methods:

# length(coef(res1))
# [1] 35
# dim(model.matrix(res1))
# [1] 1395   28
# dim(vcov(res1))
# [1] 35 35
# 
# as their dimensions differ. You need to set-up the contrasts manually or use MASS::polr.
# 
# Best, Torsten

#trying same with polr
mydata$polity2 <- factor(mydata$polity2)
mydata$physintflip <- factor(mydata$physintflip)
res2 <- polr(physintflip ~ polity2 + gdpenl + jobinsecurity_bdm +
              +lpop + diss_lag + cowwarl + lmtnest + physint_lag,data=mydata)
print(res2)
summary(res2) #very similar to lrm but not exactly the same

library(multcomp)
summary(glht(res2, mcp(polity2="Tukey"))) #error




# Re-fitting to get Hessian
# 
# 
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: polr(formula = physintflip ~ polity2 + gdpenl + jobinsecurity_bdm + 
#             +lpop + diss_lag + cowwarl + lmtnest + physint_lag, data = mydata)
# 
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# -9 - -10 == 0  0.111012   0.326772   0.340   1.0000    
# -8 - -10 == 0 -0.040347   0.318033  -0.127   1.0000    
# -7 - -10 == 0 -0.329782   0.300797  -1.096   0.9999    
# -6 - -10 == 0  0.299993   0.361516   0.830   1.0000    
# -5 - -10 == 0 -0.097113   0.387910  -0.250   1.0000    
# -4 - -10 == 0 -0.709839   0.576453  -1.231   0.9996    
# -3 - -10 == 0  0.011601   0.499648   0.023   1.0000    
# -2 - -10 == 0 -0.700438   0.416555  -1.682   0.9802    
# -1 - -10 == 0 -0.242907   0.462872  -0.525   1.0000    
# 0 - -10 == 0  -0.115153   0.663139  -0.174   1.0000    
# 1 - -10 == 0  -0.599146   0.676477  -0.886   1.0000    
# 2 - -10 == 0  -1.256817   0.606168  -2.073   0.8598    
# 3 - -10 == 0   0.323159   0.583533   0.554   1.0000    
# 4 - -10 == 0  -0.718312   0.471234  -1.524   0.9937    
# 5 - -10 == 0  -0.145112   0.415429  -0.349   1.0000    
# 6 - -10 == 0  -0.698226   0.363784  -1.919   0.9259    
# 7 - -10 == 0  -0.757672   0.398205  -1.903   0.9313    
# 8 - -10 == 0  -0.662418   0.347034  -1.909   0.9302    
# 9 - -10 == 0  -0.572909   0.323665  -1.770   0.9660    
# 10 - -10 == 0 -2.083121   0.311001  -6.698    <0.01 ***
#   -8 - -9 == 0  -0.151359   0.267312  -0.566   1.0000    
# -7 - -9 == 0  -0.440794   0.216407  -2.037   0.8777    
# -6 - -9 == 0   0.188981   0.299332   0.631   1.0000    
# -5 - -9 == 0  -0.208125   0.322044  -0.646   1.0000    
# -4 - -9 == 0  -0.820851   0.538409  -1.525   0.9936    
# -3 - -9 == 0  -0.099411   0.452906  -0.219   1.0000    
# -2 - -9 == 0  -0.811450   0.380241  -2.134   0.8256    
# -1 - -9 == 0  -0.353919   0.412853  -0.857   1.0000    
# 0 - -9 == 0   -0.226165   0.630884  -0.358   1.0000    
# 1 - -9 == 0   -0.710158   0.644537  -1.102   0.9999    
# 2 - -9 == 0   -1.367829   0.568407  -2.406   0.6324    
# 3 - -9 == 0    0.212147   0.546111   0.388   1.0000    
# 4 - -9 == 0   -0.829324   0.431973  -1.920   0.9263    
# 5 - -9 == 0   -0.256125   0.366089  -0.700   1.0000    
# 6 - -9 == 0   -0.809238   0.311785  -2.596   0.4830    
# 7 - -9 == 0   -0.868684   0.343469  -2.529   0.5350    
# 8 - -9 == 0   -0.773430   0.281713  -2.745   0.3687    
# 9 - -9 == 0   -0.683921   0.273955  -2.496   0.5611    
# 10 - -9 == 0  -2.194134   0.281746  -7.788    <0.01 ***
#   -7 - -8 == 0  -0.289435   0.226905  -1.276   0.9994    
# -6 - -8 == 0   0.340340   0.306920   1.109   0.9999    
# -5 - -8 == 0  -0.056765   0.332669  -0.171   1.0000    
# -4 - -8 == 0  -0.669492   0.543555  -1.232   0.9996    
# -3 - -8 == 0   0.051948   0.459510   0.113   1.0000    
# -2 - -8 == 0  -0.660091   0.380021  -1.737   0.9720    
# -1 - -8 == 0  -0.202559   0.419166  -0.483   1.0000    
# 0 - -8 == 0   -0.074806   0.634323  -0.118   1.0000    
# 1 - -8 == 0   -0.558799   0.646677  -0.864   1.0000    
# 2 - -8 == 0   -1.216469   0.573165  -2.122   0.8315    
# 3 - -8 == 0    0.363506   0.550207   0.661   1.0000    
# 4 - -8 == 0   -0.677965   0.434598  -1.560   0.9916    
# 5 - -8 == 0   -0.104765   0.368238  -0.285   1.0000    
# 6 - -8 == 0   -0.657879   0.313081  -2.101   0.8446    
# 7 - -8 == 0   -0.717325   0.346099  -2.073   0.8607    
# 8 - -8 == 0   -0.622070   0.285138  -2.182   0.7966    
# 9 - -8 == 0   -0.532562   0.274287  -1.942   0.9186    
# 10 - -8 == 0  -2.042774   0.270477  -7.552    <0.01 ***
#   -6 - -7 == 0   0.629775   0.263417   2.391   0.6468    
# -5 - -7 == 0   0.232669   0.284115   0.819   1.0000    
# -4 - -7 == 0  -0.380057   0.517831  -0.734   1.0000    
# -3 - -7 == 0   0.341383   0.428301   0.797   1.0000    
# -2 - -7 == 0  -0.370656   0.351694  -1.054   1.0000    
# -1 - -7 == 0   0.086875   0.385500   0.225   1.0000    
# 0 - -7 == 0    0.214629   0.614903   0.349   1.0000    
# 1 - -7 == 0   -0.269364   0.626942  -0.430   1.0000    
# 2 - -7 == 0   -0.927035   0.545358  -1.700   0.9781    
# 3 - -7 == 0    0.652941   0.525323   1.243   0.9996    
# 4 - -7 == 0   -0.388530   0.405758  -0.958   1.0000    
# 5 - -7 == 0    0.184670   0.332743   0.555   1.0000    
# 6 - -7 == 0   -0.368444   0.271996  -1.355   0.9986    
# 7 - -7 == 0   -0.427890   0.305283  -1.402   0.9979    
# 8 - -7 == 0   -0.332636   0.233325  -1.426   0.9973    
# 9 - -7 == 0   -0.243127   0.235182  -1.034   1.0000    
# 10 - -7 == 0  -1.753339   0.244373  -7.175    <0.01 ***
#   -5 - -6 == 0  -0.397106   0.356239  -1.115   0.9999    
# -4 - -6 == 0  -1.009832   0.559511  -1.805   0.9589    
# -3 - -6 == 0  -0.288392   0.476702  -0.605   1.0000    
# -2 - -6 == 0  -1.000431   0.411560  -2.431   0.6149    
# -1 - -6 == 0  -0.542900   0.440319  -1.233   0.9996    
# 0 - -6 == 0   -0.415146   0.646410  -0.642   1.0000    
# 1 - -6 == 0   -0.899139   0.659334  -1.364   0.9985    
# 2 - -6 == 0   -1.556810   0.588628  -2.645   0.4436    
# 3 - -6 == 0    0.023166   0.562557   0.041   1.0000    
# 4 - -6 == 0   -1.018305   0.456685  -2.230   0.7641    
# 5 - -6 == 0   -0.445105   0.393732  -1.130   0.9999    
# 6 - -6 == 0   -0.998219   0.344390  -2.899   0.2673    
# 7 - -6 == 0   -1.057665   0.374905  -2.821   0.3167    
# 8 - -6 == 0   -0.962411   0.317026  -3.036   0.1927    
# 9 - -6 == 0   -0.872902   0.315870  -2.763   0.3561    
# 10 - -6 == 0  -2.383114   0.325601  -7.319    <0.01 ***
#   -4 - -5 == 0  -0.612726   0.570704  -1.074   1.0000    
# -3 - -5 == 0   0.108714   0.486180   0.224   1.0000    
# -2 - -5 == 0  -0.603325   0.427224  -1.412   0.9976    
# -1 - -5 == 0  -0.145794   0.452819  -0.322   1.0000    
# 0 - -5 == 0   -0.018040   0.658425  -0.027   1.0000    
# 1 - -5 == 0   -0.502034   0.671802  -0.747   1.0000    
# 2 - -5 == 0   -1.159704   0.593160  -1.955   0.9131    
# 3 - -5 == 0    0.420271   0.577458   0.728   1.0000    
# 4 - -5 == 0   -0.621199   0.470845  -1.319   0.9991    
# 5 - -5 == 0   -0.048000   0.412053  -0.116   1.0000    
# 6 - -5 == 0   -0.601113   0.365869  -1.643   0.9848    
# 7 - -5 == 0   -0.660560   0.389615  -1.695   0.9784    
# 8 - -5 == 0   -0.565305   0.333778  -1.694   0.9785    
# 9 - -5 == 0   -0.475796   0.337211  -1.411   0.9976    
# 10 - -5 == 0  -1.986009   0.343679  -5.779    <0.01 ***
#   -3 - -4 == 0   0.721440   0.653832   1.103   0.9999    
# -2 - -4 == 0   0.009401   0.604987   0.016   1.0000    
# -1 - -4 == 0   0.466932   0.628663   0.743   1.0000    
# 0 - -4 == 0    0.594686   0.787502   0.755   1.0000    
# 1 - -4 == 0    0.110693   0.799982   0.138   1.0000    
# 2 - -4 == 0   -0.546978   0.735818  -0.743   1.0000    
# 3 - -4 == 0    1.032998   0.720168   1.434   0.9970    
# 4 - -4 == 0   -0.008473   0.637609  -0.013   1.0000    
# 5 - -4 == 0    0.564726   0.596449   0.947   1.0000    
# 6 - -4 == 0    0.011613   0.563682   0.021   1.0000    
# 7 - -4 == 0   -0.047833   0.581923  -0.082   1.0000    
# 8 - -4 == 0    0.047421   0.546638   0.087   1.0000    
# 9 - -4 == 0    0.136930   0.545620   0.251   1.0000    
# 10 - -4 == 0  -1.373283   0.548407  -2.504   0.5555    
# -2 - -3 == 0  -0.712039   0.532579  -1.337   0.9989    
# -1 - -3 == 0  -0.254508   0.555476  -0.458   1.0000    
# 0 - -3 == 0   -0.126754   0.727883  -0.174   1.0000    
# 1 - -3 == 0   -0.610747   0.741947  -0.823   1.0000    
# 2 - -3 == 0   -1.268418   0.673337  -1.884   0.9374    
# 3 - -3 == 0    0.311558   0.654482   0.476   1.0000    
# 4 - -3 == 0   -0.729913   0.565795  -1.290   0.9993    
# 5 - -3 == 0   -0.156713   0.519711  -0.302   1.0000    
# 6 - -3 == 0   -0.709827   0.484151  -1.466   0.9961    
# 7 - -3 == 0   -0.769273   0.503134  -1.529   0.9934    
# 8 - -3 == 0   -0.674019   0.457059  -1.475   0.9958    
# 9 - -3 == 0   -0.584510   0.461331  -1.267   0.9995    
# 10 - -3 == 0  -2.094722   0.461991  -4.534    <0.01 ***
#   -1 - -2 == 0   0.457531   0.500904   0.913   1.0000    
# 0 - -2 == 0    0.585285   0.690362   0.848   1.0000    
# 1 - -2 == 0    0.101292   0.705632   0.144   1.0000    
# 2 - -2 == 0   -0.556379   0.629667  -0.884   1.0000    
# 3 - -2 == 0    1.023597   0.615360   1.663   0.9823    
# 4 - -2 == 0   -0.017874   0.508332  -0.035   1.0000    
# 5 - -2 == 0    0.555325   0.457899   1.213   0.9997    
# 6 - -2 == 0    0.002212   0.413339   0.005   1.0000    
# 7 - -2 == 0   -0.057235   0.437273  -0.131   1.0000    
# 8 - -2 == 0    0.038020   0.389948   0.098   1.0000    
# 9 - -2 == 0    0.127529   0.381662   0.334   1.0000    
# 10 - -2 == 0  -1.382684   0.370178  -3.735   0.0232 *  
#   0 - -1 == 0    0.127754   0.710640   0.180   1.0000    
# 1 - -1 == 0   -0.356240   0.716448  -0.497   1.0000    
# 2 - -1 == 0   -1.013910   0.651494  -1.556   0.9918    
# 3 - -1 == 0    0.566065   0.634860   0.892   1.0000    
# 4 - -1 == 0   -0.475405   0.541297  -0.878   1.0000    
# 5 - -1 == 0    0.097794   0.486028   0.201   1.0000    
# 6 - -1 == 0   -0.455319   0.447554  -1.017   1.0000    
# 7 - -1 == 0   -0.514766   0.467802  -1.100   0.9999    
# 8 - -1 == 0   -0.419511   0.427076  -0.982   1.0000    
# 9 - -1 == 0   -0.330002   0.424658  -0.777   1.0000    
# 10 - -1 == 0  -1.840215   0.434493  -4.235    <0.01 ** 
#   1 - 0 == 0    -0.483993   0.860597  -0.562   1.0000    
# 2 - 0 == 0    -1.141663   0.808043  -1.413   0.9976    
# 3 - 0 == 0     0.438312   0.787237   0.557   1.0000    
# 4 - 0 == 0    -0.603159   0.714427  -0.844   1.0000    
# 5 - 0 == 0    -0.029959   0.680063  -0.044   1.0000    
# 6 - 0 == 0    -0.583073   0.653499  -0.892   1.0000    
# 7 - 0 == 0    -0.642519   0.669392  -0.960   1.0000    
# 8 - 0 == 0    -0.547265   0.634173  -0.863   1.0000    
# 9 - 0 == 0    -0.457756   0.638643  -0.717   1.0000    
# 10 - 0 == 0   -1.967968   0.637914  -3.085   0.1702    
# 2 - 1 == 0    -0.657670   0.818488  -0.804   1.0000    
# 3 - 1 == 0     0.922305   0.797278   1.157   0.9999    
# 4 - 1 == 0    -0.119165   0.730820  -0.163   1.0000    
# 5 - 1 == 0     0.454034   0.690111   0.658   1.0000    
# 6 - 1 == 0    -0.099080   0.664052  -0.149   1.0000    
# 7 - 1 == 0    -0.158526   0.679330  -0.233   1.0000    
# 8 - 1 == 0    -0.063271   0.648878  -0.098   1.0000    
# 9 - 1 == 0     0.026237   0.651603   0.040   1.0000    
# 10 - 1 == 0   -1.483975   0.658850  -2.252   0.7494    
# 3 - 2 == 0     1.579975   0.737212   2.143   0.8213    
# 4 - 2 == 0     0.538505   0.659196   0.817   1.0000    
# 5 - 2 == 0     1.111704   0.620337   1.792   0.9617    
# 6 - 2 == 0     0.558591   0.589157   0.948   1.0000    
# 7 - 2 == 0     0.499144   0.603543   0.827   1.0000    
# 8 - 2 == 0     0.594399   0.566664   1.049   1.0000    
# 9 - 2 == 0     0.683907   0.570638   1.198   0.9998    
# 10 - 2 == 0   -0.826305   0.569009  -1.452   0.9966    
# 4 - 3 == 0    -1.041470   0.639857  -1.628   0.9862    
# 5 - 3 == 0    -0.468271   0.597605  -0.784   1.0000    
# 6 - 3 == 0    -1.021385   0.565604  -1.806   0.9580    
# 7 - 3 == 0    -1.080831   0.586148  -1.844   0.9491    
# 8 - 3 == 0    -0.985576   0.543546  -1.813   0.9570    
# 9 - 3 == 0    -0.896068   0.549785  -1.630   0.9862    
# 10 - 3 == 0   -2.406280   0.554107  -4.343    <0.01 ** 
#   5 - 4 == 0     0.573199   0.499755   1.147   0.9999    
# 6 - 4 == 0     0.020086   0.460240   0.044   1.0000    
# 7 - 4 == 0    -0.039361   0.482689  -0.082   1.0000    
# 8 - 4 == 0     0.055894   0.434626   0.129   1.0000    
# 9 - 4 == 0     0.145403   0.436389   0.333   1.0000    
# 10 - 4 == 0   -1.364810   0.430709  -3.169   0.1379    
# 6 - 5 == 0    -0.553113   0.397468  -1.392   0.9980    
# 7 - 5 == 0    -0.612560   0.422140  -1.451   0.9966    
# 8 - 5 == 0    -0.517305   0.371207  -1.394   0.9980    
# 9 - 5 == 0    -0.427797   0.372888  -1.147   0.9999    
# 10 - 5 == 0   -1.938009   0.375222  -5.165    <0.01 ***
#   7 - 6 == 0    -0.059447   0.376168  -0.158   1.0000    
# 8 - 6 == 0     0.035808   0.318429   0.112   1.0000    
# 9 - 6 == 0     0.125317   0.316829   0.396   1.0000    
# 10 - 6 == 0   -1.384896   0.318686  -4.346    <0.01 ** 
#   8 - 7 == 0     0.095255   0.344125   0.277   1.0000    
# 9 - 7 == 0     0.184763   0.346947   0.533   1.0000    
# 10 - 7 == 0   -1.325449   0.344772  -3.844   0.0156 *  
#   9 - 8 == 0     0.089509   0.282184   0.317   1.0000    
# 10 - 8 == 0   -1.420704   0.270584  -5.251    <0.01 ***
#   10 - 9 == 0   -1.510212   0.255634  -5.908    <0.01 ***
#   ---
#   Signif. codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1
# (Adjusted p values reported -- single-step method)


# result: 0 to 7 not significantly different; 8 and 9 not significantly different from each other; 10 is significantly different from all other separately -> exactly same as Armstrong and Davenport; but not exactly same as Joe reports in his paper; he said he created a new scale from 0-5 collapsed into one, then a new category for each: 6,7,8,9,10. This would only make sense if 6,7,8,9,10 were all significantly different from each other whih they are not according to this result;









# # # # # # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # #
# If I can't get it to work with ologit, do OLS and try glht function to test for equality of coefficients there to get an approximation 
# # # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # ## # # #

# trying out glht with lm to see if it works http://www.ats.ucla.edu/stat/r/faq/testing_contrasts.htm
m2 <- lm(physintflip ~ polity2 + gdpenl + jobinsecurity_bdm +
           +lpop + diss_lag + cowwarl + lmtnest + physint_lag,data=mydata)
l2 <- glht(m2, linfct = mcp(polity2 = "Tukey"))
summary(l2)

summary(glht(res1)) #this works but shows if coefficients are nonzero (not equality between coefficients)

library(testthat)
var <- names(res2$coef)
summary(glht(res2),var))

#http://stats.stackexchange.com/questions/59085/how-to-test-for-simultaneous-equality-of-choosen-coefficients-in-logit-or-probit
var.mat <- vcov(res1)[9:10,9:10]


############ old code that didn't work ##############


# replication with zelig (http://stackoverflow.com/questions/16498849/logistic-regression-with-robust-clustered-standard-errors-in-r)
install.packages("Zelig")
library(Zelig)
mydata <- mydata[complete.cases(mydata),]
res <-zelig(factor(physintflip) ~ polity2 + gdpenl + jobinsecurity_bdm +
              +lpop + diss_lag + cowwarl + lmtnest + physint_lag,model="ologit",robust=TRUE,cluster="cowcode",data=mydata)
summary(res) #NO match with Young 2009  ; I sent an email to googlegroup for Zelig; what I did wrong was using "logit" instead of "ologit" earlier; when I use "ologit" then Zelig calls polr() and provides the exact same results as polr() -- meaning that it cannot give the country clustered se's -- I askeda a follow up question on googlegroup for Zelig; 
install.packages("Zelig.Choice") #not working with my current version of R

# try using older Zelig version where according to internet clustering still worked (google group https://groups.google.com/forum/#!topic/zelig-statistical-software/WiHx2Lp5RrE)

#Â with polr
res2 <- polr(factor(physintflip) ~ polity2 + gdpenl + jobinsecurity_bdm +
               +lpop + diss_lag + cowwarl + lmtnest + physint_lag, data=mydata, Hess = TRUE)
summary(res2)
coeftest(res2) #is quite similar to Young but not as good as lrm() and robcov()
pR2(res2) # works:  0.3035340 , but sometimes does not work 0.6242868 
