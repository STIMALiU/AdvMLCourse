# PGMs

library(bnlearn)
library(gRain)

MH<-model2network("[D][P][Ma|D:P][Mb|D:P]") # Structure
cptD = matrix(c(1/4, 1/4, 1/4, 1/4), ncol = 4, dimnames = list(NULL, c("D1", "D2", "D3", "D4"))) # Parameters
cptP = matrix(c(1/4, 1/4, 1/4, 1/4), ncol = 4, dimnames = list(NULL, c("P1", "P2", "P3", "P4")))
cptMa = c(
  0,1/3,1/3,1/3,
  0,0,1/2,1/2,
  0,1/2,0,1/2,
  0,1/2,1/2,0,
  0,0,1/2,1/2,
  1/3,0,1/3,1/3,
  1/2,0,0,1/2,
  1/2,0,1/2,0,
  0,1/2,0,1/2,
  1/2,0,0,1/2,
  1/3,1/3,0,1/3,
  1/2,1/2,0,0,
  0,1/2,1/2,0,
  1/2,0,1/2,0,
  1/2,1/2,0,0,
  1/3,1/3,1/3,0)
dim(cptMa) = c(4, 4, 4)
dimnames(cptMa) = list("Ma" = c("M1", "M2", "M3", "M4"), "D" =  c("D1", "D2", "D3", "D4"), "P" = c("P1", "P2", "P3", "P4"))
cptMb = c(
  0,1/3,1/3,1/3,
  0,0,1/2,1/2,
  0,1/2,0,1/2,
  0,1/2,1/2,0,
  0,0,1/2,1/2,
  1/3,0,1/3,1/3,
  1/2,0,0,1/2,
  1/2,0,1/2,0,
  0,1/2,0,1/2,
  1/2,0,0,1/2,
  1/3,1/3,0,1/3,
  1/2,1/2,0,0,
  0,1/2,1/2,0,
  1/2,0,1/2,0,
  1/2,1/2,0,0,
  1/3,1/3,1/3,0)
dim(cptMb) = c(4, 4, 4)
dimnames(cptMb) = list("Mb" = c("M1", "M2", "M3", "M4"), "D" =  c("D1", "D2", "D3", "D4"), "P" = c("P1", "P2", "P3", "P4"))
MHfit<-custom.fit(MH,list(D=cptD,P=cptP,Ma=cptMa,Mb=cptMb))
MHcom<-compile(as.grain(MHfit))

MHfitEv<-setFinding(MHcom,nodes=c(""),states=c("")) # Exact inference
querygrain(MHfitEv,c("P"))
MHfitEv<-setFinding(MHcom,nodes=c("D","Ma"),states=c("D1","M2"))
querygrain(MHfitEv,c("P"))
MHfitEv<-setFinding(MHcom,nodes=c("D","Ma","Mb"),states=c("D1","M2","M2"))
querygrain(MHfitEv,c("P"))
MHfitEv<-setFinding(MHcom,nodes=c("D","Ma","Mb"),states=c("D1","M2","M3"))
querygrain(MHfitEv,c("P"))

# HMMs

library(HMM)

TimeSteps=100
States=1:10 # 1:5 correspond to mode 1 and 6:10 to mode 2
Symbols=1:5
transProbs=matrix(c(.4,.4,0,0,0,.1,.1,0,0,0,
                    0,.4,.4,0,0,0,.1,.1,0,0,
                    0,0,.4,.4,0,0,0,.1,.1,0,
                    0,0,0,.4,.4,0,0,0,.1,.1,
                    .4,0,0,0,.4,.1,0,0,0,.1,
                    .06,.14,0,0,0,.24,.56,0,0,0,
                    0,.06,.14,0,0,0,.24,.56,0,0,
                    0,0,.06,.14,0,0,0,.24,.56,0,
                    0,0,0,.06,.14,0,0,0,.24,.56,
                    .14,0,0,0,.06,.56,0,0,0,.24), nrow=length(States), ncol=length(States), byrow = TRUE)
emissionProbs=matrix(c(1/3,1/3,0,0,1/3,
                       1/3,1/3,1/3,0,0,
                       0,1/3,1/3,1/3,0,
                       0,0,1/3,1/3,1/3,
                       1/3,0,0,1/3,1/3,
                       1/3,1/3,0,0,1/3,
                       1/3,1/3,1/3,0,0,
                       0,1/3,1/3,1/3,0,
                       0,0,1/3,1/3,1/3,
                       1/3,0,0,1/3,1/3), nrow=length(States), ncol=length(Symbols), byrow = TRUE)
startProbs=c(.2,.2,.2,.2,.2,0,0,0,0,0)
hmm=initHMM(States,Symbols,startProbs,transProbs,emissionProbs)

sim=simHMM(hmm,TimeSteps)
post=posterior(hmm,sim$observation[1:TimeSteps])
sum(post[6:10,TimeSteps])

# RL

# GPs

library(kernlab)
CWData <- read.table('https://raw.githubusercontent.com/STIMALiU/AdvMLCourse/master/GaussianProcess/Code/CanadianWages.dat', 
                     header = T)
logWage<-CWData$logWage
age<-CWData$age

polyFit <- lm(logWage ~  age + I(age^2) + I(age^3))
sigmaNoise = sd(polyFit$residuals)

ell <- 100
SEkernel <- rbfdot(sigma = 1/(2*ell^2)) # Note the reparametrization.

x<-age
xs<-c(age,seq(max(age)+1,100,1)) # XStar.
n <- length(x)
Kxx <- kernelMatrix(kernel = SEkernel, x = x, y = x)
Kxs <- kernelMatrix(kernel = SEkernel, x = x, y = xs)
Meanf = sin(xs)+t(Kxs)%*%solve(Kxx + sigmaNoise^2*diag(n), logWage-sin(x))

plot(age,logWage,xlim=c(min(age),100),ylim=c(0,max(logWage)))
lines(xs, Meanf, col="red", lwd = 2)

# ell=1 allows for enough non-smoothness to fit the training data. Outside the training
# data, it falls back to the prior. There is an abrupt transition between the two regimes.
# Such an abrupt transition is not allowed by ell=100, since it implies more smoothness.
# In this case, the posterior mean combines the prior and the data, i.e., it has the sine
# form of the prior mean but translated so as to fit the training data.
