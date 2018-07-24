library(devtools)
library(Biobase)
library(broom)

data(ALL)
edata=exprs(ALL)
pdata=pData(ALL)


### linear model with sex(gender)
edata=as.matrix(edata)
lm1=lm(edata[1,]~pdata$sex)
lm1
tidy(lm1)
plot(pdata$age,edata[1,], col=1)
pdata$sex
table(pdata$sex)
boxplot(edata[1,]~pdata$sex)
points(edata[1,]~jitter(as.numeric(pdata$sex)),col=as.numeric(pdata$sex))
par(mfrow=c(2,2))
plot(lm1)

### linear model with age
edata=as.matrix(edata)
lm2=lm(edata[1,]~pdata$age)
lm2
tidy(lm2)
par(mfrow=c(1,1))
plot(pdata$age,edata[1,], col=1)
pdata$age
table(pdata$age)
boxplot(edata[1,]~pdata$age)
points(edata[1,]~jitter(as.numeric(pdata$age)),col=as.numeric(pdata$age))
par(mfrow=c(2,2))
plot(lm2)


### linear model with age + sex
edata=as.matrix(edata)
lm3=lm(edata[1,]~pdata$sex + pdata$age)
lm3
tidy(lm3)
table(pdata$age, pdata$sex)
par(mfrow=c(2,2))
plot(lm3)


### linear model with age * sex (interaction)
edata=as.matrix(edata)
lm4=lm(edata[1,]~pdata$sex * pdata$age)
lm4
tidy(lm4)
table(pdata$age, pdata$sex)
par(mfrow=c(2,2))
plot(lm4)
