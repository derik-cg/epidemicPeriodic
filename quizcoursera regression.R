data<-c(0.18, -1.54, 0.42, 0.95)
weights<-c(2, 1, 3, 1)

sum(weights*(data-0.0025)^2)

sum(weights*(data-1.077)^2)

sum(weights*(data-0.1471)^2)

sum(weights*(data-0.3)^2)

#2

x <- c(0.8, 0.47, 0.51, 0.73, 0.36, 0.58, 0.57, 0.85, 0.44, 0.42)
y <- c(1.39, 0.72, 1.55, 0.48, 1.19, -1.59, 1.23, -0.65, 1.49, 0.05)
reg<-lm(y~x-1)
summary(reg )

#3

data("mtcars")
y<-mtcars$mpg
x<-mtcars$wt
reg<-lm(y~x)
summary(reg)
#6
x <- c(8.58,
       10.46, 9.01, 9.64, 8.86)
(x-mean(x))/sd(x)
#7
x <- c(0.8, 0.47, 0.51, 0.73, 0.36, 0.58, 0.57, 0.85, 0.44, 0.42)
y <- c(1.39, 0.72, 1.55, 0.48, 1.19, -1.59, 1.23, -0.65, 1.49, 0.05)
reg<-lm(y~x)
summary(reg )
#9
x <- c(0.8, 0.47, 0.51, 0.73, 0.36, 0.58, 0.57, 0.85, 0.44, 0.42)
mean(x)


######################3
x <- c(0.61, 0.93, 0.83, 0.35, 0.54, 0.16, 0.91, 0.62, 0.62)
y <- c(0.67, 0.84, 0.6, 0.18, 0.85, 0.47, 1.1, 0.65, 0.36)
reg<-lm(y~x)
summary(reg)
1-pt(2.325,7)
###
summary(reg)$sigma
###
reg<-lm(mpg~wt,data=mtcars)
summary(reg)
(37.2851-5.3445*mean(mtcars$wt))-qt(0.975,28)*0.5591
###
attach(mtcars)
predict(reg,newdata=data.frame(wt=mean(wt)),interval="conf")
####
predict(reg,newdata=data.frame(wt=3),interval="pred")
###
reg<-lm(mpg~I(wt*0.5))
summary(reg)
-10.689-qt(0.975,28)*1.118
####
#model with intercept only
regio<-lm(mpg~1)
srsio<-sum(resid(regio)^2)
srs<-sum(resid(reg)^2)
srs/srsio

