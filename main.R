# Written by Robin Lee

library(dplyr)
library(tidyr)
library(tidyverse)
library(geoR)
library(gstat)
library(sp)
library(maps)
library(corrplot)
library(ggplot2)
library(gridExtra)

data <- read.csv("epa_data.csv")
data <- data %>% 
  group_by(TRI.Facility.ID) %>%
  pivot_wider(names_from = Chemical, values_from = Releases..lb.) %>%
  rename("lead" = "Lead compounds (N420)") %>%
  mutate(Waste.Managed..lb. = as.numeric(Waste.Managed..lb.),
         lead = as.numeric(lead)) %>%
  summarize(longitude = mean(Longitude),
            latitude = mean(Latitude),
            mean_color = mean(People.of.Color.Percentile),
            mean_lowincome = mean(Low.Income.Percentile),
            mean_unemployment = mean(Unemployment.Rate.Percentile),
            mean_waste = mean(Waste.Managed..lb., na.rm=T),
            mean_lead = mean(lead, na.rm=T)) %>%
  filter(!is.na(mean_lead), !is.na(mean_waste))

map("county", "ca", main = "California Map of County Seats")
points(data$longitude, data$latitude)

unemployment <- read.csv("Unemployment.csv") %>% 
  filter(State == "CA", Area_Name != "California") %>%
  spread(Attribute, Value) %>%
  select(Area_Name, Med_HH_Income_Percent_of_State_Total_2021,
         Unemployment_rate_2021) %>%
  rename("PercentStateMedianIncome" = 
           "Med_HH_Income_Percent_of_State_Total_2021") %>%
  mutate(Area_Name = str_replace(Area_Name, " County, CA", "")) %>%
  mutate(Area_Name = if_else(Area_Name == "San Francisco County/city, CA", 
                               "San Francisco", as.character(Area_Name))) 
poverty <- read.csv("PovertyEstimates.csv") %>% 
  filter(Stabr == "CA", Area_name != "California") %>%
  spread(Attribute, Value) %>%
  rename("Area_Name" = "Area_name") %>%
  select(Area_Name, PCTPOVALL_2021) %>%
  mutate(Area_Name = str_replace(Area_Name, " County", ""))
education <- read.csv("Education.csv") %>%
  filter(State == "CA", Area.name != "California") %>%
  spread(Attribute, Value)  %>%
  rename("Area_Name" = "Area.name", 
         "Bachelors_2021" = 
           "Percent of adults with a bachelor's degree or higher, 2017-21",
         "Highschool_2021" =
           "Percent of adults with a high school diploma only, 2017-21") %>%
  select(Area_Name, Bachelors_2021, Highschool_2021) %>%
  mutate(Area_Name = str_replace(Area_Name, " County", ""))
population <- read.csv("PopulationEstimates.csv") %>% 
  filter(State == "CA", Area_Name != "California") %>%
  spread(Attribute, Value) %>% 
  select(Area_Name, R_NET_MIG_2021, R_DEATH_2021) %>%
  mutate(Area_Name = str_replace(Area_Name, " County", ""))
county_loc <- 
  read.table("http://www.stat.ucla.edu/~nchristo/statistics_c173_c273/ca_seats_coord.txt", 
                header=TRUE)  %>%
  rename("Area_Name" = "county")
adj <- 
  read.table("http://www.stat.ucla.edu/~nchristo/statistics_c173_c273/county_adjacency.txt", 
                sep="\t", fill=FALSE, strip.white=TRUE)[,c(1,3)]
data <- unemployment %>% 
  left_join(poverty, by="Area_Name") %>%
  left_join(population, by="Area_Name") %>%
  left_join(education, by="Area_Name") %>%
  mutate(longitude = county_loc$longitude,
         latitude = county_loc$latitude) %>%
  select(longitude, latitude, Unemployment_rate_2021, PercentStateMedianIncome, 
         PCTPOVALL_2021, Bachelors_2021, Highschool_2021, R_NET_MIG_2021, R_DEATH_2021)
head(data)

par(mfrow=c(1,4))
map("county", "ca", main = "California Map of County Seats", xlim = c(-124.3,-114.2))
#title("California County Seat Locations", cex.main = 1)
points(data$longitude, data$latitude)

# Green: Low, Orange: Moderate, Red: High
colors <- c("green", "orange", "red")
# Low: 0-6, Moderate: 6-10, High: 10+
levels <- as.numeric(cut(data$Unemployment_rate_2021, 
                         c(0, 6, 10, max(data$Unemployment_rate_2021))))

map("county", "ca", main = "Map of County Seats", xlim = c(-124.3,-114.2))
title("Bubble Plot of Unemployment Rate", cex.main = 1)
points(data$longitude, data$latitude, 
       cex=(data$Unemployment_rate_2021/mean(data$Unemployment_rate_2021)), 
       col=colors[levels], pch=19)


par(mfrow=c(1,3))
hist(data$Unemployment_rate_2021,main="County Unemployment Rate",xlab ="Percent")
hist(data$PercentStateMedianIncome, 
     main = "County Median Income against State Median", xlab = "Percent")
hist(data$PCTPOVALL_2021,main ="Percent of Individuals below Poverty",xlab ="Percent")

hist(data$Bachelors_2021,main="Percent with Bachelor's or Higher",xlab ="Percent")
hist(data$Highschool_2021,main ="Percent with Only High School Diploma",xlab="Percent")
hist(data$R_NET_MIG_2021,main = "Rate of County Migration",xlab = "Percent")

pairs(~ Unemployment_rate_2021 + PercentStateMedianIncome + PCTPOVALL_2021 + 
        Bachelors_2021 + latitude + longitude, data=data)

data_h <- data
coordinates(data_h) <- ~longitude+latitude
qq <- hscat((Unemployment_rate_2021)~1, data=data_h, c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1))
plot(qq, main="h-scatterplots")

correlation_matrix <- cor(data[, 3:9])
corrplot(correlation_matrix, method = "color")

par(mfrow=c(1,3))
boxplot(data$Unemployment_rate_2021,main="Unemployment Rate")
boxplot(data$PercentStateMedianIncome,main="Percent of State Median Income")
boxplot(data$Bachelors_2021,main="Percent with Bachelor's Degree")

par(mfrow=c(1,3))
plot(ecdf(data$Unemployment_rate_2021), main = "Unemployment Rate")
plot(ecdf(data$PCTPOVALL_2021), main = "Poverty")
plot(ecdf(data$Highschool_2021), main = "Percent with High School Only")

par(mfrow=c(1,4))
map("county", "ca", main = "Map of County Seats")
points(county_loc$longitude, county_loc$latitude)
text(county_loc$longitude, county_loc$latitude-0.1, labels=county_loc$Area_Name)

map("county", "ca", main = "Map of County Seats", xlim=c(-123.4,-121),ylim=c(36,39))
points(county_loc$longitude, county_loc$latitude)
text(county_loc$longitude, county_loc$latitude-0.1, labels=county_loc$Area_Name)

adj <- adj %>% mutate(V1 = replace(V1, V1 == "", NA))
adj1 <- adj %>% fill(V1,.direction = "down")
adj1 <- adj1 %>% filter(grepl("CA", V1, useBytes = TRUE)) %>%
  filter(grepl("CA", V3, useBytes = TRUE)) %>% 
  mutate(V1 = replace(V1, V1==V3, NA))
adj_df <- as.matrix(table(adj1$V1, adj1$V3))
head(adj_df,c(6,3))

ones <- rep(1, nrow(data))
X <- as.matrix(cbind(ones, data$PercentStateMedianIncome, data$PCTPOVALL_2021, 
                     data$Bachelors_2021, data$Highschool_2021, data$R_NET_MIG_2021,
                     data$R_DEATH_2021, data$longitude))
Q <- t(X) %*% X
beta_hat <- solve(t(X) %*% X) %*% t(X) %*% data$Unemployment_rate_2021
Yhat <- X %*% beta_hat
e <- data$Unemployment_rate_2021 - Yhat
H <- X %*% solve(t(X) %*% X) %*% t(X)
n <- 58
ones <- rep(1, nrow(adj_df))
S0 <- t(ones) %*% adj_df %*% ones
I <- (n / S0) * (t(e) %*% adj_df %*% e) / (t(e) %*% e)
I

k <- 7
n <- 58
E_I <- -n/S0 * (sum(diag((adj_df %*% H))))/(n-k-1)
E_I

data_unemployment <- data %>% select(longitude, latitude, Unemployment_rate_2021)
b <- as.geodata(data_unemployment)

var1 <- variog(b, max.dist=5)
var2 <- variog(b, max.dist=5, estimator.type="modulus")

plot(var1, main = "empirical semivariogram")
points(var2$u, var2$v, col="blue")
legend("topleft", legend = c("Classical Estimator", "Robust Estimator"), 
       col = c("black", "blue"), pch=21)

par(mfrow=c(1,3))
plot(var1, main = "Default Spherical Semivariogram"); points(var2$u, var2$v, col="blue")
default <- variofit(var1,cov.model="sph",ini.cov.pars=c(3,1.5),fix.nugget=FALSE,nugget=0)
lines(default, lty=1, col = "magenta")
plot(var1, main = "Cressie Spherical Semivariogram"); points(var2$u, var2$v, col="blue")
cressie <- variofit(var1, cov.model="sph", weights="cressie", ini.cov.pars=c(3,1.5),
       fix.nugget=FALSE, nugget=0)
lines(cressie, lty=1, col="red")
plot(var1, main = "Equal Spherical Semivariogram"); points(var2$u, var2$v, col="blue")
equal <- variofit(var1, cov.model="sph", ini.cov.pars=c(3,1.5), weights="equal",
       fix.nugget=FALSE, nugget=0)
lines(equal, lty=1, col="darkgreen")

par(mfrow=c(1,3))
plot(var1, main = "Default Exponential Semivariogram"); points(var2$u, var2$v, col="blue")
default_exp <- variofit(var1,cov.model="exp",ini.cov.pars=c(3,1.5),fix.nugget=FALSE,nugget=0)
lines(default_exp, lty=1, col = "magenta")
plot(var1, main = "Cressie Exponential Semivariogram"); points(var2$u, var2$v, col="blue")
cressie_exp <- variofit(var1, cov.model="exp", weights="cressie", ini.cov.pars=c(3,1.5),
       fix.nugget=FALSE, nugget=0)
lines(cressie_exp, lty=1, col="red")
plot(var1, main = "Equal Exponential Semivariogram"); points(var2$u, var2$v, col="blue")
equal_exp <- variofit(var1, cov.model="exp", ini.cov.pars=c(3,1.5), weights="equal",
       fix.nugget=FALSE, nugget=0)
lines(equal_exp, lty=1, col="darkgreen")

## Fit the spherical variogram to the sample variogram 
fit1 <- variofit(var1, cov.model="sph",ini.cov.pars=c(3,1.5),fix.nugget=FALSE,nugget=0)
x_val1 <- xvalid(b, model=fit1,reest=TRUE,variog.obj = var1)
press1 <- sum(x_val1$error^2)
press1/58

## Fit the exponential variogram to the sample variogram 
fit2 <- variofit(var1, cov.model="exp",ini.cov.pars=c(3,1.5),fix.nugget=FALSE,nugget=0)
x_val2 <- xvalid(b, model=fit2,reest=TRUE,variog.obj = var1)
press2 <- sum(x_val2$error^2)
press2/58

g_unemployment <- gstat(id="unemployment_log", formula = log(Unemployment_rate_2021)~1, 
                  locations = ~longitude+latitude, data = data_unemployment)
var_unemployment <- variogram(g_unemployment) 
g_unemployment1 <- gstat(id="unemployment_log", formula = 
                        log(Unemployment_rate_2021)~longitude+latitude, 
                      locations = ~longitude+latitude, data = data_unemployment) 
var_unemployment1 <- variogram(g_unemployment1)
plot1 <- plot(var_unemployment, main = "Empirical Variogram")
plot2 <- plot(var_unemployment1, main = "De-Trended Empirical Variogram") 
grid.arrange(plot1, plot2, nrow = 1)

var_fit1 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",0.5,0),fit.method=1)
plot1 <- plot(var_unemployment1, var_fit1, main = "N_h Weights") 
var_fit2 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",1,0),fit.method=2) 
plot2 <- plot(var_unemployment1, var_fit2, main = "Cressie's Weights")
grid.arrange(plot1, plot2, nrow = 1)

var_fit6 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",1.5,0),fit.method=6) 
plot3 <- plot(var_unemployment1, var_fit6, main = "OLS (No Weights)")
var_fit7 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",1.2,0),fit.method=7) 
plot4 <- plot(var_unemployment1, var_fit7, main = "Default Weights")
grid.arrange(plot3, plot4, nrow = 1)

var_fit1_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",0.5,0),fit.method=1)
plot1 <- plot(var_unemployment1, var_fit1_exp, main = "N_h Weights") 
var_fit2_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",1,0),fit.method=2) 
plot2 <- plot(var_unemployment1, var_fit2_exp, main = "Cressie's Weights")
grid.arrange(plot1, plot2, nrow = 1)

var_fit6_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",1.5, 0),fit.method=6) 
plot3 <- plot(var_unemployment1, var_fit6_exp, main = "OLS (No Weights)")
var_fit7_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",1.2,0),fit.method=7) 
plot4 <- plot(var_unemployment1, var_fit7_exp, main = "Default Weights")
grid.arrange(plot3, plot4, nrow = 1)

#Using spherical: 
var_fit1 <- fit.variogram(variogram(g_unemployment1), vgm(0.045,"Sph",0.5,0), 
                          fit.method=1)
cv_pr <- krige.cv(log(Unemployment_rate_2021)~1,data=data_unemployment, 
                  locations=~longitude+latitude, 
                  model=var_fit1,nfold=nrow(data_unemployment))
#names(cv_pr)
press_cv <- sum(cv_pr$residual^2)
summary(cv_pr)[,1:4]
press_cv

#Using exponential:
var_fit1_exp <- fit.variogram(variogram(g_unemployment1), vgm(0.045,"Exp",0.5,0), 
                          fit.method=1)
cv_pr1 <- krige.cv(log(Unemployment_rate_2021)~1,data=data_unemployment, 
                   locations=~longitude+latitude, 
                   model=var_fit1_exp,nfold=nrow(data_unemployment))
#names(cv_pr1)
press_cv1 <- sum(cv_pr1$residual^2)
summary(cv_pr1)[,1:4]
press_cv1

# Kriging Predictions using Exponential Semivariogram
data_logunemployment <- data %>% 
  select(longitude, latitude, Unemployment_rate_2021) %>%
  mutate(Unemployment_rate_2021 = log(Unemployment_rate_2021)) %>%
  rename("x" = "longitude", "y" = "latitude")
data_unemployment <- data %>% select(longitude, latitude, Unemployment_rate_2021)
b <- as.geodata(data_unemployment)
x.range <- as.integer(range(data_logunemployment[,1])) 
y.range <- as.integer(range(data_logunemployment[,2])) 
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=0.1), 
                   y=seq(from=y.range[1], to=y.range[2], by=0.1)) 

## Ordinary Kriging (gstat)
var_fit1_exp <- fit.variogram(variogram(g_unemployment1), vgm(0.05,"Exp",0.5,0), 
                          fit.method=1)

#Perform ordinary kriging predictions:
pr_ok <- krige(id="logunemployment", Unemployment_rate_2021~1, 
               locations=~x+y, model=var_fit1_exp, 
               data=data_logunemployment, newdata=grd) 

#Collapse the vector of the predicted values into a matrix: 
qqq <- matrix(pr_ok$logunemployment.pred, 
              length(seq(from=x.range[1], to=x.range[2], by=0.1)), 
              length(seq(from=y.range[1], to=y.range[2], by=0.1)) ) 

par(mfrow=c(1,2))
image.orig <- image
image.orig(seq(from=x.range[1],to=x.range[2],by=0.1),
           seq(from=y.range[1],to=y.range[2], by=0.1),qqq,
           xlab="West to East",ylab="South to North", main="Predicted values")
points(data_logunemployment)
#map("county", "ca",add=TRUE)
#Identify the location of each point in the grid:
in.what.state <- map.where(database="state", x=grd$x, y=grd$y)
#Find the points of the grid that belong to California:
in.ca <- which(in.what.state=="california")
pred <- pr_ok$logunemployment.pred
#Assign NA values to all the points outside California:
pred[-in.ca] <- NA
qqq.ca <- matrix(pred, length(seq(from=x.range[1], to=x.range[2], by=0.1)),
                 length(seq(from=y.range[1], to=y.range[2], by=0.1)))
image.orig(seq(from=x.range[1],to=x.range[2],by=0.1),
           seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
           xlab="West to East",ylab="South to North", main="Predicted values")
contour(seq(from=x.range[1], to=x.range[2], by=0.1),
        seq(from=y.range[1],to=y.range[2], by=0.1), qqq.ca, add=TRUE,
        col="black", labcex=1)

filled.contour(seq(from=x.range[1],to=x.range[2],by=0.1),
               seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
               xlab="West to East",ylab="South to North", main="Predicted values",
               col=heat.colors(10))

## Universal Kriging (geoR)

fit2 <- variofit(var2,cov.model="exp",ini.cov.pars=c(3,1.5),fix.nugget=FALSE,nugget=0)
kc <- krige.conv(b, locations=grd, krige=krige.control(obj.model=fit2),
                 nugget=0, trend.l="1st", trend.d="1st")

#univ_kr <- krige(id="logunemployment", Unemployment_rate_2021~x+y, 
               #locations=~x+y, model=var_fit1, 
               #data=data_logunemployment, newdata=grd) 
#qqq <- matrix(univ_kr$logunemployment.pred, 
              #length(seq(from=x.range[1], to=x.range[2], by=0.1)), 
              #length(seq(from=y.range[1], to=y.range[2], by=0.1)) ) 

par(mfrow=c(1,2))
image.orig <- image
image.orig(seq(from=x.range[1],to=x.range[2],by=0.1),seq(from=y.range[1],to=y.range[2],by=0.1),
           qqq,xlab="West to East",ylab="South to North", main="Predicted values")
points(data_unemployment)
in.what.state <- map.where(database="state", x=grd$x, y=grd$y)
in.ca <- which(in.what.state=="california")
pred <- kc$predict; pred[-in.ca] <- NA
qqq.ca <- matrix(pred, length(seq(from=x.range[1], to=x.range[2], by=0.1)),
                 length(seq(from=y.range[1], to=y.range[2], by=0.1)))
image.orig(seq(from=x.range[1],to=x.range[2],by=0.1),
           seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
           xlab="West to East",ylab="South to North", main="Predicted values")
contour(seq(from=x.range[1], to=x.range[2], by=0.1),
        seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,add=TRUE,col="black",labcex=1)

filled.contour(seq(from=x.range[1],to=x.range[2],by=0.1),
               seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
               xlab="West to East",ylab="South to North", main="Predicted values",
               col=heat.colors(10))

var_fit1_exp <- fit.variogram(variogram(g_unemployment1), vgm(0.05,"Exp",0.5,0), 
                          fit.method=1)

data_cokrig <- data %>%
  select(longitude, latitude, Unemployment_rate_2021,
         PercentStateMedianIncome, PCTPOVALL_2021, Bachelors_2021,Highschool_2021) %>% 
  mutate(Unemployment_rate_2021 = log(Unemployment_rate_2021),
         PercentStateMedianIncome = log(PercentStateMedianIncome),
         PCTPOVALL_2021 = log(PCTPOVALL_2021),
         Bachelors_2021 = log(Bachelors_2021),
         Highschool_2021 = log(Highschool_2021)) %>% 
  rename("x" = "longitude", "y" = "latitude")
x.range <- as.integer(range(data_cokrig[,1])) 
y.range <- as.integer(range(data_cokrig[,2])) 
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=0.1), 
                   y=seq(from=y.range[1], to=y.range[2], by=0.1)) 
#target variable
g1 <- gstat(id="log_unemployment", formula = (Unemployment_rate_2021)~1, 
            locations = ~x+y, data = data_cokrig) 
g1 <- gstat(g1,id="log_percentmedian", formula = (PercentStateMedianIncome)~1, 
            locations = ~x+y, data = data_cokrig)
g1 <- gstat(g1,id="log_poverty", formula = (PCTPOVALL_2021)~1, 
            locations = ~x+y, data = data_cokrig)
g1 <- gstat(g1,id="log_bachelors", formula = (Bachelors_2021)~1, 
            locations = ~x+y, data = data_cokrig)
g1 <- gstat(g1,id="log_highschool", formula = (Highschool_2021)~1, 
            locations = ~x+y, data = data_cokrig)
vm <- variogram(g1) 
vm.fit_exp <- fit.lmc(vm, g1, model=var_fit1_exp,correct.diagonal=1.01)
ck <- predict(vm.fit_exp, grd)
qqq <- matrix(ck$log_unemployment.pred,
              length(seq(from=x.range[1], to=x.range[2], by=0.1)),
              length(seq(from=y.range[1], to=y.range[2], by=0.1)))

par(mfrow=c(1,2))
image.orig <- image
image.orig(seq(from=x.range[1],to=x.range[2],by=0.1),
           seq(from=y.range[1],to=y.range[2], by=0.1),qqq,
           xlab="West to East",ylab="South to North", main="Predicted values")
points(data_cokrig)
in.what.state <- map.where(database="state", x=grd$x, y=grd$y)
in.ca <- which(in.what.state=="california")
pred <- ck$log_unemployment.pred
pred[-in.ca] <- NA
qqq.ca <- matrix(pred, length(seq(from=x.range[1], to=x.range[2], by=0.1)),
                 length(seq(from=y.range[1], to=y.range[2], by=0.1)))
image.orig(seq(from=x.range[1],to=x.range[2],by=0.1),
           seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
           xlab="West to East",ylab="South to North", main="Predicted values")
contour(seq(from=x.range[1], to=x.range[2], by=0.1),
        seq(from=y.range[1],to=y.range[2], by=0.1), qqq.ca, add=TRUE,
        col="black", labcex=1)

filled.contour(seq(from=x.range[1],to=x.range[2],by=0.1),
               seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
               xlab="West to East",ylab="South to North", main="Predicted values",
               col=heat.colors(10))

var_fit1 <- fit.variogram(variogram(g_unemployment1), vgm(0.045,"Sph",0.5,0), 
                          fit.method=1)

pr_ok <- krige.cv(Unemployment_rate_2021~1, data=data_logunemployment,
               locations=~x+y, model=var_fit1, nfold = nrow(data_logunemployment)) 
PRESS_ok <- sum(pr_ok$residual^2) / nrow(data_logunemployment)

PRESS_ok

pr_uk <- krige.cv(Unemployment_rate_2021~x+y, data=data_logunemployment,
               locations=~x+y, model=var_fit1, nfold = nrow(data_logunemployment)) 
PRESS_uk <- sum(pr_uk$residual^2) / nrow(data_logunemployment)

PRESS_uk

var_fit1 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",0.4,0),fit.method=1)
vm <- variogram(g1) 
vm.fit_sph <- fit.lmc(vm, g1, model=var_fit1,correct.diagonal=1.01)
cv_ck_sph <- gstat.cv(vm.fit_sph)
PRESS_ck_sph <- sum(cv_ck_sph$residual^2)/nrow(data_cokrig)

vm.fit_sph <- fit.lmc(vm, g1, model=var_fit1,correct.diagonal=1.01)
cv_ck_sph <- gstat.cv(vm.fit_sph)
PRESS_ck_sph <- sum(cv_ck_sph$residual^2)/nrow(data_cokrig)

PRESS_ck_sph

var_fit1_exp <- fit.variogram(variogram(g_unemployment1), vgm(0.05,"Exp",0.5,0), 
                          fit.method=1)

#Compute PRESS for universal kriging with the exponential semivariogram
pr_ok_exp <- krige.cv(Unemployment_rate_2021~1, data=data_logunemployment,
               locations=~x+y, model=var_fit1_exp, nfold = nrow(data_logunemployment)) 
PRESS_ok_exp <- sum(pr_ok_exp$residual^2) / nrow(data_logunemployment)

PRESS_ok_exp

pr_uk_exp <- krige.cv(Unemployment_rate_2021~x+y, data=data_logunemployment,
               locations=~x+y, model=var_fit1_exp, nfold = nrow(data_logunemployment)) 
PRESS_uk_exp <- sum(pr_uk_exp$residual^2) / nrow(data_logunemployment)

PRESS_uk_exp

var_fit1_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",0.5,0),fit.method=1)
vm <- variogram(g1) 
vm.fit_exp <- fit.lmc(vm, g1, model=var_fit1_exp,correct.diagonal=1.01)
cv_ck_exp <- gstat.cv(vm.fit_exp)
PRESS_ck_exp <- sum(cv_ck_exp$residual^2)/nrow(data_cokrig)

vm.fit_exp <- fit.lmc(vm, g1, model=var_fit1_exp,correct.diagonal=1.01)
cv_ck_exp <- gstat.cv(vm.fit_exp)
PRESS_ck_exp <- sum(cv_ck_exp$residual^2)/nrow(data_cokrig)

PRESS_ck_exp
