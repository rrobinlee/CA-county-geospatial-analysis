# Geo-Spatial Analysis of County-Level Unemployment Rates and Socio-Economic Factors in California

This project utilizes data from the U.S. Department of Agriculture's Economic Research Service (ERS). The ERS compiles the latest statistics on socioeconomic indicators—like poverty rates, population change, unemployment rates, and education levels—and provides maps and data for U.S. States and counties/county equivalents. 

The data can be accessed at: https://www.ers.usda.gov/data-products/county-level-data-sets/. 

![image](https://github.com/user-attachments/assets/f9416e68-cb79-48c2-9c17-91e8a030b0d0)

## Lattice Data Analysis

**1. Moran’s I statistic**

Because the mean is not constant, we compute the Moran’s $I$ test statistic using the residuals. We construct and calculate the following: $\boldsymbol{X, X'X, \hat{\beta}, H, \hat{Y}, e}$, obtaining: 

$$I = \dfrac{n}{S_0} \dfrac{\boldsymbol{e'we}}{\boldsymbol{e'e}}$$

We compute the Moran's $I$ statistic as follows:
```{r}
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
```
* Because we have a negative Moran's $I$ statistic (`-0.03992056`), this suggests that we have spatial dispersion rather than spatial clustering; in other words, similar values are spread apart between counties. 
* Because the value is extremely close to 0, we are confident that there is spatial randomness or no spatial autocorrelation.

**2. Expected Moran's I Statistic**

Given the mean is not constant, we define the expected Moran's $I$ statistic as: 

$$E[I] = -\dfrac{n}{S_0} \dfrac{tr(\boldsymbol{WH})}{n-k-1}$$

where $k$ is the number of predictors, $\boldsymbol{W}$ is the adjacency matrix, and $n$ is the number of counties.
```{r}
k <- 7
n <- 58
E_I <- -n/S0 * (sum(diag((adj_df %*% H))))/(n-k-1)
```
Because our calculated Moran's $I$ statistic (-0.05198904) is neither significantly above or below the expected Moran's $I$ statistic, this suggests that there is no significant spatial autocorrelation in the data as we mentioned above. Therefore, this implies that the distribution of values across space is essentially random, with no discernible clustering or dispersion tendencies in the data. 

## Variogram Calculations using `GeoR`

Spatial statistics computations can be done in `R` using the package `geoR`. To use the `geoR` package we need to convert our data to a `geodata` object. We initialize our data as such:

```{r}
data_unemployment <- data %>% select(longitude, latitude, Unemployment_rate_2021)
b <- as.geodata(data_unemployment)
```

**1. Computing the Empirical Variogram**

The variogram can be computed using the function `variog`. The robust estimator can be used with the argument `estimator.type=“modulus”`, and is robust to outliers compared to the classical estimator. We compute and plot both as such:
```{r, fig.height=6, fig.width=8, out.height="70%", out.width="70%", message=FALSE}
var1 <- variog(b, max.dist=5)
var2 <- variog(b, max.dist=5, estimator.type="modulus")

plot(var1, main = "empirical semivariogram")
points(var2$u, var2$v, col="blue")
legend("topleft", legend = c("Classical Estimator", "Robust Estimator"), 
       col = c("black", "blue"), pch=21)
```
![image](https://github.com/user-attachments/assets/28a3ac7f-2ee3-491c-b876-26ac16635cdb)

Because the variogram function is increasing but does not appear to reach a saturation point, we are confident that the data does not have strong spatial correlation. This is consistent with our prior analysis. 

#### 4.1.2 Fitting Spherical Semivariograms to the Empirical Variogram

Having computed the variogram, we fit a function to it to compute what the variogram graph would look like if we had the entire population of all possible pairs. As such, we select the spherical model and compute the semivariogram using the following three weights: Default (npairs), Cressie, and Equal. 
```{r, fig.height=4, fig.width=10, message=FALSE}
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
```
![image](https://github.com/user-attachments/assets/bedf134e-9784-4b5e-a28e-2987e668b378)

Because the variogram function is increasing but does not appear to reach a saturation point, we are confident that the data does not have strong spatial correlation. This is consistent with our prior analysis. 

#### 4.1.3 Fitting Spherical Semivariograms to the Empirical Variogram

Having computed the variogram, we fit a function to it to compute what the variogram graph would look like if we had the entire population of all possible pairs. As such, we select the spherical model and compute the semivariogram using the following three weights: Default (npairs), Cressie, and Equal. 
```{r, fig.height=4, fig.width=10, message=FALSE}
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
```

![image](https://github.com/user-attachments/assets/b5bf8b27-a135-448a-821f-0b972f1cc3ae)


Using the spherical semivariogram, the default weights seem to fit the data the best. 

#### 4.1.4 Fitting Exponential Semivariograms to the Empirical Variogram

Having computed the variogram, we fit a function to it to compute what the variogram graph would look like if we had the entire population of all possible pairs. As such, we select the exponential model and compute the semivariogram using the following three weights: Default (npairs), Cressie, and Equal. 

```{r, fig.height=4, fig.width=10, message=FALSE}
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
```

![image](https://github.com/user-attachments/assets/a1c96c5f-4e73-44df-bc8e-8ca73a04b5a1)

For the exponential semivariogram, default weights appear to fit the data the best.

#### 4.1.5 PRESS Calculation between Spherical and Exponential

Because we are using the `geoR` package, we can implement the `xvalid` function with `reest=TRUE` in order to compute the PRESS. As such, we then compare the PRESS for the two variograms in order to identify whether the spherical or exponential model is a better fit. **We will use the default weights (npairs).**

```{r}
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
```

Using the `geoR` package and the default weights, the **exponential model has a lower PRESS value**, so we will select this model to make predictions. We will also compare the PRESS for exponential and spherical using the `gstat` package. The `gstat` package utilizes a different method for spatial statistics computations. Thus, we will re-compute the empirical variogram, and fit both spherical and exponential semivariograms. The results should be similar, but will not be exactly the same. 

### 4.2 Variogram Calculations using `Gstat`

Spatial statistics computations can also be done in `R` using the package `Gstat`. Unlike the `geoR` package, we do not need to convert our data to a `geodata` object. We continue using the `data_unemployment` dataset. We expect that the default weight and exponential semivariogram will best fit the data for both packages.

#### 4.2.1 Computing the Empirical Variogram

We will have to take trend into account when computing the variogram. We can fit a linear surface to the data by regressing the data against the Longitude and Latitude coordinates.

```{r, fig.height=4, fig.width=10}
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
```

![image](https://github.com/user-attachments/assets/c6ea14ab-1612-4fd1-b1d1-982167b1ce88)

By de-trending the data, we can see that the variogram appears to reach a saturation point. However, there does not appear to be strong spatial correlation. This is consistent with our prior analysis. 

In order to fit a semivariogram to our empirical variogram above, we can apply certain weights to better fit the data. The default is $\frac{N_h}{h^2}$ where $N_h$ is the number of pairs and $h$ is the separation distance. The other semivariogram weights are as follows:

![image](https://github.com/user-attachments/assets/510e5dc0-91f5-4e56-a7f2-fe1eb8314d10)

We will test these weights with both the spherical and exponential semi-variograms. 

#### 4.2.2 Fitting Spherical Semivariograms to the Empirical Variogram

Having computed our sample empirical variogram, we fit a function to it to compute what the variogram graph would look like if we had the entire population of all possible pairs. As such, we first select the spherical semivariogram with the four different weights: 

```{r, fig.height=4, fig.width=10, warning=FALSE}
var_fit1 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",0.5,0),fit.method=1)
plot1 <- plot(var_unemployment1, var_fit1, main = "N_h Weights") 
var_fit2 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",1,0),fit.method=2) 
plot2 <- plot(var_unemployment1, var_fit2, main = "Cressie's Weights")
grid.arrange(plot1, plot2, nrow = 1)
```

![image](https://github.com/user-attachments/assets/67af46e1-d3bc-47a2-95db-df07003b58cb)

```{r, fig.height=4, fig.width=10, warning=FALSE}
var_fit6 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",1.5,0),fit.method=6) 
plot3 <- plot(var_unemployment1, var_fit6, main = "OLS (No Weights)")
var_fit7 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",1.2,0),fit.method=7) 
plot4 <- plot(var_unemployment1, var_fit7, main = "Default Weights")
grid.arrange(plot3, plot4, nrow = 1)
```

![image](https://github.com/user-attachments/assets/b82f03d6-f6f9-41a9-92fe-acbe8099241b)

For the spherical semivariograms, the default weight ($N_h$) appears to fit the data the best, so we will use this weight for our future analysis. This is the same conclusion as when using the `geoR` package. 

#### 4.2.3 Fitting Exponential Semivariograms to the Empirical Variogram

Having computed the sample empirical variogram, we fit a function to it to compute what the variogram graph would look like if we had the entire population of all possible pairs. As such, we compute the exponential semivariogram with the four weights.
```{r, fig.height=4, fig.width=10, warning=FALSE}
var_fit1_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",0.5,0),fit.method=1)
plot1 <- plot(var_unemployment1, var_fit1_exp, main = "N_h Weights") 
var_fit2_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",1,0),fit.method=2) 
plot2 <- plot(var_unemployment1, var_fit2_exp, main = "Cressie's Weights")
grid.arrange(plot1, plot2, nrow = 1)
```

![image](https://github.com/user-attachments/assets/f13ab046-965f-4952-ae51-11b1eb3bdd04)

```{r, fig.height=4, fig.width=10, warning=FALSE}
var_fit6_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",1.5, 0),fit.method=6) 
plot3 <- plot(var_unemployment1, var_fit6_exp, main = "OLS (No Weights)")
var_fit7_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",1.2,0),fit.method=7) 
plot4 <- plot(var_unemployment1, var_fit7_exp, main = "Default Weights")
grid.arrange(plot3, plot4, nrow = 1)
```

![image](https://github.com/user-attachments/assets/b89bb05d-5477-470c-b08a-98ed18ca838a)

For the exponential semivariograms, the default weight ($N_h$) still appears to fit the data the best, so we will use this weight for our future analysis. This is the same conclusion as when using the `geoR` package. 

#### 4.2.4 PRESS Calculation between Spherical and Exponential

Because we are using the `gstat` package, we can use the `krige.cv` function to conduct leave one out cross validation. The function automatically deletes one point at a time and uses the remaining n-1 points to predict it. We then compare the PRESS for the two variograms in order to identify whether the spherical or exponential model is a better fit. **We will use the default weights (npairs).**

```{r}
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
```

Using the `gstat` package and the default weights, the **exponential model has a lower PRESS value**, so we will select this model to make predictions. This is the same result as with the `geoR` package.  

## 5. Kriging Predictions using Exponential Semivariogram

```{r}
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
```

### 5.1 Ordinary Kriging (`gstat`)

Performing ordinary kriging using `gstat` package and exponential semivariogram:

```{r}
#Perform ordinary kriging predictions:
pr_ok <- krige(id="logunemployment", Unemployment_rate_2021~1, 
               locations=~x+y, model=var_fit1_exp, 
               data=data_logunemployment, newdata=grd) 

#Collapse the vector of the predicted values into a matrix: 
qqq <- matrix(pr_ok$logunemployment.pred, 
              length(seq(from=x.range[1], to=x.range[2], by=0.1)), 
              length(seq(from=y.range[1], to=y.range[2], by=0.1)) ) 
```
#### 5.1.1 Raster Map using the Predicted Values
```{r, fig.height=5, fig.width=8, out.height="70%", out.width="70%"}
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
```
![image](https://github.com/user-attachments/assets/ef1e1baa-92e5-4d2f-b4f1-863a8ae0ae79)

```{r, fig.height=5, fig.width=5, out.height="55%", out.width="55%", message=FALSE}
filled.contour(seq(from=x.range[1],to=x.range[2],by=0.1),
               seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
               xlab="West to East",ylab="South to North", main="Predicted values",
               col=heat.colors(10))
```

![image](https://github.com/user-attachments/assets/a4de9da8-a770-4f2b-b154-08b8a49e18d6)

We can see that the Ordinary Kriging raster map is fairly similar to our bubble plot from the Explanatory Data Analysis. Counties within Central California and nearby the Mexico border have a more extreme unemployment rate. We will also analyze the Universal Kriging and Co-Kriging models with the exponential semivariogram. We will use the Ordinary Kriging for our Co-Kriging model.

### 5.2 Universal Kriging (`geoR`)

Performing universal kriging using `geoR` package and exponential semivariogram:

```{r,include=TRUE, eval=TRUE}
fit2 <- variofit(var2,cov.model="exp",ini.cov.pars=c(3,1.5),fix.nugget=FALSE,nugget=0)
kc <- krige.conv(b, locations=grd, krige=krige.control(obj.model=fit2),
                 nugget=0, trend.l="1st", trend.d="1st")
```
#### 5.2.1 Raster Map using the Predicted Values

```{r, fig.height=5, fig.width=8, out.height="70%", out.width="70%"}
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
```

![image](https://github.com/user-attachments/assets/404ef9f0-6fd5-47c9-9678-e3dda0dbcdba)

```{r, fig.height=5, fig.width=5, out.height="55%", out.width="55%", message=FALSE}
filled.contour(seq(from=x.range[1],to=x.range[2],by=0.1),
               seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
               xlab="West to East",ylab="South to North", main="Predicted values",
               col=heat.colors(10))
```

![image](https://github.com/user-attachments/assets/a5d34c40-e72f-4981-8d33-74faec63a59e)

Analyzing the raster map for Universal Kriging with the exponential semivariogram, we can see that there are two primary areas with an extreme unemployment rate: Central California and the border with Mexico. However, unlike with our Ordinary Kriging raster map, this does not fully encapsulate the region, and we do not believe it is as accurate as the Ordinary model. Thus, we will implement our Co-Kriging model using the Ordinary Kriging. This difference may also be due to the fact that we are using the `geoR` package.

### 5.3 Co-Kriging (`gstat`)

Performing co-kriging using `gstat` package and **Ordinary Kriging**:

```{r}
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
```

For our Co-Kriging model, we have decided to implement the following predictors that we believe are associated with our target variable `Unemployment_rate_2021`. We have log-transformed all variables. 

- `log_percentmedian`: the log median income as a percentage of total median income for CA
- `log_poverty`: the log percentage of the total population living below the poverty line
- `log_bachelors`: the log percentage of adults with a bachelor’s degree or higher
- `log_highschool`: the log county-level percentage of adults with a high school diploma only

Having selected these four predictors, we seek to construct a predictive model that best represents the unemployment rate throughout California counties. We believe that this model will be the most accurate, as we incorporate more factors to boost our prediction.

#### 5.3.1 Raster Map using the Predicted Values

```{r, fig.height=5, fig.width=8, out.height="70%", out.width="70%"}
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
```

![image](https://github.com/user-attachments/assets/89275035-19c3-40f4-8d79-3545323877d1)

```{r, fig.height=5, fig.width=5, out.height="55%", out.width="55%", message=FALSE}
filled.contour(seq(from=x.range[1],to=x.range[2],by=0.1),
               seq(from=y.range[1],to=y.range[2], by=0.1),qqq.ca,
               xlab="West to East",ylab="South to North", main="Predicted values",
               col=heat.colors(10))
```

![image](https://github.com/user-attachments/assets/c99935bc-8a8e-497f-b1c2-6a2937728463)

Comparing the three Kriging models—ordinary, universal, and co-kriging—we can see that the co-kriging model best emulates the bubble plot from our explanatory data analysis when using the exponential semivariogram. We can see that areas in Central California and around the Mexico border tend to have a higher unemployment rate. Furthermore, parts of Northern California also have higher rates compared to the rest of the state.

Therefore, we believe that co-kriging is the best predictive model, with ordinary and universal kriging as the second and third best models. We will utilize cross-validation in order to enumerate their accuracy.

## 6 Kriging Cross-Validation

For our final section, we compute the sample variogram and fit the spherical and exponential variograms to it. While we exclusively utilize the exponential variogram in the previous section, we will compute all kriging PRESS values for both to validate our assumption. Then, using leave one out cross-validation (`krige.cv`), we predict the points and compare the prediction sum of squares (PRESS) for each variogram. 

### 6.1 PRESS Calculations for Kriging using Spherical Variogram

We first focus on computing the cross-validation metrics for the following Kriging models: Ordinary, Universal, and Co-Kriging. We will use the `gstat` package and leave-one-out cross-validation (`krige.cv`).  

#### 6.1.1 Ordinary Kriging PRESS for Spherical
```{r}
pr_ok <- krige.cv(Unemployment_rate_2021~1, data=data_logunemployment,
               locations=~x+y, model=var_fit1, nfold = nrow(data_logunemployment)) 
PRESS_ok <- sum(pr_ok$residual^2) / nrow(data_logunemployment)
```

The PRESS for ordinary kriging is: 0.04684231


#### 6.1.2 Universal Kriging PRESS for Spherical
```{r}
pr_uk <- krige.cv(Unemployment_rate_2021~x+y, data=data_logunemployment,
               locations=~x+y, model=var_fit1, nfold = nrow(data_logunemployment)) 
PRESS_uk <- sum(pr_uk$residual^2) / nrow(data_logunemployment)
```

The PRESS for universal kriging is: 0.04780679


#### 6.1.3 Co-Kriging PRESS for Spherical
```{r, message=FALSE,include=FALSE}
var_fit1 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",0.4,0),fit.method=1)
vm <- variogram(g1) 
vm.fit_sph <- fit.lmc(vm, g1, model=var_fit1,correct.diagonal=1.01)
cv_ck_sph <- gstat.cv(vm.fit_sph)
PRESS_ck_sph <- sum(cv_ck_sph$residual^2)/nrow(data_cokrig)
```

```{r, eval=FALSE}
vm.fit_sph <- fit.lmc(vm, g1, model=var_fit1,correct.diagonal=1.01)
cv_ck_sph <- gstat.cv(vm.fit_sph)
PRESS_ck_sph <- sum(cv_ck_sph$residual^2)/nrow(data_cokrig)
```

The PRESS for co-kriging is: 0.04259223

Because it has the lowest PRESS value, we would recommend choosing the co-kriging model specifically when using the spherical semivariogram. However, because we believe that the exponential variogram will be more accurate, so we will need to calculate the PRESS for all three models again when using this model. We expect the PRESS calculations for kriging using the exponential variogram to be lower.






### 6.2 PRESS Calculations for Kriging using Exponential Variogram

Having caluclated the PRESS for Kriging using the spherical variogram, we compute the cross-validation metrics for the following Kriging models using the exponential variogram: Ordinary, Universal, and Co-Kriging. We expect these values to be lower.   

#### 6.2.1 Ordinary Kriging PRESS for Exponential
```{r}
#Compute PRESS for universal kriging with the exponential semivariogram
pr_ok_exp <- krige.cv(Unemployment_rate_2021~1, data=data_logunemployment,
               locations=~x+y, model=var_fit1_exp, nfold = nrow(data_logunemployment)) 
PRESS_ok_exp <- sum(pr_ok_exp$residual^2) / nrow(data_logunemployment)
```

The PRESS for ordinary kriging is: 0.04654859

#### 6.2.2 Universal Kriging PRESS for Exponential
```{r}
pr_uk_exp <- krige.cv(Unemployment_rate_2021~x+y, data=data_logunemployment,
               locations=~x+y, model=var_fit1_exp, nfold = nrow(data_logunemployment)) 
PRESS_uk_exp <- sum(pr_uk_exp$residual^2) / nrow(data_logunemployment)
```

The PRESS for universal kriging is: 0.04709442

#### 6.2.3 Co-Kriging PRESS for Exponential
```{r, message=FALSE,include=FALSE}
var_fit1_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",0.5,0),fit.method=1)
vm <- variogram(g1) 
vm.fit_exp <- fit.lmc(vm, g1, model=var_fit1_exp,correct.diagonal=1.01)
cv_ck_exp <- gstat.cv(vm.fit_exp)
PRESS_ck_exp <- sum(cv_ck_exp$residual^2)/nrow(data_cokrig)
```

```{r, eval=FALSE}
vm.fit_exp <- fit.lmc(vm, g1, model=var_fit1_exp,correct.diagonal=1.01)
cv_ck_exp <- gstat.cv(vm.fit_exp)
PRESS_ck_exp <- sum(cv_ck_exp$residual^2)/nrow(data_cokrig)
```

The PRESS for co-kriging is: 0.03860466

Because it has the lowest PRESS value, we would recommend choosing the co-kriging model when using the exponential semivariogram. For ordinary and universal kriging, the exponential model is better than the spherical model. Furthermore, the exponential co-kriging is better than the spherical one. 




