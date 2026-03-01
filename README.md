# Geo-Spatial Analysis of County-Level Unemployment Rates and Socio-Economic Factors in California

This project utilizes data from the U.S. Department of Agriculture's Economic Research Service (ERS). The ERS compiles the latest statistics on socioeconomic indicators—like poverty rates, population change, unemployment rates, and education levels—and provides maps and data for U.S. States and counties/county equivalents. 

The data can be accessed at: https://www.ers.usda.gov/data-products/county-level-data-sets/. 

<p align="center">
  <img src="https://github.com/user-attachments/assets/f9416e68-cb79-48c2-9c17-91e8a030b0d0" width="600"/>
    </p>

## Lattice Analysis

**1. Moran’s $\mathcal{I}$ statistic**

Because the mean is not constant, I calculate Moran’s $\mathcal{I}$ test statistic using residuals. Leveraging $\boldsymbol{X, X'X, \hat{\beta}, H, \hat{Y}, e}$, the formula is: 

$$\mathcal{I} = \dfrac{n}{S_0} \dfrac{\boldsymbol{e'we}}{\boldsymbol{e'e}}$$

<details><summary>R Code:</summary>
  
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
</details>

Because Moran's $\mathcal{I}$ statistic (`-0.03992056`) is negative and close to 0, I am confident there is no spatial autocorrelation; in other words, similar values are spread apart between counties.

**2. Expected Moran's I Statistic**

The expected Moran's $\mathcal{I}$ statistic is: 

$$E[\mathcal{I}] = -\dfrac{n}{S_0} \dfrac{tr(\boldsymbol{WH})}{n-k-1}$$

where $k$ is the number of predictors, $\boldsymbol{W}$ is the adjacency matrix, and $n$ is the number of counties.

<details><summary>R Code:</summary>
  
```{r}
k <- 7
n <- 58
E_I <- -n/S0 * (sum(diag((adj_df %*% H))))/(n-k-1)
```
</details>
  
Because the Moran's $I$ statistic (`-0.03992056`) is neither significantly above nor below the expected Moran's $I$ statistic (`-0.05198904`), there is no significant spatial autocorrelation in the data. Therefore, the distribution of values across space is random, with no discernible clustering or dispersion tendencies in the data. 

## Variograms and Semivariograms in Gstat

**1. Computing the Empirical Variogram**

To account for trend, I fit a linear surface to the data by regressing the data against the longitude and latitude coordinates.

<details><summary>R Code:</summary>
  
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
</details>

<p align="center">
  <img src="https://github.com/user-attachments/assets/c6ea14ab-1612-4fd1-b1d1-982167b1ce88" width="600"/>
</p>

While the de-trended variogram appears to reach a saturation point, there does not appear to be strong spatial correlation. In order to improve the semivariogram's fit to the empirical variogram, the data is transformed using the following weights:  

<p align="center">
  <img src="https://github.com/user-attachments/assets/e2d87f96-b781-4d46-b89d-89913aac6aef"  width="300"/>
</p>

**2. Fitting Spherical Semivariograms to the Empirical Variogram**

<details><summary>R Code:</summary>

```{r, fig.height=4, fig.width=10, warning=FALSE}
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
```

</details>

<p align="center">
  <img src="https://github.com/user-attachments/assets/67af46e1-d3bc-47a2-95db-df07003b58cb" width="600"/>
  <img src="https://github.com/user-attachments/assets/b82f03d6-f6f9-41a9-92fe-acbe8099241b" width="600"/>
</p>

For the spherical semivariogram, default weights ($N_h$) appear to best fit the data. 

**3. Fitting Exponential Semivariograms to the Empirical Variogram**

<details><summary>R Code:</summary>

```{r, fig.height=4, fig.width=10, warning=FALSE}
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
```
</details>

<p align="center">
  <img src="https://github.com/user-attachments/assets/f13ab046-965f-4952-ae51-11b1eb3bdd04" width="600"/>
  <img src="https://github.com/user-attachments/assets/b89bb05d-5477-470c-b08a-98ed18ca838a" width="600"/>
</p>


For the exponential semivariogram, default weights ($N_h$) also appear to best fit the data. 

**4. PRESS Statistic between Spherical and Exponential Semivariograms**

The `krige.cv` function conducts LOOCV (leave-one-out cross-validation), which automatically deletes one point at a time and uses the remaining $n-1$ points to predict future values. Comparing the PRESS for the two variograms:

<details><summary>R Code:</summary>
  
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
</details>

The exponential model (`2.699818`) has a lower PRESS value than the spherical model (`2.716854`), so we will select this model to make predictions. This is the same result as with the `geoR` package.  

## Kriging Predictions using the Exponential Semivariogram

**1. Ordinary Kriging Predictions using gstat**

Performing ordinary kriging using the `gstat` package and exponential semivariogram:

<details><summary>R Code:</summary>
  
```{r}
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
```
</details>

<p align="center">
  <img src="https://github.com/user-attachments/assets/ef1e1baa-92e5-4d2f-b4f1-863a8ae0ae79" width="600"/>
  <img src="https://github.com/user-attachments/assets/a4de9da8-a770-4f2b-b154-08b8a49e18d6" width="600"/>
</p>

Counties within Central California and nearby the Mexico border have a more extreme unemployment rate. 

**2. Co-Kriging Predictions using gstat**

Performing Co-Kriging using the `gstat` package and previous Ordinary Kriging:

<details><summary>R Code:</summary>
  
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
```
</details>

Using the following log-transformed predictors: 

- `log_percentmedian`: the log median income as a percentage of total median income for CA
- `log_poverty`: the log percentage of the total population living below the poverty line
- `log_bachelors`: the log percentage of adults with a bachelor’s degree or higher
- `log_highschool`: the log county-level percentage of adults with a high school diploma only


<p align="center">
  <img src="https://github.com/user-attachments/assets/89275035-19c3-40f4-8d79-3545323877d1" width="600"/>
  <img src="https://github.com/user-attachments/assets/c99935bc-8a8e-497f-b1c2-6a2937728463" width="600"/>
</p>

Counties within Central California and around the Mexico border tend to have a higher unemployment rate. Furthermore, parts of Northern California also have higher rates compared to the rest of the state.

## Kriging Cross-Validation

Using LOOCV (`krige.cv`) to compare the prediction sum of squares (PRESS) for each variogram:

**1. PRESS Statistics for the Kriging Predictions using the Spherical Variogram**

Calculating the cross-validation metrics for the following Kriging models: Ordinary, Universal, and Co-Kriging:  

<details><summary>R Code:</summary>
  
```{r}
pr_ok <- krige.cv(Unemployment_rate_2021~1, data=data_logunemployment,
               locations=~x+y, model=var_fit1, nfold = nrow(data_logunemployment)) 
PRESS_ok <- sum(pr_ok$residual^2) / nrow(data_logunemployment)

pr_uk <- krige.cv(Unemployment_rate_2021~x+y, data=data_logunemployment,
               locations=~x+y, model=var_fit1, nfold = nrow(data_logunemployment)) 
PRESS_uk <- sum(pr_uk$residual^2) / nrow(data_logunemployment)

var_fit1 <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Sph",0.4,0),fit.method=1)
vm <- variogram(g1) 
vm.fit_sph <- fit.lmc(vm, g1, model=var_fit1,correct.diagonal=1.01)
cv_ck_sph <- gstat.cv(vm.fit_sph)
PRESS_ck_sph <- sum(cv_ck_sph$residual^2)/nrow(data_cokrig)
vm.fit_sph <- fit.lmc(vm, g1, model=var_fit1,correct.diagonal=1.01)
cv_ck_sph <- gstat.cv(vm.fit_sph)
PRESS_ck_sph <- sum(cv_ck_sph$residual^2)/nrow(data_cokrig)
```
</details>
  
* The PRESS for ordinary kriging is: 0.04684231
* The PRESS for universal kriging is: 0.04780679
* **The PRESS for co-kriging is: 0.04259223**

Because it has the lowest PRESS value, the co-kriging model with the spherical semivariogram best fits the data.

**2. PRESS Calculations for the Kriging Predictions using Exponential Variogram**

Calculating the cross-validation metrics for the following Kriging models: Ordinary, Universal, and Co-Kriging:     

<details><summary>R Code:</summary>
  
```{r}
#Compute PRESS for universal kriging with the exponential semivariogram
pr_ok_exp <- krige.cv(Unemployment_rate_2021~1, data=data_logunemployment,
               locations=~x+y, model=var_fit1_exp, nfold = nrow(data_logunemployment)) 
PRESS_ok_exp <- sum(pr_ok_exp$residual^2) / nrow(data_logunemployment)

pr_uk_exp <- krige.cv(Unemployment_rate_2021~x+y, data=data_logunemployment,
               locations=~x+y, model=var_fit1_exp, nfold = nrow(data_logunemployment)) 
PRESS_uk_exp <- sum(pr_uk_exp$residual^2) / nrow(data_logunemployment)

var_fit1_exp <- fit.variogram(variogram(g_unemployment1),vgm(0.045,"Exp",0.5,0),fit.method=1)
vm <- variogram(g1) 
vm.fit_exp <- fit.lmc(vm, g1, model=var_fit1_exp,correct.diagonal=1.01)
cv_ck_exp <- gstat.cv(vm.fit_exp)
PRESS_ck_exp <- sum(cv_ck_exp$residual^2)/nrow(data_cokrig)
vm.fit_exp <- fit.lmc(vm, g1, model=var_fit1_exp,correct.diagonal=1.01)
cv_ck_exp <- gstat.cv(vm.fit_exp)
PRESS_ck_exp <- sum(cv_ck_exp$residual^2)/nrow(data_cokrig)
```
</details>

* The PRESS for ordinary kriging is: 0.04654859
* The PRESS for universal kriging is: 0.04709442
* **The PRESS for co-kriging is: 0.03860466**

Because it has the lowest PRESS value, I would recommend choosing the co-kriging model when using the exponential semivariogram. For ordinary and universal kriging, the exponential model is better than the spherical model. Furthermore, the exponential co-kriging is better than the spherical one. 




