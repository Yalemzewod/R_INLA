
###Spatio-temporal mapping of malaria incidence in Ethiopia, 2019 - 2021###

remove(list=ls())#clear working directory

#Set working directory
my.dir <- paste(getwd(),"/",sep="")

#################Load packages if required########################
#remotes::install_github("tidyverse/tidyverse")#lubridate included in the tidyverse
library(tidyverse, quietly=T) # data wrangling and visualization
library(hablar, quietly = T) #rename
#to read and clean spatial data
library(spdep, quietly = T) 
library(sp, quietly = T) 
library(sf, quietly = T) 
library(spData, quietly = T) 
require(geoR, quietly = T)
require(Matrix, quietly = T)

library(rgdal, quietly=T)
library(rgeos, quietly=T)

#for ploting
library(tmap, quietly = T)
library(viridis, quietly = T)
library(classInt, quietly=T)

##For modeling 
#remove.packages("INLA")
#install.packages("INLA", repos="https://protect-au.mimecast.com/s/WePPCmOxBAH5Aq0GC3A4pk?domain=inla.r-inla-download.org")
#remotes::install_version("INLA", version="22.05.03",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
require(INLA)
#Inla.upgrade( testing=T)
#remotes::install_github("tidyverse/tidyverse")#lubridate included in the tidyverse


#####################PREPARE DATA##########################
#'Load shapefile
nameshp <- system.file("Z:/Yalem/Ethiopia_malaria/ethiopia_megregd_shp/V3/edited_v3.shp",
                       package = "SpatialEpiApp")

map.eth <- readOGR("Z:/Yalem/Ethiopia_malaria/ethiopia_megregd_shp/V3/edited_v3.shp",
                   verbose = TRUE)

map.eth@data$ADM3_EN <- recode(map.eth@data$ADM3_EN,"Hawi Gudina\n" = "Hawi Gudina") #name mismatch with the dataset

map.eth@data$district <- map.eth@data$ADM3_EN #rename `ADM3_EN` to district

admin1_shp <- st_read("Z:/master_geometries/Admin_Units/Country_Specific/Ethiopia/CSA_BoFED/eth_admbnda_adm1_csa_bofed_20201008.shp")

#projection (if required) - this will give us lat and long for each district
proj4string(map.eth)
tx_ll <- spTransform(map.eth, CRS("+proj=longlat +datum=WGS84 +no_defs"))
sf_eth2 <- st_as_sf(tx_ll)


sf_eth<- st_as_sf(map.eth) #convert to sf for easy data management

#' read case data
#' Due to the conflict in the northern Ethiopia no data in Tigray region for 2021
d4 <- readRDS("Z:/Yalem/Ethiopia_malaria/Ahmed_data/tes_case_year.rds") |> 
  rename("year" = "Year") |> 
  mutate(confirmed_all = case_when(region%in%"Tigray" & year%in%2021 &
                                     !district%in%c("Raya Alamata", "Raya Azebo",
                                                                       "Kafta Humera", "Setit Humera",
                                                                       "Tsegede (TG)","Welkait")~NA_real_,
                                   confirmed_all>test_performed~test_performed,
                                   
                                   TRUE~confirmed_all),
         confirmed_u5 = case_when(region%in%"Tigray" & year%in%2021 &
                                    !district%in%c("Raya Alamata", "Raya Azebo",
                                                   "Kafta Humera", "Setit Humera",
                                                   "Tsegede (TG)","Welkait")~NA_real_,
                                  TRUE~confirmed_u5),
          confirmed_5_14 = case_when(region%in%"Tigray" & year%in%2021 &
                                       !district%in%c("Raya Alamata", "Raya Azebo",
                                                      "Kafta Humera", "Setit Humera",
                                                      "Tsegede (TG)","Welkait")~NA_real_,
                                     TRUE~confirmed_5_14),
         confirmed_15 = case_when(region%in%"Tigray" & year%in%2021 &
                                    !district%in%c("Raya Alamata", "Raya Azebo",
                                                   "Kafta Humera", "Setit Humera",
                                                   "Tsegede (TG)","Welkait")~NA_real_,
                                  TRUE~confirmed_15),
         test_performed = case_when(region%in%"Tigray" & year%in%2021 &
                                     !district%in%c("Raya Alamata", "Raya Azebo",
                                                    "Kafta Humera", "Setit Humera",
                                                    "Tsegede (TG)","Welkait")~NA_real_,
                                   TRUE~test_performed))
#' inspect the data
d4 |> 
  glimpse()


d4 <- d4 %>% 
  mutate(r = sum(confirmed_all)/sum(population), # overall rate stratified by year
         ir = (confirmed_all/population)*1000, #rate per 1000 population
         tp = (confirmed_all/test_performed)*100, #test positivity
         Ex = r*population,
         sir = confirmed_all/Ex)
##explore the data before analysis##
#subset the data by years and inspect the data: the dataset doesn't seems to be well behaved
#2019

d19 <- d4 |> 
  filter(year==2019) 
  

d20 <- d4 |> 
  filter(year==2020) #2020

d21 <- d4 |> 
  filter(year==2021) #2021

##Histogram
hist(d19$test_performed, main=NULL)
hist(d20$test_performed, main=NULL)
hist(d21$test_performed, main=NULL)

##Box-plot
boxplot(d19$test_performed, horizontal = TRUE)
boxplot(d20$test_performed, horizontal = TRUE)
boxplot(d21$test_performed, horizontal = TRUE)


#Computing the Moran’s I statistic
#Recall that the Moran’s I value is the slope of the line that best fits the relationship between neighboring malaria cases
#and each polygon’s number of malaria case in the dataset.
I <- moran(eth_eth9$confirmed_all, lw, length(nb2), Szero(lw))[1]
I
#Analytical method
#Performing a hypothesis test
#There are two methods to testing this hypothesis: an analytical method and a Monte Carlo method.
#To run the Moran’s I analysis using the analytical method, use the `moran.test` function.
moran.test(eth_eth9$confirmed_all,lw, alternative="greater")

#' Note that ArcMap adopts this analytical approach to its hypothesis test however,
#'  it implements a two-sided test as opposed to the one-sided test adopted in the above example (i.e. alternative = "greater"). 
#' A two-sided p-value is nothing more than twice the one-sided p-value.
#'  Unfortunately, ArcMap does not seem to make this important distinction in any of its documentation. 
#`Monte Carlo method`
MC<- moran.mc(eth_eth9$confirmed_all, lw, nsim=999, alternative="greater")

# View results (including p-value)
MC

# Plot the Null distribution (note that this is a density plot instead of a histogram)
plot(MC)

# The p-value is >0.05 suggesting that there is no global spatial correlations of malaria incidence. 
# ## Expected value without considering the age strata
#  <- d3 %>% 
#   group_by(Year.ID) %>% 
#   mutate(r = sum(O.malaria)/sum(Pop), # overall rate stratified by year
#          E2 = r*Pop,
#          sir = Y/Ex)

##Covariate data##
#'load covariate data: static:AI, TWI, PET, Elevation, Slope, access to cities, access to health facilities, distance to weight,ratio of urban population
#' dynamic(annual): IRS, ITN, VIIRS
#' dynamic(monthly): EVI, day light temperature, light time temperature, TCB, TCW, Rainfall
covariate <- readRDS("Z:/Yalem/Ethiopia_malaria/Ahmed_data/covariate_scaled.rds") 

#saveRDS(covariate, file = "Z:/Yalem/Ethiopia_malaria/Ahmed_data/covariate_scaled.rds")  # unscaled

####Modeling############
## Inference using INLA
## INCIDENCE INTENSITY RATE 
#' In any Bayesian analysis, the choice of the prior may have a considerable impact on the results.
#' By default, minimally informative priors are specified on:
#' (i) the log of the unstructured effect precision log(????) ∼ logGamma(1, 0.0005), and
#' (ii) the log of the structured effect precision log(????u) ∼ logGamma(1, 0.0005).
#' However, for intrinsic models such as iCAR, the precision matrix depending on the random parameter and the neighboring structure)
#' In this model no prior information included in the model.

## I read two options of model fitting 
#' option 1: a non-parametric model for the time trend 
#' option 2: parametric model

## OPTION 1#
#' Data preparation
##Merge covariates with case data

d4 <- left_join(d4, covariate)
d4$year <- as.numeric(d4$year) #year as numeric
d4$district <- as.factor(d4$district) #district as factor

##create the index vectors for the counties and years that will be used to specify the random effects of the model
d4$idtime <- 1 + d4$year - min(d4$year)
d4$idarea <- as.numeric(d4$district) # is the vector with the indices of districts 1 to 967
d4$idarea1 <- d4$idarea
d4$idtime2 <- d4$idtime
d4$idtime1 <- d4$idtime# is the vector with the indices of years and has 1 to 3 number of years
d4$ID.area.year = seq(1,length(d4$district)) #interaction

##Assign weights to the neighbors
#create adjacency matrix using `spdep` package
nb2 <- poly2nb(map.eth, queen = TRUE, row.names = map.eth$district)
head(nb2)

##Convert the adjacency matrix into a file in the INLA format and save

nb2INLA("ethiopi.map.adj", nb2)

#Visualize the graph and get summary - if necessary 
g<-inla.read.graph(filename ="ethiopi.map.adj")
summary(g)


# equal neighboring polygon will be assigned equal weight when computing the neighboring mean income values.
lw <- nb2listw(nb2, style="W", zero.policy=TRUE)
##To see the weight of the first polygon’s neighbors type:
lw$weights[1]

##read the graph
eth.map2.graph <- inla.read.graph(filename = "ethiopi.map.adj")
##Adjacency matrix: rows and columns identify areas; dotes identify neighbors
image(inla.graph2matrix(g), xlab = "", ylab = "")

##' Exploratory analysis : Choropleth map
#Plot raw incidence rate

# let's find a natural interval with quantile breaks
ni = classIntervals(d4$ir,
                    n = 6,
                    style = 'quantile')$brks

# this function uses above intervals to create categories
labels <- c()
for(i in 1:length(ni)){
  labels <- c(labels, paste0(round(ni[i], 0),
                             "-",
                             round(ni[i + 1], 0)))
}

labels <- labels[1:length(labels)-1]

d4$cat <- cut(d4$ir,
 
                         breaks = ni,
             labels = labels,
             include.lowest = T)
levels(d4$cat) # let's check how many levels it has (6)

# label NAs, too
lvl <- levels(d4$cat)
lvl[length(lvl) + 1] <- "No data"
d4$cat <- factor(d4$cat, levels = lvl)
d4$cat[is.na(d4$cat)] <- "No data"
levels(d4$cat)


#test positivity
#check outlier
d4 |> 
  filter(confirmed_all>test_performed) |> 
  View()

np = classIntervals(d4$tp,
                    n = 6,
                    style = 'quantile')$brks

# this function uses above intervals to create categories
labels <- c()
for(i in 1:length(np)){
  labels <- c(labels, paste0(round(np[i], 0),
                             "-",
                             round(np[i + 1], 0)))
}

labels <- labels[1:length(labels)-1]

d4$tp_cat <- cut(d4$tp,
              
              breaks = np,
              labels = labels,
              include.lowest = T)
levels(d4$tp_cat) # let's check how many levels it has (6)

# label NAs, too
lvl <- levels(d4$tp_cat)
lvl[length(lvl) + 1] <- "No data"
d4$tp_cat <- factor(d4$tp_cat, levels = lvl)
d4$tp_cat[is.na(d4$tp_cat)] <- "No data"
levels(d4$tp_cat)

#Merge the shapefile with the dataset
map_incidence <- left_join(sf_eth, d4)

#Incidence
p_incidence <- ggplot() +
    geom_sf(data= map_incidence, aes(fill = cat),
          color = NA) +
  geom_sf(data= admin1_shp, colour = 'black', fill = NA, size = .15 ) + 
  scale_fill_manual(name= "Incidence rate per 1000 population",
                    values=c('#ffd080', '#f4a77a', '#e18079','#c35e7d', '#9b4283','#6c2d83',"grey80"),
                    labels=c("0-0","0-1","1-4", "4-10","10-29","29-708","No data"),
                    drop=F) +
  facet_wrap(~year,
             dir = "h",
             ncol = 3) + 
  labs(fill = "Malaria IIR/1000 population") +
  theme_bw() + 
  theme(line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  ) +
  theme(legend.position = "bottom", 
        legend.key.width = unit(1.5, "cm"), 
        legend.key.height = unit(0.2, "cm"), 
        strip.text.y = element_text(angle = 180))

  #ggsave(paste0(my.dir, '/raw_incidence.png'))
#Test positivity
p_pt <- ggplot() +
  geom_sf(data= map_incidence, aes(fill = tp_cat),
          color = NA) +
  geom_sf(data= admin1_shp, colour = 'black', fill = NA, size = .15 ) + 
  scale_fill_manual(name= "Malaria test positivity rate (%)",
                    values=c('#ffd080', '#f4a77a', '#e18079','#c35e7d', '#9b4283','#6c2d83',"grey80"),
                    labels=c("0-6","6-12","12-20", "20-31","31-52","52-100","No data"),
                    drop=F) +
  facet_wrap(~year,
             dir = "h",
             ncol = 3) + 
  labs(fill = "Malaria test positivity rate (%)") +
  theme_bw() + 
  theme(line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  ) +
  theme(legend.position = "bottom", 
        legend.key.width = unit(1.5, "cm"), 
        legend.key.height = unit(0.2, "cm"), 
        strip.text.y = element_text(angle = 180))

library(ggpubr)
ggarrange(p_incidence, p_pt,
          # labels = c('Incidence rate', 'Test positivity rate'),
          ncol = 1, nrow = 2)
#ggsave(paste0(my.dir, '/incidence_tp.png'))
## summarize reported cases by year for each region
d4 |> 
  group_by(district, year) |> 
  summarize(all_age = sum(confirmed_all, na.rm = T),
            u5 = sum(confirmed_u5, na.rm = T),
            reported5_14 = sum(confirmed_5_14, na.rm = T),
            ov15 = sum(confirmed_15, na.rm = T)) |> 
          knitr::kable()
            

## summarize reported cases by year 
d4 |> 
  group_by(year) |> 
  summarize(all_age_test = sum(test_performed, na.rm = T),
            all_cases = sum(confirmed_all, na.rm = T),
            u5 = sum(confirmed_u5, na.rm = T),
            reported5_14 = sum(confirmed_5_14, na.rm = T),
            ov15 = sum(confirmed_15, na.rm = T)) |> 
  knitr::kable()

#average confirmed cases reported across districts throughout the study periods
d4 |> 
  group_by(district) |> 
  summarize(all_age_test = sum(test_performed, na.rm = T),
            all_cases = sum(confirmed_all, na.rm = T),
            u5 = sum(confirmed_u5, na.rm = T),
            reported5_14 = sum(confirmed_5_14, na.rm = T),
            ov15 = sum(confirmed_15, na.rm = T)) |> 
  knitr::kable()

#check coliniarity
library("olsrr")


d_cov=d4 %>% 
  ungroup() |> 
  select(confirmed_all,AI,TWI, PET, Elevation , Slope , access_cities , AccessHF ,DTW, Urban , AM , IRS , ITN ,VIIRS , EVI , Day , Night , TCB , TCW , Rainfall)

#fit linear regression model
#the goal of this model is to predict the incidence 
my_model <- lm(confirmed_all~., data = d_cov)

## 1. CORRELATION MATRIX
library("corrplot")
data_cor <- d_cov |> 
  select(-confirmed_all) 

corrplot(cor(data_cor), method = "number")

##2. Test for Multicollinearity with Variance Inflation Factors (VIF)

ols_vif_tol(my_model) |> 
  knitr::kable()

#we need to remove variables with VIF>10 one by one
#AI, PET,EVI,Nigtht(ITN), TWI has a correlation/VIF higher values

data_cora <- d_cov |> 
select(-AI,-PET,-EVI,-TCW,-VIIRS, -access_cities, -Night, -Elevation, -Slope)

my_model2 <- lm(confirmed_all~., data = data_cora)


data_cor2 <- data_cora |> 
  select(-confirmed_all) 

corrplot(cor(data_cor2), method = "number")

##2. Test for Multicollinearity with Variance Inflation Factors (VIF)

ols_vif_tol(my_model2) |> 
  knitr::kable()
##FIT THE MODEL
#'model1 - covariate only
#'The model includes a Besag-York-Mollié (BYM) model

formula1<-  confirmed_all ~ + TWI +  AccessHF + DTW + Urban +
  AM + IRS + ITN + Day + TCB +  Rainfall 


model1<- inla(formula1,
              family = "poisson", 
              data = d4, E = population,
              control.predictor = list(link=1),
              control.compute=list(dic=TRUE,cpo=TRUE,  waic = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose = TRUE)

#output
summary(model1)
round(model1$summary.fixed[,1:5],3)
#model diagnostic
plot(d4$confirmed_all,model1$summary.fitted.values$mean)

plot(log(d4$confirmed_all),log(model1$summary.fitted.values$mean),
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")##not dispersed

plot(model1, prior = T)

#'Model2: include spatially structured component with index variable `idarea`  f(idarea,model='bym',graph=eth.map.graph)
#'spatial random effect defined `ui` as follows a CAR distribution and `vi`is an independent and identically distributed normal variable

 formula2<- confirmed_all ~ + TWI +  AccessHF + DTW + Urban +
   AM + IRS + ITN + Day + TCB +  Rainfall +
   f(idarea, model='bym',graph=g)



model2<- inla(formula2,
              family = "poisson", 
              data = d4, E = population,
              control.predictor = list(link=1),
              control.compute=list(dic=TRUE,config = TRUE, cpo=TRUE,  waic = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose = TRUE)

#Out put
summary(model2)

#summary fixed effect
round(model2$summary.fixed[,1:5],3)

#The summary statistics of the predicted values can be accessed by
summary(model2$summary.fitted.values$mean)

#Inspect the model output
plot(d4$confirmed_all,model2$summary.fitted.values$mean)
plot(log(d4$confirmed_all),log(model2$summary.fitted.values$mean),
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")##not dispersed 


#marginals.hyperpar
plot(model2$marginals.hyperpar$`Precision for idarea (spatial component)`) #someting wrong here


#Model3- adding the time component 'rw2' 
#'temporal random effect a random walk in time of first order (RW) or a random walk in time of second order (RW2)
formula3<-  confirmed_all ~ + TWI +  AccessHF + DTW + Urban +
  AM + IRS + ITN + Day + TCB +  Rainfall +
  f(idarea, model='bym',graph=g) + 
  f(idtime, model='rw2') 



model3<- inla(formula3,
              family = "poisson", 
              data = d4, E = population,
              control.predictor = list(link=1),
              control.compute=list(dic=TRUE,config = TRUE, cpo=TRUE,  waic = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose = TRUE)

#output
summary(model1)
summary(model2)
summary(model3)
round(model3$summary.fixed[,1:5],3)

##The summary statistics of the predicted values can be accessed by
summary(model3$summary.fitted.values$mean)
.
##Inspect the model output
plot(d4$confirmed_all,model3$summary.fitted.values$mean)
plot(log(d4$confirmed_all),log(model3$summary.fitted.values$mean),
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")##a bit dispersed 

#marginals.hyperparameter
plot(model3$marginals.hyperpar$`Precision for idarea (spatial component)`) #still?
plot(model3$marginals.hyperpar$`Precision for idarea (iid component)`) 
plot(model3$marginals.hyperpar$`Precision for idtime`) 

##Model 4: include a temporal random effect,`iid`
#'unstructured temporal effect that is modeled with an independent and identically distributed normal variable: `f(idtime1,model = 'iid')`
 
formula4<-  confirmed_all ~ + TWI +  AccessHF + DTW + Urban + AM + IRS + ITN + Day + TCB +  Rainfall +
  f(idarea, model='bym',graph=g) + 
  f(idtime, model='rw2') + 
  f(idtime1,model = 'iid') 



model4<- inla(formula4,
              family = "poisson", 
              data = d4, E = population,
              control.predictor = list(link=1),
              control.compute=list(dic=TRUE,config = TRUE, cpo=TRUE,  waic = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose = TRUE)

#Output
summary(model3)
summary(model4)
round(model4$summary.fixed[,1:5],3)

summary(model4$summary.fitted.values$mean)

plot(d4$confirmed_all,model4$summary.fitted.values$mean)
plot(log(d4$confirmed_all),log(model4$summary.fitted.values$mean),
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")##a bit dispersed 


##Model5: adding the interaction#
#'an interaction between space and time that can be specified the interaction between `vi (i.i.d)` and an unstructured temporal effect `ϕj(i.i.d)`
#'assumes no spatial  or temporal structure on `f(ID.area.year,model = 'iid')`
formula5<- confirmed_all ~ + Night +  AccessHF + DTW + Urban + IRS +
  Elevation + Day + TCB + Rainfall + Slope +
  f(idarea, model='bym',graph=g) + 
  f(idtime, model='rw2') + 
  # f(idtime1,model = 'iid') +
  f(ID.area.year, model = 'iid')




model5<- inla(formula5,
              family = "poisson", 
              data = d4, E = population,
              control.predictor = list(link=1),
              control.compute=list(dic=TRUE,config = TRUE, cpo=TRUE,  waic = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose = TRUE)

summary(model5)

# ##opt2
# formula5_2<-  cconfirmed_all ~ + Night +  AccessHF + DTW + Urban +
#   AM + IRS + Elevation + Day + TCB +  Rainfall + Slope +
#   f(idarea, model='bym',graph=g) + 
#   f(idtime, model='rw2') + 
#   # f(idtime1,model = 'iid') +
#   f(idtime2, idarea, model = 'iid')
# 
# model5_2<- inla(formula5,
#               family = "poisson", 
#               data = d4, E = population,
#               control.predictor = list(link=1),
#               control.compute=list(dic=TRUE,config = TRUE, cpo=TRUE,  waic = TRUE),
#               control.inla = list(int.strategy = "eb"),
#               verbose = TRUE)
# 
# summary(model5_2)
# summary(model5)

output <- round(model5$summary.fixed[,1:5],3) |> 
 data.frame() 
  write.csv(output, file = "model5.csv")

##The summary statistics of the predicted values can be accessed by
summary(model5$summary.fitted.values$mean)

plot(d4$confirmed_all,model5$summary.fitted.values$mean)

hist(model5$summary.fitted.values$mean)

plot(log(d4$confirmed_all),log(model5$summary.fitted.values$mean),
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")##a bit dispersed 



temporal.malaria.CAR <- lapply(model5$marginals.random$idarea, 
                               function(X){
                                 marg <- inla.tmarginal(function(x) exp(x), X)
                                 inla.emarginal(mean, marg)
                               })


temporal.malaria <- lapply(model5$marginals.random$idtime, 
                           function(X){
                             marg <- inla.tmarginal(function(x) exp(x), X)
                             inla.emarginal(mean, marg)
                           })



plot(seq(1,3),seq(0.5,5,length=3),type="n",xlab="Year",ylab=expression(exp(phi[Year])))

lines(unlist(temporal.malaria), col = "red")
lines(unlist(temporal.malaria.CAR),lty=2)
abline(h=1,lty=1)

#marginals.hyperparameter 
plot(model5$marginals.hyperpar$`Precision for idarea (spatial component)`)
plot(model5$marginals.hyperpar$`Precision for idarea (iid component)`) 
plot(model5$marginals.hyperpar$`Precision for idtime`) 
#plot(model5$marginals.hyperpar$`Precision for idtime1`) 
plot(model5$marginals.hyperpar$`Precision for ID.area.year`)


####Residual
ob = d4$confirmed_all
fitted = 10**(10**model5$summary.fitted.values$mean)
resid = ob - fitted
hist(resid)
plot(ob,fitted)
y1_index <- d4$year == d4$year[1]
resid_y1 <- resid[y1_index]
sf_eth$resid_y1 <- resid_y1

plot(sf_eth["resid_y1"])

head(d4)

plot(d4$confirmed_all, fitted)
pit <- model5$cpo$pit
pit <- pmax(pit, 1e-3)
pit <- pmin(pit, 1 - 1e-3)

pit_z <- qnorm(pit)
hist(pit_z,breaks=10,main="",xlab="PIT")

sf_eth$pit_z_y1 <- pit_z[y1_index]

plot(sf_eth["pit_z_y1"])


#Model6: an interaction assumes interaction between unstructured spatial component `vi(CAR)` and `yt(i.i.d)`
#'an interaction between space and time that can be specified the interaction between `vi (i.i.d)` and structured temporal effect `????t(rw2)`
#'assumes no spatial  or temporal structure on `f(ID.area.year,model = 'iid')`

formula6<-  confirmed_all ~ + TWI +  AccessHF + DTW + Urban + AM + IRS + ITN + Day + TCB +  Rainfall +
  f(idarea, model='bym',graph=g) + 
  f(idtime, model='rw2') + 
  #f(idtime1,model = 'iid') +
  f(idarea1,
    model = 'iid',
    group = idtime2,
    control.group = list(model = "rw2")
  )


model6<- inla(formula6,
              family = "poisson", 
              data = d4, E = population,
              control.predictor = list(link=1),
              control.compute=list(dic=TRUE,
                                   config = TRUE,
                                   cpo=TRUE,
                                   waic = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose = TRUE)
summary(model5)
summary(model6)
##Inspect the model output
round(model6$summary.fixed[,1:5],3)

##The summary statistics of the predicted values can be accessed by
summary(model6$summary.fitted.values$mean)

plot(d4$confirmed_all,model6$summary.fitted.values$mean)

plot(log(d4$confirmed_all),log(model6$summary.fitted.values$mean) ,
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")##a bit dispersed 


##Model 7: assumes interaction between  unstructured temporal effect and the spatially structured main effect
formula7<-  confirmed_all ~ + TWI +  AccessHF + DTW + Urban +  AM + IRS + ITN + Day + TCB +  Rainfall +
  f(idarea, model='bym',graph=g) + 
  f(idtime, model='rw2') + 
#f(idtime1,model = 'iid') +
  f(idtime2,
    model = 'iid',
    group = idarea1,
    control.group = list(model = "besag",
                         graph=g)
  )


model7<- inla(formula7,
              family = "poisson", 
              data = d4, E = population,
              control.predictor = list(link=1),
              control.compute=list(dic=TRUE,config = TRUE, 
                                   cpo=TRUE,  waic = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose = TRUE)

summary(model6)
summary(model7)

round(model7$summary.fixed[,1:5],3)

hist(model7$summary.fitted.values$mean)
summary(d4$confirmed_all)

##Model 8: spatially and temporally structured
formula8<-  confirmed_all ~ + TWI +  AccessHF + DTW + Urban +  AM + IRS + ITN + Day + TCB +  Rainfall +
  f(idarea, model='bym',graph=g) + 
  f(idtime, model='rw2') + 
  #f(idtime1, model = 'iid') +
  f(idarea1,
    model = 'besag',
    graph=g,
    group = idtime2,
    control.group = list(model = 'rw2'))


model8<- inla(formula8,
              family = "poisson", 
              data = d4, E = population,
              control.predictor = list(link = 1),
              control.compute=list(dic=TRUE,config = TRUE,
                                   cpo=TRUE,  waic = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose = TRUE)
#Output
summary(model5)
summary(model8)

round(model8$summary.fixed[,1:5],3)

hist(model8$summary.fitted.values$mean)
summary(d4$confirmed_all)

fixed_effect <-  
  model8$summary.fixed |> 
  data.frame() |> 
  view()

plot(d4$confirmed_all,model8$summary.fitted.values$mean)
plot(log(d4$confirmed_all),log(model8$summary.fitted.values$mean),
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")##a bit dispersed 


####################################################### OPTION2 #####################################################
#### The model assumes a linear time trend in each area.
##Spatio-temporal model with parametric time trends that expresses the logarithm of the relative risks (OR IRR)
#'model1 - area random effect `ui +vi` 
f_besag<-  confirmed_all ~ + AI + TWI + PET + Elevation + Slope + access_cities + AccessHF + DTW + Urban +
  AM + IRS + ITN + VIIRS + EVI + Day + Night + TCB + TCW + Rainfall +
  f(idarea,  model = 'bym',  graph = g)  
#f(idarea1, idtime, model = 'iid', constr=TRUE) + idtime

m_besag<- inla(f_besag,
               family = "poisson", 
               data = d4, E = population,
               control.predictor = list(link=1),
               control.compute=list(dic=TRUE,config = TRUE,
                                    cpo=TRUE,  waic = TRUE),
               control.inla = list(int.strategy = "eb"),
               verbose = TRUE)

summary(m_besag)
##Inspect the model output
round(m_besag$summary.fixed[,1:5],3)

##The summary statistics of the predicted values can be accessed by
summary(m_besag$summary.fitted.values$mean)

plot(d4$confirmed_all,m_besag$summary.fitted.values$mean)

plot(log(d4$confirmed_all),log(m_besag$summary.fitted.values$mean),
     xlab="Observed Values",ylab="Mean Post. Pred. Distr.")##a bit dispersed 


#marginals.hyperparameter 
plot(m_besag$marginals.hyperpar$`Precision for idarea (spatial component)`)
plot(m_besag$marginals.hyperpar$`Precision for idarea (iid component)`)

####Residual
ob = d4$confirmed_all
fitted = m_besag$summary.fitted.values$mean
resid = ob - fitted
hist(resid)

y1_index <- d4$idtime ==d4$idtime[1]
resid_y1 <- resid[y1_index]
sf_eth$resid_y1 <- resid_y1

plot(sf_eth["resid_y1"])

head(d4)

plot(d4$confirmed_all, fitted)
pit <- m_besag$cpo$pit
pit <- pmax(pit, 1e-3)
pit <- pmin(pit, 1 - 1e-3)

pit_z <- qnorm(pit)
hist(pit_z,breaks=10,main="",xlab="PIT")

sf_eth$pit_z_y1 <- pit_z[y1_index]

plot(sf_eth["pit_z_y1"])




#'include the time effect - `global linear trend effect`- **idtime** denotes the global trend
#' **f(idarea1, idtime, model = "iid")** is the differential time trend

f_besag2<- confirmed_all ~ + AI + TWI + PET + Elevation + Slope + access_cities + AccessHF + DTW + Urban +
  AM + IRS + ITN + VIIRS + EVI + Day + Night + TCB + TCW + Rainfall +
  f(idarea,  model = 'bym',  graph = g)+  
 f(idarea1, idtime, model = 'iid', constr=TRUE) + idtime


######   Yalem to measure the posterior exceeeding probability you need to return return.marginals.predictor=TRUE
m_besag1<- inla(f_besag2,
                family = "poisson", 
                data = d4, E = population,
                control.predictor = list(link=1),
                control.compute=list(dic=TRUE,config = TRUE,
                                     cpo=TRUE,  waic = TRUE,return.marginals.predictor=TRUE),
                control.inla = list(int.strategy = "eb"),
                verbose = TRUE)

##Output
summary(m_besag1)
summary(m_besag1$summary.fitted.values$mean) #posterioir mean (IIR)
plot(m_besag1$marginals.hyperpar$`Precision for idarea (iid component)`)
plot(m_besag1$marginals.hyperpar$`Precision for idarea (spatial component)`)

plot(log(d4$confirmed_all), log(m_besag1$summary.fitted.values$mean))



  p <- autoplot(m_besag1$summary.fixed$mean)
  cowplot::plot_grid(plotlist = p)
#Global linear temporal trend for malaria incidence ratio in Ethiopia.
x <- seq(1,6) # Years
plot(x,m_besag1$summary.fixed[10,1]*x, type="l", main="",xlab="t",ylab=expression(beta*t))
lines(m_besag1$summary.fixed[10,3]*x,lty=2) # lower bound - 0.025quant
lines(m_besag1$summary.fixed[10,5]*x,lty=2) # upper bound - 0.975quant
lines((m_besag1$summary.random$idarea1[10,2] + m_besag1$summary.fixed[10, 1])*x, col = "red")
lines((m_besag1$summary.random$idarea1[1,2] + m_besag1$summary.fixed[10, 1])*x, col = "red")
lines((m_besag1$summary.random$idarea1[101,2] + m_besag1$summary.fixed[10, 1])*x, col = "red")

##Residual
ob = d4$confirmed_all
fitted = m_besag1$summary.fitted.values$mean
resid = ob - fitted
hist(resid)

y1_index <- d4$idtime == d4$idtime[1]
resid_y1 <- resid[y1_index]
sf_eth$resid_y1 <- resid_y1

plot(sf_eth["resid_y1"])

head(d4)

plot(d4$confirmed_all, fitted)
pit <- m_besag1$cpo$pit
pit <- pmax(pit, 1e-3)
pit <- pmin(pit, 1 - 1e-3)

pit_z <- qnorm(pit)
hist(pit_z,breaks=10,main="",xlab="PIT")

sf_eth$pit_z_y1 <- pit_z[y1_index]

plot(sf_eth["pit_z_y1"])

##the best fit model is model 5
####Check the model
model_data = d4
selected_model = model5 
summary(selected_model)

#add the model in the data
preds <- as.data.frame(model_data) %>%
 dplyr:: mutate(pred_mean = selected_model$summary.fitted.values$mean,
         pred_lower = selected_model$summary.fitted.values$`0.025quant`,
         pred_upper = selected_model$summary.fitted.values$`0.975quant`) %>%
 dplyr:: select(district, year, zone, pred_mean, pred_lower, pred_upper)

##########
#add the predicted value in the data

# model_data$pred_mean <- round(exp(log(m_besag1$summary.fitted.values[, "mean"])), digits = 3)
# model_data$pred_lower <- round(exp(log(m_besag1$summary.fitted.values[, "0.025quant"])), digits = 3)
# model_data$Pred_upper <- round(exp(log(m_besag1$summary.fitted.values[, "0.975quant"])), digits = 3)

library(hablar)
## To make maps of the relative risks and Incidence rate, we first merge map_sf with d so map_sf has columns RR, LL and UL.
# predic_plot <- model_data |> 
#   select(region, zone, district, year, pred_mean, pred_lower, Pred_upper) |> 
#   convert(chr(year))
# 
# predic_plot <- predic_plot |> 
#   pivot_wider(names_from = "year",
#               names_sep = "_",
#               values_from = c("pred_mean", "pred_lower", "Pred_upper"))

mapredict <- left_join(sf_eth, preds)

# mapredict@data[1:2, ]
# eth_sfpredict <- st_as_sf(mapredict)
# 
# library(tidyr)
# map_mean <- gather(eth_sfpredict, year, pred_mean, paste0("pred_mean_", 2019:2021)) #mean
# map_lower <- gather(eth_sfpredict, year, pred_lower, paste0("pred_lower_", 2019:2021)) # lower bound
# map_upper <- gather(eth_sfpredict, year, Pred_upper, paste0("Pred_upper_", 2019:2021))#upper bound

# #Column year of map_sf contains the values SIR.2019, . . . , SIR.2021.
# eth_sfpredict$year <- as.integer(substring(eth_sfpredict$year, 11, 14))
# eth_sfpredict$year <- as.integer(substring(eth_sfpredict$year, 5, 8))
# eth_sfpredict$year <- as.integer(substring(eth_sfpredict$year, 6, 9))


## Incidence rate - something is wrong in the IR plot  ---- You need to load viridis for the code to work :-)

p_mean <- ggplot() +
  geom_sf(data= mapredict, aes(fill = pred_mean*1000), colour = NA) +
  geom_sf(data = admin1_shp, colour = 'black', fill = NA, size = .15) +
    scale_fill_viridis(option = "magma", direction = -1) +
  facet_wrap(~year,
             dir = "h",
             ncol = 3) +
  # guides(fill = guide_colorbar(title = "IIR/1000 population", 
  #                              title.position = "top")) +
#labs(fill = "Malaria IIR/1000 population") +
  theme_bw() + 
  theme(line = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  theme(legend.position = "none")
# ggsave(paste0(my.dir, '/IIR_preds_mean.png'))

#Lower
p_lower <- ggplot() +
  geom_sf(data= mapredict, aes(fill = pred_lower*1000), 
          color = NA) +
  geom_sf(data = admin1_shp, colour = 'black', fill = NA, size = .15) +
  scale_fill_viridis(option = "magma", direction = -1) +
  facet_wrap(~year,
             dir = "h",
             ncol = 3) +
  # guides(fill = guide_colorbar(title = "Lower bound CrI/1000 population", 
  #                              title.position = "top")) +
  # labs(fill = "Malaria IIR/1000 population") +
  theme_bw() + 
  theme(line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  ) +
  theme(legend.position = "none")
# ggsave(paste0(my.dir, '/IIR_preds_lower.png'))

#upper
p_upper <- ggplot() +
  geom_sf(data= mapredict, aes(fill = pred_upper*1000), 
          color = NA) +
  geom_sf(data = admin1_shp, colour = 'black', fill = NA, size = .15) +
  scale_fill_viridis(option = "magma", direction = -1) +
  facet_wrap(~year,
             dir = "h",
             ncol = 3) +
  # guides(fill = guide_colorbar(title = "Upper bound CrI/1000 population", 
  #                              title.position = "top")) +
  labs(fill = "Malaria IIR/1000 population") +
  theme_bw() + 
  theme(line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  ) +
  theme(legend.position = "bottom", 
        legend.key.width = unit(1.5, "cm"), 
        legend.key.height = unit(0.2, "cm"), 
        strip.text.y = element_text(angle = 180))



ggarrange(p_lower, p_mean, p_upper,
          labels = c('Lower', 'Mean', 'Upper'),
          ncol = 1, nrow = 3)

#ggsave(paste0(my.dir, '/uncertaiinity.png'))
################### Exceeding posterior probability
#Here we calculate our exceedence probabilities, in this case, I calculate the probability that a county has 
# a value greater than 0.17.
d4 <- as.data.frame(d4)
a2<-0.07
d4$inlaprob<-unlist(lapply(model5$marginals.fitted.values, function(X){
  1-inla.pmarginal(a2, X)
}))

hist(d4$inlaprob)

d4$exceed.cut<-cut(d4$inlaprob,breaks =  c(0,0.1,0.2,0.5,0.8,.95,1), include.lowest = T)


## Now I think the plot should not be a problem because you have the probability and the counties saved in d4
# dataframe. I have bot made the plot because It will be up to you to chose the cut off limit for the upper 
#or lower probability 
map_probab <- left_join(sf_eth, d4)
p <- ggplot() +
  geom_sf(data= map_probab, aes(fill = exceed.cut)) +
  scale_fill_manual(name= "Incidence rate per 1000 population",
                    values=c('#ffd080', '#f4a77a', '#e18079', '#c35e7d', '#9b4283', '#6c2d83'),
                    labels=c("[0,0.1]",     "(0.1,0.2]",     "(0.2,0.5]",     "(0.5,0.8]",   "(0.8,0.95]", "(0.95,1]"),
                    drop=F) +
  facet_wrap(~year,
             dir = "h",
             ncol = 3) +
  guides(fill=guide_legend(
    direction = "horizontal",
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )
  )

p1 <- p + theme_bw() +
  theme(panel.background = element_blank(),
        legend.background = element_blank(),
        legend.position = c(.45, .04),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "white", size = 0),
        plot.title = element_text(size=20, color="#6c2d83", hjust=0.5, vjust=-10),
        plot.subtitle = element_text(size=14, color="#bd5288", hjust=0.5, vjust=-15, face="bold"),
        plot.caption = element_text(size=9, color="grey60", hjust=0.5, vjust=9),
        axis.title.x = element_text(size=7, color="grey60", hjust=0.5, vjust=5),
        legend.text = element_text(size=10, color="grey20"),
        legend.title = element_text(size=11, color="grey20"),
        strip.text = element_text(size=12),
        plot.margin = unit(c(t=-2, r=-2, b=-2, l=-2),"lines"), #added these narrower margins to enlarge map
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) 
ggsave(paste0(my.dir, '/marginal_prob.png'))

#####################   The most positive and negative covariates for that since you are measuring the values of incidence per county it will be the 
# covariates at each county times it respected beta value. 

d_cov=d4 %>% select( Night, AccessHF, DTW , Urban , IRS ,  Elevation, Day , TCB , Rainfall , Slope )

#d_cov$AI=m_besag1$summary.fixed$mean[2]*d_cov$AI
#d_cov$TWI=m_besag1$summary.fixed$mean[3]*d_cov$TWI
#d_cov$PET=m_besag1$summary.fixed$mean[4]*d_cov$PET
d_cov$Elevation=model5$summary.fixed$mean[5]*d_cov$Elevation
d_cov$Slope=model5$summary.fixed$mean[6]*d_cov$Slope
#d_cov$access_cities=m_besag1$summary.fixed$mean[7]*d_cov$access_cities
d_cov$AccessHF=model5$summary.fixed$mean[8]*d_cov$AccessHF
d_cov$DTW=model5$summary.fixed$mean[9]*d_cov$DTW
d_cov$Urban=model5$summary.fixed$mean[10]*d_cov$Urban
#d_cov$AM=m_besag1$summary.fixed$mean[11]*d_cov$AM
d_cov$IRS=model5$summary.fixed$mean[12]*d_cov$IRS
#d_cov$ITN=m_besag1$summary.fixed$mean[13]*d_cov$ITN
#d_cov$VIIRS=m_besag1$summary.fixed$mean[14]*d_cov$VIIRS
#d_cov$EVI=m_besag1$summary.fixed$mean[15]*d_cov$EVI
d_cov$Day=model5$summary.fixed$mean[16]*d_cov$Day
d_cov$Night=model5$summary.fixed$mean[17]*d_cov$Night
d_cov$TCB=model5$summary.fixed$mean[18]*d_cov$TCB
#d_cov$TCW=m_besag1$summary.fixed$mean[19]*d_cov$TCW
d_cov$Rainfall=model5$summary.fixed$mean[20]*d_cov$Rainfall


maxColumnNames <- apply(d_cov,1,function(row) colnames(d_cov)[which.max(row)])     ### most positive covariate
minColumnNames <- apply(d_cov,1,function(row) colnames(d_cov)[which.min(row)])    ### most negative covariate

resDf <- cbind(data.frame(d_cov),data.frame(maxColumnNames = maxColumnNames))
resDf <- cbind(data.frame(resDf),data.frame(minolumnNames = minColumnNames))

cov_plot <- left_join(sf_eth, resDf)

#plot
ggplot() +
  geom_sf(data= cov_plot, aes(fill = maxColumnNames), 
          color = "transparent") +
  scale_fill_viridis(option = "magma") +
  guides(fill = guide_colorbar(title = "minimu", 
                               title.position = "top")) +
  labs(fill = "least important") +
  theme_bw() + 
  theme(line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  ) +
  theme(legend.position = "bottom", 
        legend.key.width = unit(1.5, "cm"), 
        legend.key.height = unit(0.2, "cm"), 
        strip.text.y = element_text(angle = 180))

ggsave(paste0(my.dir, '/cov_min.png'))
#############################
#############################
### Model for monthly data





dm <- readRDS("Z:/Yalem/Ethiopia_malaria/Ahmed_data/test_case_monthlydata_update.rds") |> 
  rename("year" = "Year")
#observed cases


dm <- dm %>% 
  mutate(r = sum(!is.na(d4$confirmed_all))/sum(population), # overall rate stratified by year
         Ex = r*population,
         sir = confirmed_all[!is.na(confirmed_all)]/Ex)
##explore the data before analysis##
#subset the data by years and inspect the data: the dataset doesn't seems to be well behaved
d4_2019 <- d4 |> filter(year==2019)
d4_2020 <- d4 |> filter(year==2020)
d4_2021 <- d4 |> filter(year==2021)

##map : SIR 
#date is a character format in the `sf_eth` dataframe
library(lubridate)
sf_eth  <- sf_eth |> 
  mutate(date = as.Date(date))  
eth_eth9 <- right_join(sf_eth, d4_2019,by= c("district"))



eth_eth9 <-eth_eth9[1:11603,]



#left_join(d4_2019, covariate,by= c("month","year","region","zone","district"),all.x=TRUE, all.y=TRUE)
##map : SIR 

tm_shape(eth_eth9[eth_eth9$Month=="June",]) +
  tm_fill(col="sir", style="quantile", n=8, palette="Reds") +
  tm_legend(outside=TRUE)


#Computing the Moran’s I statistic
#Recall that the Moran’s I value is the slope of the line that best fits the relationship between neighboring malaria cases
#and each polygon’s number of malaria case in the dataset.
I <- moran(eth_eth9$confirmed_all, lw, length(nb2), Szero(lw))[1]
I
#Analytical method
#Performing a hypothesis test
#There are two methods to testing this hypothesis: an analytical method and a Monte Carlo method.
#To run the Moran’s I analysis using the analytical method, use the `moran.test` function.
moran.test(eth_eth9$confirmed_all,lw, alternative="greater")

#' Note that ArcMap adopts this analytical approach to its hypothesis test however,
#'  it implements a two-sided test as opposed to the one-sided test adopted in the above example (i.e. alternative = "greater"). 
#' A two-sided p-value is nothing more than twice the one-sided p-value.
#'  Unfortunately, ArcMap does not seem to make this important distinction in any of its documentation. 
#`Monte Carlo method`
MC<- moran.mc(eth_eth9$confirmed_all, lw, nsim=999, alternative="greater")

# View results (including p-value)
MC

# Plot the Null distribution (note that this is a density plot instead of a histogram)
plot(MC)

# The p-value is >0.05 suggesting that there is no global spatial correlations of malaria incidence. 
# ## Expected value without considering the age strata
#  <- d3 %>% 
#   group_by(Year.ID) %>% 
#   mutate(r = sum(O.malaria)/sum(Pop), # overall rate stratified by year
#          E2 = r*Pop,
#          sir = Y/Ex)

##'load covariate data: static:AI, TWI, PET, Elevation, Slope, access to cities, access to health facilities, distance to weight,ratio of urban population
##' dynamic(annual): IRS, ITN, VIIRS
##' dynamic(monthly): EVI, day light temperature, light time temperature, TCB, TCW, Rainfall
covariate <- readRDS("C:/Users/aelagali/OneDrive - The Peepingee Trust/Desktop/Yalem/month.covariate.rds") 
#saveRDS(covariate, file = "Z:/Yalem/Ethiopia_malaria/Ahmed_data/covariate_scaled.rds")  # unscaled

##Merge covariates with case data

d4_2019$month <-  d4_2019$Month
d4_2020$month <-  d4_2020$Month
d4_2021$month <-  d4_2021$Month

covariate$region <- covariate$ADM1_EN
covariate$zone <- covariate$ADM2_EN
covariate$year <- as.character(covariate$year)

d4_2019 <- right_join(d4_2019, covariate,by= c("month","year","region","zone","district"),all.x=TRUE, all.y=TRUE)# & ADM1_EN == region &  year  == year & month== Month))
d4_2020 <- right_join(d4_2020, covariate,by= c("month","year","region","zone","district"),all.x=TRUE, all.y=TRUE)# & ADM1_EN == region &  year  == year & month== Month))
d4_2021 <- right_join(d4_2021, covariate,by= c("month","year","region","zone","district"),all.x=TRUE, all.y=TRUE)# & ADM1_EN == region &  year  == year & month== Month))




d4_2019$year <- as.numeric(d4_2019$year)
d4_2020$year <- as.numeric(d4_2020$year)
d4_2021$year <- as.numeric(d4_2021$year)

d4_2019$district <- as.factor(d4_2019$district)
d4_2020$district <- as.factor(d4_2020$district)
d4_2021$district <- as.factor(d4_2021$district)

####Model
## Inference using INLA
#######################  INCIDENCE INTENSITY RATE #####################################################################
#' In any Bayesian analysis, the choice of the prior may have a considerable impact on the results.
#' By default, minimally informative priors are specified on:
#' (i) the log of the unstructured effect precision log(????) ∼ logGamma33(1, 0.0005), and
#' (ii) the log of the structured effect precision log(????u) ∼ logGamma(1, 0.0005).
#' However, for intrinsic models such as iCAR, the precision matrix depending on the random parameter and the neighboring structure)

d4_2019$yearm[d4_2019$month=="Jan"] <- 2019.083
d4_2019$yearm[d4_2019$month=="Feb"] <- 2019.167
d4_2019$yearm[d4_2019$month=="Mar"] <- 2019.250
d4_2019$yearm[d4_2019$month=="Apr"] <- 2019.333 
d4_2019$yearm[d4_2019$month=="May"] <- 2019.417
d4_2019$yearm[d4_2019$month=="Jun"] <- 2019.500
d4_2019$yearm[d4_2019$month=="Jul"] <- 2019.583
d4_2019$yearm[d4_2019$month=="Aug"] <- 2019.667
d4_2019$yearm[d4_2019$month=="Sep"] <- 2019.750
d4_2019$yearm[d4_2019$month=="Oct"] <- 2019.833
d4_2019$yearm[d4_2019$month=="Nov"] <- 2019.917
d4_2019$yearm[d4_2019$month=="Dec"] <- 2020.000

d4_2020$yearm[d4_2020$month=="Jan"] <- 2020.083
d4_2020$yearm[d4_2020$month=="Feb"] <- 2020.167
d4_2020$yearm[d4_2020$month=="Mar"] <- 2020.250
d4_2020$yearm[d4_2020$month=="Apr"] <- 2020.333 
d4_2020$yearm[d4_2020$month=="May"] <- 2020.417
d4_2020$yearm[d4_2020$month=="Jun"] <- 2020.500
d4_2020$yearm[d4_2020$month=="Jul"] <- 2020.583
d4_2020$yearm[d4_2020$month=="Aug"] <- 2020.667
d4_2020$yearm[d4_2020$month=="Sep"] <- 2020.750
d4_2020$yearm[d4_2020$month=="Oct"] <- 2020.833
d4_2020$yearm[d4_2020$month=="Nov"] <- 2020.917
d4_2020$yearm[d4_2020$month=="Dec"] <- 2021.000

d4_2021$yearm[d4_2021$month=="Jan"] <- 2021.083
d4_2021$yearm[d4_2021$month=="Feb"] <- 2021.167
d4_2021$yearm[d4_2021$month=="Mar"] <- 2021.250
d4_2021$yearm[d4_2021$month=="Apr"] <- 2021.333 
d4_2021$yearm[d4_2021$month=="May"] <- 2021.417
d4_2021$yearm[d4_2021$month=="Jun"] <- 2021.500
d4_2021$yearm[d4_2021$month=="Jul"] <- 2021.583
d4_2021$yearm[d4_2021$month=="Aug"] <- 2021.667
d4_2021$yearm[d4_2021$month=="Sep"] <- 2021.750
d4_2021$yearm[d4_2021$month=="Oct"] <- 2021.833
d4_2021$yearm[d4_2021$month=="Nov"] <- 2021.917
d4_2021$yearm[d4_2021$month=="Dec"] <- 2021.000



d4_tot=rbind(d4_2019,d4_2020,d4_2021)
################################################OPTION 1##############################################################
###The model assume a non-parametric model for the time trend 
##create the index vectors for the counties and years that will be used to specify the random effects of the model
d4_tot$idtime <- 1 + d4_tot$yearm - min(d4_tot$yearm)
d4_tot$idarea <- as.numeric(d4_tot$district) # is the vector with the indices of districts 1 to 967
d4_tot$idarea1 <- d4_tot$idarea
d4_tot$idtime1 <- d4_tot$idtime# is the vector with the indices of years and has 1 to 3 number of years
d4_tot$ID.area.year = seq(1,length(d4_tot$district)) #interaction

##Assign weights to the neighbors 
nb2 <- poly2nb(map.eth)
head(nb2)

nb2INLA("ethiopi.map.adj", nb2)

# equal neighboring polygon will be assigned equal weight when computing the neighboring mean income values.
lw <- nb2listw(nb2, style="W", zero.policy=TRUE)
##To see the weight of the first polygon’s neighbors type:
lw$weights[1]

##read the graph
eth.map2.graph <- inla.read.graph(filename = "ethiopi.map.adj")
##Adjacency matrix: rows and columns identify areas; dotes identify neighbors
image(inla.graph2matrix(eth.map2.graph), xlab = "", ylab = "")





f_besag2<- confirmed_all ~  EVI + LST_Day + LST_Night + TCB + TCW + Rainfall+
  f(idarea,  model = 'besag',  graph = eth.map2.graph) + 
  f(idarea1, idtime, model = 'iid', constr=TRUE) + idtime


######   Yalem to measure the posterior exceeeding probability you need to return return.marginals.predictor=TRUE
m_besag1<- inla(f_besag2,
                family = "poisson", 
                data = d4_tot, E = population,
                control.predictor = list(compute = TRUE),
                control.compute=list(dic=TRUE,config = TRUE,
                                     cpo=TRUE,  waic = TRUE,return.marginals.predictor=TRUE),
                control.inla = list(int.strategy = "eb"),
                verbose = TRUE)




##Output
summary(m_besag1)


d4_tot$RR <- exp(log(m_besag1$summary.fitted.values[, "mean"]))
d4_tot$RRL <- exp(log(m_besag1$summary.fitted.values[, "0.025quant"]))
d4_tot$RRUL <- exp(log(m_besag1$summary.fitted.values[, "0.975quant"]))





