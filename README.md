# Bayesian Spatio-Temporal Analysis of Malaria Risk in Ethiopia 🦟🗺️

[![Language](https://img.shields.io/badge/Language-R-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Framework](https://img.shields.io/badge/Framework-R--INLA-blue)](https://www.r-inla.org/)
[![Topic](https://img.shields.io/badge/Topic-Spatial%20Epidemiology-brightgreen)]()
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## 📖 About

This repository contains R code for **Bayesian spatio-temporal modelling of malaria incidence** across districts of **Ethiopia (2019–2021)** using the [R-INLA](https://www.r-inla.org/) (Integrated Nested Laplace Approximations) framework.

The analysis models confirmed malaria cases at the district level by combining:
- **Spatial random effects** (Besag-York-Mollié / CAR models) to capture geographic clustering
- **Temporal random effects** (random walk models) to capture year-to-year trends
- **Space–time interaction terms** to capture heterogeneous local trends
- **Environmental and intervention covariates** to identify key drivers of malaria transmission

> ⚠️ Note: Data from Tigray region for 2021 are excluded due to conflict-related reporting gaps.

---

## 🎯 Objectives

1. Estimate district-level malaria **incidence intensity rates (IIR)** per 1,000 population across Ethiopia.
2. Identify environmental and programmatic **predictors** of malaria risk.
3. Produce **uncertainty-quantified risk maps** (posterior mean, lower/upper credible intervals).
4. Calculate **posterior exceedance probabilities** to flag high-risk districts.

---

## 📦 Dependencies

Install the required R packages before running the analysis:

```r
# Core spatial packages
install.packages(c("sf", "sp", "spdep", "spData", "rgdal", "rgeos", "geoR"))

# Data wrangling & visualisation
install.packages(c("tidyverse", "hablar", "tmap", "viridis", "classInt", "ggpubr"))

# Multicollinearity diagnostics
install.packages(c("olsrr", "corrplot"))

# Bayesian spatial modelling (INLA)
install.packages("INLA", repos = "https://inla.r-inla-download.org/R/stable")
```

---

## 🗂️ Repository Structure

```
R_INLA/
├── INLA_code_final.R   # Main analysis script (data prep → EDA → modelling → mapping)
└── README.md
```

---

## ⚙️ Analysis Workflow

### 1. Data Preparation

```r
# Load Ethiopia district shapefile and case data
map.eth <- readOGR("edited_v3.shp", verbose = TRUE)
d4 <- readRDS("tes_case_year.rds") |>
  rename("year" = "Year") |>
  mutate(
    ir  = (confirmed_all / population) * 1000,  # Incidence rate per 1,000
    tp  = (confirmed_all / test_performed) * 100, # Test positivity rate
    Ex  = r * population,                         # Expected cases
    sir = confirmed_all / Ex                      # Standardised Incidence Ratio
  )
```

### 2. Spatial Autocorrelation (Moran's I)

```r
# Build spatial weights and test for clustering
nb2 <- poly2nb(map.eth, queen = TRUE, row.names = map.eth$district)
lw  <- nb2listw(nb2, style = "W", zero.policy = TRUE)

moran.test(eth_eth9$confirmed_all, lw, alternative = "greater")
MC <- moran.mc(eth_eth9$confirmed_all, lw, nsim = 999, alternative = "greater")
plot(MC)
```

### 3. Covariate Selection (VIF / Correlation)

```r
library(olsrr)
my_model <- lm(confirmed_all ~ ., data = d_cov)
ols_vif_tol(my_model)   # Remove variables with VIF > 10

library(corrplot)
corrplot(cor(data_cor), method = "number")
```

### 4. Bayesian Spatio-Temporal Modelling with INLA

Eight nested models are compared (DIC / WAIC / CPO):

| Model | Components |
|-------|-----------|
| `model1` | Covariates only (Poisson) |
| `model2` | + BYM spatial random effect |
| `model3` | + RW2 temporal random effect |
| `model4` | + IID temporal random effect |
| `model5` ⭐ | + Space–time IID interaction (**best fit**) |
| `model6` | + IID spatial × RW2 temporal interaction |
| `model7` | + IID temporal × Besag spatial interaction |
| `model8` | + Besag spatial × RW2 temporal interaction |

```r
# Best-fitting model (Model 5) — space-time IID interaction
formula5 <- confirmed_all ~ Night + AccessHF + DTW + Urban + IRS +
  Elevation + Day + TCB + Rainfall + Slope +
  f(idarea,      model = 'bym', graph = g) +
  f(idtime,      model = 'rw2') +
  f(ID.area.year, model = 'iid')

model5 <- inla(
  formula5,
  family  = "poisson",
  data    = d4,
  E       = population,
  control.predictor = list(link = 1),
  control.compute   = list(dic = TRUE, cpo = TRUE, waic = TRUE)
)

summary(model5)
round(model5$summary.fixed[, 1:5], 3)
```

### 5. Risk Mapping

```r
# Posterior mean IIR map (per 1,000 population)
p_mean <- ggplot() +
  geom_sf(data = mapredict, aes(fill = pred_mean * 1000), colour = NA) +
  geom_sf(data = admin1_shp, colour = 'black', fill = NA) +
  scale_fill_viridis(option = "magma", direction = -1) +
  facet_wrap(~year, ncol = 3) +
  theme_bw()
```

### 6. Posterior Exceedance Probability

```r
# Probability that district-level IIR exceeds threshold (e.g., 0.07)
a2 <- 0.07
d4$inlaprob <- unlist(lapply(model5$marginals.fitted.values, function(X) {
  1 - inla.pmarginal(a2, X)
}))

d4$exceed.cut <- cut(
  d4$inlaprob,
  breaks = c(0, 0.1, 0.2, 0.5, 0.8, 0.95, 1),
  include.lowest = TRUE
)
```

---

## 📊 Key Covariates

| Category | Variables |
|----------|-----------|
| **Climate** | Rainfall, Day/Night land surface temperature, Tasseled Cap Brightness (TCB) |
| **Environment** | Elevation, Slope, Distance to water bodies (DTW), EVI |
| **Urbanisation** | Urban population ratio, VIIRS night-time lights |
| **Healthcare access** | Access to health facilities (AccessHF), Access to cities |
| **Intervention** | Indoor residual spraying (IRS), Insecticide-treated nets (ITN) |

---

## 📌 Key Outputs

- 🗺️ **Choropleth maps** of raw incidence rate and test positivity rate by year
- 📈 **Temporal trend plots** of the spatial and global random effects
- 🔬 **Fixed-effect tables** (posterior mean, SD, 95% credible intervals)
- ⚖️ **Model comparison** via DIC, WAIC, and CPO
- 🎯 **Exceedance probability maps** highlighting high-risk districts
- 📉 **Residual diagnostics** and PIT (Probability Integral Transform) histograms

---

## 📚 References

- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian models using integrated nested Laplace approximations. *Journal of the Royal Statistical Society: Series B*, 71(2), 319–392.
- Besag, J., York, J., & Mollié, A. (1991). Bayesian image restoration with two applications in spatial statistics. *Annals of the Institute of Statistical Mathematics*, 43, 1–20.
- Blangiardo, M., & Cameletti, M. (2015). *Spatial and Spatio-temporal Bayesian Models with R-INLA*. Wiley.

---

## ✉️ Author

**Yalemzewod Gelaw**  
Bayesian spatial epidemiologist | Malaria & infectious disease modelling  
🔗 [GitHub Profile](https://github.com/Yalemzewod)
