library(dplyr)
# library(lubridate)
library(fastmatrix)
library(graphicalVAR)

source("src/FUN/helpers.R")

# Dataframe to gather estimates parameters
df_param = data.frame()
list_param = list()






###################################
############ Epskamp 2018
###################################

# Load the data:
Data <- read.csv("data/raw/epskampsupplementary2_data.csv",header=TRUE, stringsAsFactors = FALSE)

# Variables to include in analysis:
Vars <- c("relaxed","sad","nervous","concentration","tired","rumination","bodily.discomfort")

# Time varible:
Data$time <- as.POSIXct(Data$time)

# Day variable:
Data$date <- as.numeric(as.Date(Data$time))

# Detrending significant linear trends:
for (v in seq_along(Vars)){
  ff <- as.formula(paste0(Vars[[v]]," ~ time"))
  fit <- lm(ff, data = Data)
  if (anova(fit)$P[1] < 0.05){
    message(paste("Detrending variable",v))
    Data[[Vars[v]]][!is.na(Data[[Vars[[v]]]])] <- residuals(fit)
  }
}

# First prepare a dataset with lagged values
Data_var <- Data %>%
  mutate(intercept = rep(1, 70),
         relaxed = scale(relaxed),
         sad = scale(sad),
         nervous = scale(nervous),
         concentration = scale(concentration),
         tired = scale(tired),
         rumination = scale(rumination),
         bodily.discomfort = scale(bodily.discomfort)) %>%
  mutate(relaxed_l1 = lag(relaxed),
         sad_l1 = lag(sad),
         nervous_l1 = lag(nervous),
         concentration_l1 = lag(concentration),
         tired_l1 = lag(tired),
         rumination_l1 = lag(rumination),
         bodily.discomfort_l1 = lag(bodily.discomfort)) %>%
  mutate(firstbeep = rep(c(1, 0, 0, 0, 0), 14)) %>%
  na.omit() %>%
  filter(firstbeep == 0) %>%
  dplyr::select(-c(time, date, firstbeep))

Res <- graphicalVAR(Data, gamma = 0, vars = Vars, dayvar = "date")

# Now calculate the in-sample residual terms
residual_mat <- as.matrix(Data_var[,1:7], ncol = 7) - t(as.matrix(Res$beta, ncol = 8) %*% t(as.matrix(Data_var[, 8:15], ncol = 8)))


VAR_phi <- as.matrix(Res$beta[,-1], ncol = 7)
VAR_delta <- as.matrix(Res$beta[,1], ncol = 1)
VAR_sigma <- var(residual_mat)



list_param[['epskamp']] <- list("vars" = 7, "delta" = VAR_delta, "phi" = VAR_phi, "sigma" = VAR_sigma)
df_param <- bind_rows(df_param, estimates_to_df("epskamp", VAR_delta, VAR_phi, VAR_sigma))




###################################
############ Bak2016 Impending relapse
###################################

# Estimates acquired from Marjan Drukker
int_mat <- as.vector(c(2.604700, 
                       2.209149, 
                       1.639116,
                       3.175556,
                       1.130181))

reg_mat <- matrix(c(0.08951318, -0.23055938, 0.2978268, 0.06635082, -0.06652466,
                    0.29872946, 0.34655150, -0.1367041, -0.03469113, 0.22448471,
                    0.36220131, -0.21625369, 0.2106809, 0.16662747, -0.20135875,
                    -0.15787115, 0.05234745, 0.4998060, 0.07528141, -0.11408725,
                    0.14749443, -0.10265810, 0.1986619, 0.01222618, 0.04620206),
                 nrow = 5, byrow = T)

SE_mat <- matrix(c(0.1571154, 0.1372617, 0.12022895, 0.11406118, 0.1745014,
                   0.1292835, 0.1129467, 0.09893123, 0.09385604, 0.1435897,
                   0.1888545, 0.1649901, 0.14451658, 0.13710285, 0.2097526,
                   0.1547560, 0.1352004, 0.11842348, 0.11234833, 0.1718809,
                   0.1362937, 0.1190711, 0.10429562, 0.09894523, 0.1513756),
                 nrow = 5, byrow = T)

p_mat <- matrix(c(0.57047809, 0.096965160, 0.0153780350, 0.5624166, 0.7040581,
                  0.02346124, 0.002949086, 0.1709255236, 0.7126541, 0.1219613,
                  0.05873932, 0.193754355, 0.1488507784, 0.2278536, 0.3399933,
                  0.31078092, 0.699661047, 0.0000645791, 0.5047663, 0.5087779,
                  0.28246428, 0.391210086, 0.0604461771, 0.9019735, 0.7610053),
                nrow = 5, byrow = T)

pos_neg_mat <- 2*((reg_mat > 0) - .5)

eigen_1SElarger <- as.numeric(eigen(reg_mat + hadamard(pos_neg_mat, SE_mat))$values[1])

# Rescale the +1SE phi matrix if its eigenvalue is larger than 1
reg_mat_1SElarger <- (reg_mat + hadamard(pos_neg_mat, SE_mat)) * ifelse(eigen_1SElarger < 1, 1, .99/eigen_1SElarger)

Mod(eigen(reg_mat_1SElarger)$values)

reg_mat_1SEsmaller <- ifelse(hadamard(reg_mat, reg_mat - hadamard(pos_neg_mat, SE_mat))<0, 0, reg_mat - hadamard(pos_neg_mat, SE_mat))

cov_mat <- matrix(c(2.0309436, -0.47198788, 1.4410971, 0.74049942, 0.7134376,
                    -0.4719879, 1.37513898, -0.3968026, -0.05529972, 0.1656011,
                    1.4410971, -0.39680261, 2.9343717, 0.89314648, 0.4681241,
                    0.7404994, -0.05529972, 0.8931465, 1.97040458, 0.2578435,
                    0.7134376, 0.16560106, 0.4681241, 0.25784346, 1.5283116),
                  nrow = 5, byrow = T)


# Estimates
list_param[['bak2016_relapse']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_relapse", int_mat, reg_mat, cov_mat))

# Estimates at +1 se
list_param[['bak2016_relapse_plus1se']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat_1SElarger, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_relapse_plus1se", int_mat, reg_mat_1SElarger, cov_mat))

# Estimates at-+1 se
list_param[['bak2016_relapse_minus1se']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat_1SEsmaller, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_relapse_minus1se", int_mat, reg_mat_1SEsmaller, cov_mat))





###################################
############ Bak2016 full relapse
###################################

# Estimates acquired from Marjan Drukker
int_mat <- as.vector(c(3.0307550, 
                       3.0599007, 
                       1.5954037,
                       2.3062241,
                       0.7724772))

reg_mat <- matrix(c(0.18787336, -0.09466088, 0.4441990, -0.34945708, 0.01713231,
                    0.07034037, 0.17178755, -0.3351465, 0.09114686, -0.04456589,
                    -0.09401958, -0.24303723, 0.3484971, 0.23436851, 0.13096188,
                    -0.04573019, 0.08841517, 0.6017965, 0.11153617, -0.10839307,
                    0.32784438, 0.06905052, 0.3980946, -0.15269345, 0.14769753),
                 nrow = 5, byrow = T)

SE_mat <- matrix(c(0.1524883, 0.2237966, 0.1837789, 0.1420838, 0.13197875,
                   0.1117048, 0.1639414, 0.1346266, 0.1040830, 0.09668057,
                   0.1651632, 0.2484964, 0.1930229, 0.1490116, 0.13895543,
                   0.1542511, 0.2331640, 0.1810691, 0.1402464, 0.13042514,
                   0.1703194, 0.2574527, 0.1999311, 0.1548559, 0.14401151),
                 nrow = 5, byrow = T)

p_mat <- matrix(c(0.22307827, 0.6739317, 0.018929876, 0.01702647, 0.8971810,
                  0.53145245, 0.2992049, 0.015790701, 0.38492659, 0.6466125,
                  0.57150231, 0.3323414, 0.076475331, 0.12149600, 0.3500716,
                  0.76801165, 0.7060284, 0.001600725, 0.42992845, 0.4095910,
                  0.05951908, 0.7895622, 0.051532795, 0.32851316, 0.3096548),
                nrow = 5, byrow = T)

pos_neg_mat <- 2*((reg_mat > 0) - .5)

eigen_1SElarger <- as.numeric(eigen(reg_mat + hadamard(pos_neg_mat, SE_mat))$values[1])

# Rescale the +1SE phi matrix if its eigenvalue is larger than 1
reg_mat_1SElarger <- (reg_mat + hadamard(pos_neg_mat, SE_mat)) * ifelse(eigen_1SElarger < 1, 1, .99/eigen_1SElarger)

Mod(eigen(reg_mat_1SElarger)$values)

reg_mat_1SEsmaller <- ifelse(hadamard(reg_mat, reg_mat - hadamard(pos_neg_mat, SE_mat))<0, 0, reg_mat - hadamard(pos_neg_mat, SE_mat))

cov_mat <- matrix(c(1.8408029, -0.5822931, 0.3546145, -0.13246672, -0.18202570,
                    -0.5822931, 0.9878197, -0.4677148, -0.17506083, -0.26498470,
                    0.3546145, -0.4677148, 2.1931548, -0.30282461, 0.45406298,
                    -0.1324667, -0.1750608, -0.3028246, 1.84308908, 0.03489562,
                    -0.1820257, -0.2649847, 0.4540630, 0.03489562, 2.31123798),
                  nrow = 5, byrow = T)

# Export
list_param[['bak2016_full']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_full", int_mat, reg_mat, cov_mat))

# Estimates at +1 se
list_param[['bak2016_full_plus1se']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat_1SElarger, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_full_plus1se", int_mat, reg_mat_1SElarger, cov_mat))

# Estimates at-+1 se
list_param[['bak2016_full_minus1se']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat_1SEsmaller, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_full_minus1se", int_mat, reg_mat_1SEsmaller, cov_mat))





###################################
############ Bak2016 stable
###################################

# Estimates acquired from Marjan Drukker
int_mat <- as.vector(c(1.8298033, 
                       2.8341405, 
                       2.1711898,
                       4.1862781,
                       0.6600582))

reg_mat <- matrix(c(0.29352932, -0.16396087, 0.160042647, -0.055128536, 0.08153145,
                    -0.03618744, 0.39163342, -0.144333033, 0.001253406, 0.02120229,
                    0.27018456, -0.21596785, 0.237039617, -0.107376167, 0.06191572,
                    0.08420763, -0.15027115, 0.056606422, 0.170175697, -0.02705356,
                    0.19325018, -0.08521795, 0.004444648, 0.054694766, 0.11197789),
                 nrow = 5, byrow = T)

SE_mat <- matrix(c(0.05168159, 0.05514216, 0.04589414, 0.04956212, 0.06082448,
                   0.04613447, 0.04922360, 0.04096820, 0.04424248, 0.05429603,
                   0.06387909, 0.06815639, 0.05672573, 0.06125939, 0.07517982,
                   0.05542543, 0.05913668, 0.04921873, 0.05315241, 0.06523063,
                   0.05129353, 0.05472811, 0.04554954, 0.04918997, 0.06036777),
                 nrow = 5, byrow = T)

p_mat <- matrix(c(2.776127e-08, 3.142271e-03, 5.481198e-04, 0.266742164, 0.18094359,
                  4.333240e-01, 2.298778e-14, 4.812252e-04, 0.977414270, 0.69640081,
                  2.968770e-05, 1.661471e-03, 3.679403e-05, 0.080480884, 0.41072712,
                  1.295609e-01, 1.146768e-02, 2.508622e-01, 0.001487345, 0.67858007,
                  1.923942e-04, 1.203170e-01, 9.223211e-01, 0.266915228, 0.06441933),
                nrow = 5, byrow = T)

pos_neg_mat <- 2*((reg_mat > 0) - .5)

library(fastmatrix)
eigen_1SElarger <- as.numeric(eigen(reg_mat + hadamard(pos_neg_mat, SE_mat))$values[1])

# Rescale the +1SE phi matrix if its eigenvalue is larger than 1
reg_mat_1SElarger <- (reg_mat + hadamard(pos_neg_mat, SE_mat)) * ifelse(eigen_1SElarger < 1, 1, .99/eigen_1SElarger)

Mod(eigen(reg_mat_1SElarger)$values)

reg_mat_1SEsmaller <- ifelse(hadamard(reg_mat, reg_mat - hadamard(pos_neg_mat, SE_mat))<0, 0, reg_mat - hadamard(pos_neg_mat, SE_mat))

cov_mat <- matrix(c(1.8679824, -0.3608477, 0.5418325, 0.4365137, 0.3071496,
                    -0.3608477, 1.4885110, -0.3513998, -0.3334411, -0.1629183,
                    0.5418325, -0.3513998, 2.8537668, 0.8841644, 0.5358375,
                    0.4365137, -0.3334411, 0.8841644, 2.1484194, 0.2967760,
                    0.3071496, -0.1629183, 0.5358375, 0.2967760, 1.8400356),
                  nrow = 5, byrow = T)



# Export
list_param[['bak2016_stable']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_stable", int_mat, reg_mat, cov_mat))

# Estimates at +1 se
list_param[['bak2016_stable_plus1se']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat_1SElarger, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_stable_plus1se", int_mat, reg_mat_1SElarger, cov_mat))

# Estimates at-+1 se
list_param[['bak2016_stable_minus1se']] <- list("vars" = 5, "delta" = as.matrix(int_mat), "phi" = reg_mat_1SEsmaller, "sigma" = cov_mat, "p_mat"=p_mat)
df_param <- bind_rows(df_param, estimates_to_df("bak2016_stable_minus1se", int_mat, reg_mat_1SEsmaller, cov_mat))




###################################
############ Illustration 
###################################

VAR_phi <- matrix(c(.5, .2, .3, .4), nrow = 2, byrow = F)
VAR_delta <- matrix(c(0, 2), nrow = 2, byrow = F)
VAR_sigma <- matrix(c(10, 4, 4, 10), nrow = 2, byrow = F)

N_variable <- nrow(VAR_phi)


list_param[['illustration']] <- list("vars" = N_variable, "delta" = VAR_delta, "phi" = VAR_phi, "sigma" = VAR_sigma)
df_param <- bind_rows(df_param, estimates_to_df("illustration", VAR_delta, VAR_phi, VAR_sigma))



###################################
############ Export everything in Robject and csv
###################################
# save(list_param, file="data/list_param.rda")
# write.csv(df_param, "data/df_param.csv", row.names=FALSE)



