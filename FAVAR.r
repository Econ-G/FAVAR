  rm(list=ls())
  # Load required packages
  library(readxl)   
  library(FAVAR)    
  library(vars)     
  library(patchwork)  
  library(repr)
  library(boot)
  library(tsDyn)
  getwd()
  data_all <- read_excel("Your path",sheet = 1)
  
  # Remove the date column (assumed to be the first column)
  data_s = scale(data_all[,2:187], center = TRUE, scale = TRUE)
  colnames(data_s)

  # Regress Principal Components for all variables
  pc_all = prcomp(data_s, center=FALSE, scale.=FALSE, rank. = 3)
  C = pc_all$x # saving the principal components

  # Define slow-moving variables
  slow_vars_indices <- 1:114
  slow_vars <- colnames(data_all[,2:187])[slow_vars_indices]
  head(slow_vars)
  data_slow = data_s[,1:114]
  head(data_slow)

  # Regress factors for slow-moving variables
  pc_slow = prcomp(data_slow, center=FALSE, scale.=FALSE, rank. = 3)
  F_slow=pc_slow$x

  # Equation 3, extract comtemporenous impact from factors
  reg = lm(C ~ F_slow + data_s[,"log_ffr"])

  # Regress Fitted Factor
  F_hat = C - data.matrix(data_s[,"log_ffr"])%*%reg$coefficients[5,] 

  #Create Data frame for VAR
  data_var = data.frame(F_hat, "log_ffr" = data_s[,"log_ffr"])

  #Check the optimal lag
  lag_selection <- VARselect(data_var, lag.max = 15, type = "none")
  print(lag_selection$selection)
  #Run VAR Analysis
  var = VAR(data_var, p = 3)
  summary(var)
  irf_point = irf(var, n.ahead = 48, impulse = "log_ffr", response = "log_ffr", boot = FALSE)
  
  # Shock size of 25 basis points
  impulse_sd = 0.25/sd(data_all$log_ffr)
  scale = impulse_sd/(irf_point$irf$log_ffr[1]) # position of log_ffr response at step 0
  
  
  # Computing Loading Factors
  reg_loadings = lm(data_s ~ 0 + F_hat + data_s[,"log_ffr"])
  loadings = reg_loadings$coefficients
  head(reg_loadings$coefficients)
  summary(reg_loadings)
  
  
  #### BOOTSTRAPING ########
  
  R = 5000 # Number of simulations
  nvars = 186 # Number of variables
  nsteps = 49 # numbers of steps
  
  IRFs = array(c(0,0,0), dim = c(nsteps,nvars,R))
  
  var = lineVar(data_var, lag = 3, include = "const")
  for(j in 1:R){    
    data_boot = VAR.boot(var, boot.scheme ="resample")
    var_boot = VAR(data_boot, lag = 3)
    irf1 = irf(var_boot, n.ahead = 48, impulse = "log_ffr", boot = FALSE)
    for(i in 1:nvars){
      IRFs[,i,j] = (irf1$irf$log_ffr %*% matrix(loadings[1:4, i]))*scale
    }
  } ## Boot simulations done
  # Extract the quantiles of IRFs we are interested: 90% confidence intervals in BBE
  Upper = array(c(0,0), dim = c(nsteps, nvars))
  for(k in 1:nsteps){
    for(i in 1:nvars){
      Upper[k,i] = quantile(IRFs[k,i,], probs = c(0.90))[1]
    }
  }
  Lower = array(c(0,0), dim = c(nsteps, nvars))
  for(k in 1:nsteps){
    for(i in 1:nvars){
      Lower[k,i] = quantile(IRFs[k,i,], probs = c(0.10))[1]
    }
  }
  IRF = array(c(0,0), dim = c(nsteps, nvars))
  for(k in 1:nsteps){
    for(i in 1:nvars){
      IRF[k,i] = quantile(IRFs[k,i,], probs = c(0.5))[1]
    }
  }
  rm(var_boot)
  rm(IRFs)
  # Select the Variables you are Interested in
  # List of variables we are interested: FYFF, IP, CPI
  head(data_s)
  variables = c(grep("^Total_GDP$",colnames(data_s)), grep("^CBAS$",colnames(data_s)), grep("^export_all$",colnames(data_s)),grep("^import_all$",colnames(data_s)),
                grep("^SOG$",colnames(data_s)), grep("^CDA_L$",colnames(data_s)), grep("^RM$",colnames(data_s)), grep("^Overnight_rt$",colnames(data_s)),
                grep("^TCLH$",colnames(data_s)),grep("^CBP_CD$",colnames(data_s)),grep("^MT$",colnames(data_s)), 
                grep("^CI$",colnames(data_s)),grep("^TI$",colnames(data_s)),grep("^UR$",colnames(data_s)), 
                grep("^ER$",colnames(data_s)),grep("^EER$",colnames(data_s)))
  length(variables)
  transf_code = c(5,5,5,5,
                  5,5,5,1,
                  5,5,5,
                  5,5,1,
                  1,1)
  variable_names = c("Total GDP", "Canadian Bond All Sectors", "Export","Import",
                    "Sales of Goods Manufactured (shipments)", "Canadian Dollar Asset", "Residential Mortgage","Overnight Rate",
                    "Total Credit Liabilities of Households","Chartered Bank Deposits","Total Mortgage",
                    "Core Inflation", "Total Inflation","Unemployment Rate",
                    "Employment Rate","Effective Exchange Rate")


 
  # Change plot size to 15 x 10
  options(repr.plot.width=12, repr.plot.height=8)
  
  par(mfrow=c(4,4), 
      mar = c(2, 2, 2, 2))
  
  # 1) define your single theme color
  theme_col <- "#9976F6"
  
  # 2) build two lighter “band” colors from it
  band_light  <- adjustcolor(theme_col, alpha.f = 0.8)  # pale
  band_lighter <- adjustcolor(theme_col, alpha.f = 0.8) # extra pale
  
  # 3) updated plotting loop
  for (i in variables) {
    idx <- which(variables == i)
    
    if (transf_code[idx] == 5) {
      # cumulative IRF
      plot(
        cumsum(IRF[,i]), type = "l", lwd = 2,
        col      = theme_col,
        main     = variable_names[idx],
        col.main = theme_col,
        ylab     = "", xlab = "Steps",
        ylim     = range(cumsum(Lower[,i]), cumsum(Upper[,i])),
        cex.main = 1.8, cex.axis = 1.3
      )
      lines(cumsum(Upper[,i]), lty = 2, col = band_light)
      lines(cumsum(Lower[,i]), lty = 2, col = band_light)
      abline(h = 0, col = theme_col, lwd = 1)
      
    } else {
      # non‐cumulative IRF
      plot(
        IRF[,i], type = "l", lwd = 2,
        col      = theme_col,
        main     = variable_names[idx],
        col.main = theme_col,
        ylab     = "", xlab = "Steps",
        ylim     = range(Lower[,i], Upper[,i]),
        cex.main = 1.8, cex.axis = 1.3
      )
      lines(Upper[,i], lty = 2, col = band_light)
      lines(Lower[,i], lty = 2, col = band_light)
      abline(h = 0, col = theme_col, lwd = 1)
    }
  }
  options(repr.plot.width=12, repr.plot.height=6)

  # Get the VAR point estimates
  hor = 60
  var = VAR(data_var, p = 3)
  irf_point = irf(var, n.ahead = hor, boot = FALSE)
  # Get IRFs for all of X we are interested in, Dimensions: (hor, key_nvars)
  # Find loadings
  results = summary(reg_loadings) # the warning comes because of FYFF
  
  key_nvars = length(variables)
  irf_X_pc1 = array(c(0,0), dim=c(hor+1, key_nvars))
  irf_X_pc2 = array(c(0,0), dim=c(hor+1, key_nvars))
  irf_X_pc3 = array(c(0,0), dim=c(hor+1, key_nvars))
  irf_X_log_ffr = array(c(0,0), dim=c(hor+1, key_nvars))
  
  for(i in 1:key_nvars){
    irf_X_pc1[,i] = irf_point$irf$PC1 %*% matrix(loadings[1:4, variables[i]])
    irf_X_pc2[,i] = irf_point$irf$PC2 %*% matrix(loadings[1:4, variables[i]])
    irf_X_pc3[,i] = irf_point$irf$PC3 %*% matrix(loadings[1:4, variables[i]])
    irf_X_log_ffr[,i] = (irf_point$irf$log_ffr) %*% matrix(loadings[1:4, variables[i]])
  }
  # Get the IRFs squared and accumulate them
  psi2_pc1 = array(0, dim = key_nvars)
  psi2_pc2 = array(0, dim = key_nvars)
  psi2_pc3 = array(0, dim = key_nvars)
  psi2_log_ffr= array(0, dim = key_nvars)
  
  for(i in 1:key_nvars){
    for(j in 1:hor){
      psi2_pc1[i] = psi2_pc1[i] + irf_X_pc1[j,i]^2
      psi2_pc2[i] = psi2_pc2[i] + irf_X_pc2[j,i]^2
      psi2_pc3[i] = psi2_pc3[i] + irf_X_pc3[j,i]^2
      psi2_log_ffr[i] = psi2_log_ffr[i] + irf_X_log_ffr[j,i]^2    
    }
  }
  var_total= array(0, dim = key_nvars)
  var_fac= array(0, dim = key_nvars)
  var_e= array(0, dim = key_nvars)
  
  for(i in 1:key_nvars){
    var_fac[i] = psi2_pc1[i] + psi2_pc2[i] + psi2_pc3[i] + psi2_log_ffr[i]
    var_total[i] = psi2_pc1[i] + psi2_pc2[i] + psi2_pc3[i] + psi2_log_ffr[i] + results[[variables[i]]]$sigma^2
    var_e[i] = results[[variables[i]]]$sigma^2
  }
  table = data.frame("PC1" = round((psi2_pc1),3), "PC2" = round((psi2_pc2),3), "PC3" = round((psi2_pc3),3), "FFR" = round((psi2_log_ffr),3),
                     "Factor_Y_total" = round(var_fac,3) ,"e" = round((var_e),3), "Total" = round(var_total,3))
  library(stargazer)
  row.names(table) = variable_names
  table
  print(variable_names)
  capture.output(print(table), file = "table.txt")
  # Replicating Table I in BBE (2005) - 3 Factors and Y = FYFF
  r2 = array(0, dim = key_nvars)
  for(i in 1:key_nvars){
    r2[i] = results[[variables[i]]]$r.squared
  }
  table2 = data.frame("Variables" = variable_names, "Contribution" = round((psi2_log_ffr/var_total),3), "R-squared" = round(r2,3))

  
  table2
  capture.output(print(table2), file = "table2.txt")
  
