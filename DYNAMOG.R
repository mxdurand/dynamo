#### Original: G is for stomatal conductance
DYNAMOG <- function(df = df, id = NULL, norm = FALSE, timestep = 30, PltLength = 100, ITER = 100, relPlt = NULL, plotall = TRUE, savePlot = FALSE)
{
  allresD <- data.frame("id" = character(0), "Tau" = numeric(0), "k" = numeric(0), "l" = numeric(0),
                        "G_EXP" = numeric(0), "g0_EXP" = numeric(0), "G_SIG" = numeric(0), "g0_SIG" = numeric(0),
                        "R2_EXP" = numeric(0), "R2_SIG" = numeric(0), "R2_EXP_ALL" = numeric(0), "R2_SIG_ALL" = numeric(0), stringsAsFactors=FALSE)
  
  allresU <- data.frame("id" = character(0), "Tau" = numeric(0), "k" = numeric(0), "l" = numeric(0),
                        "G_EXP" = numeric(0), "g0_EXP" = numeric(0), "G_SIG" = numeric(0), "g0_SIG" = numeric(0),
                        "R2_EXP" = numeric(0), "R2_SIG" = numeric(0), "R2_EXP_ALL" = numeric(0), "R2_SIG_ALL" = numeric(0), stringsAsFactors=FALSE)
  
  mresD <- data.frame("id" = character(0), "N_EXP" = integer(0), "N_SIG" = integer(0), "mTau" = numeric(0), "sdTau" = numeric(0), "mk" = numeric(0), "sdk" = numeric(0),
                      "ml" = numeric(0), "sdl" = numeric(0), "LC" = numeric(0), "mlt" = numeric(0), "mg0_EXP" = numeric(0), "sdg0_EXP" = numeric(0),
                      "mG_EXP" = numeric(0), "sdG_EXP" = numeric(0), "mg0_SIG" = numeric(0), "sdg0_SIG" = numeric(0),
                      "mG_SIG" = numeric(0), "sdG_SIG" = numeric(0), "R2quant_EXP" = numeric(0), "R2quant_SIG" = numeric(0), "R2quant_EXP_ALL" = numeric(0), "R2quant_SIG_ALL" = numeric(0),
                      "pmax_EXP" = numeric(0), "pmax_SIG" = numeric(0), "sdpmax_SIG" = numeric(0), "int_EXP" = numeric(0), "int_SIG" = numeric(0),
                      "A0" = numeric(0), "A1" = numeric(0), stringsAsFactors=FALSE)
  
  mresU <- data.frame("id" = character(0), "N_EXP" = integer(0), "N_SIG" = integer(0), "mTau" = numeric(0), "sdTau" = numeric(0), "mk" = numeric(0), "sdk" = numeric(0),
                      "ml" = numeric(0), "sdl" = numeric(0), "LC" = numeric(0), "mlt" = numeric(0), "mg0_EXP" = numeric(0), "sdg0_EXP" = numeric(0),
                      "mG_EXP" = numeric(0), "sdG_EXP" = numeric(0), "mg0_SIG" = numeric(0), "sdg0_SIG" = numeric(0),
                      "mG_SIG" = numeric(0), "sdG_SIG" = numeric(0), "R2quant_EXP" = numeric(0), "R2quant_SIG" = numeric(0), "R2quant_EXP_ALL" = numeric(0), "R2quant_SIG_ALL" = numeric(0),
                      "pmax_EXP" = numeric(0), "pmax_SIG" = numeric(0), "sdpmax_SIG" = numeric(0), "int_EXP" = numeric(0), "int_SIG" = numeric(0),
                      "A0" = numeric(0), "A1" = numeric(0), stringsAsFactors=FALSE)
  
  fit_EXP <- function(x){
    Tau <- x[1] # Time constant
    G <- x[2] # Target GS
    g0 <- x[3] # Initial GS
    EXPfit <- G+(g0-G)*exp(-(secondes_EXP/Tau))  # k en secondes et b represent temps au poitns d'inflection
    nash <- (sum((obs_EXP-EXPfit)^2, na.rm = T) / sum((obs_EXP-mean(obs_EXP, na.rm = T))^2, na.rm = T))
    return(nash)
  }
  
  fit_SIG <- function(x){
    k <- x[1]
    l <- x[2]
    G <- x[3]
    g0 <- x[4]
    #gmin <- Gstart
    #gmax <- Gend
    SIGfit <- g0+(G-g0)*exp(-exp((l-secondes_SIG)/k))  # k en secondes et b represent temps au poitns d'inflection
    nash <- (sum((obs_SIG-SIGfit)^2, na.rm = T) / sum((obs_SIG-mean(obs_SIG, na.rm = T))^2, na.rm = T))
    return(nash)
  }
  
  i = unique(df[,id])[1]
  for (i in unique(df[,id]))
  {
    if (savePlot == TRUE) {
      png(paste0("Fits/CurveFit_DYNAMOG_", i, ".png"), height = 8, width = 12, units = "in", res = 300)
    }
    
    # For later
    if (plotall == TRUE){
      split.screen(c(2,2))
    }
    
    # Subset data
    idColNo <- which(colnames(df) == id)
    sub.i <-  df[df[,idColNo] == i,]
    
    # Get down et up curves
    down.i <- sub.i[sub.i$DownG == "SS0" | sub.i$DownG == "DYN" | sub.i$DownG == "SS1",]
    up.i <- sub.i[sub.i$UpG == "SS0" | sub.i$UpG == "DYN" | sub.i$UpG == "SS1",]
    
    # Remove NAs
    down.i <- down.i[!is.na(down.i$DownG),]
    up.i <- up.i[!is.na(up.i$UpG),]
    
    # Check if data is present, otherwise go to next curve
    WhatSens <- c("DOWN", "UP")
    if(length(down.i$Cond)==0) {WhatSens <- "UP"}
    if(length(up.i$Cond)==0) {WhatSens <- "DOWN"}
    # if(sum(is.na(sub.i$Cond)) > 0) {next}
    cat("Started curve number",i, "|", WhatSens,"\n")
    
    # Run the thing twice (one for up and one for down)
    SENS = WhatSens[1]
    for (SENS in WhatSens)
    {
      if (SENS == "DOWN") {
        # If curve goes down
        start <- mean(down.i[down.i$DownG == "SS0",]$Cond, na.rm=T)
        end <- mean(down.i[down.i$DownG == "SS1",]$Cond, na.rm=T)
        curve <- down.i[down.i$DownG == "DYN",]
        A0d <- mean(down.i[down.i$DownG == "SS0",]$Photo, na.rm=T)
        A1d <- mean(down.i[down.i$DownG == "SS1",]$Photo, na.rm=T)
        
        # If curve normalized between 0 and 1
        if(norm == TRUE)
        {
          curve$Cond <- (curve$Cond - end) / (start - end)
          start <- 1
          end <- 0
        }
      } else if (SENS == "UP") {
        # If curve goes up
        start <- mean(up.i[up.i$UpG == "SS0",]$Cond, na.rm=T)
        end <- mean(up.i[up.i$UpG == "SS1",]$Cond, na.rm=T)
        curve <- up.i[up.i$UpG == "DYN",]
        A0u <- mean(up.i[up.i$UpG == "SS0",]$Photo, na.rm=T)
        A1u <- mean(up.i[up.i$UpG == "SS1",]$Photo, na.rm=T)
        
        # If curve normalized between 0 and 1
        if(norm == TRUE)
        {
          curve$Cond <- (curve$Cond - start) / (end - start)
          start <- 0
          end <- 1
        }
      }
      
      # Manipulate data to increase plateau length and adjust time
      # Exponential shouldn't have a starting plateau
      if (is.null(relPlt) == FALSE){
        PltLength = round(relPlt * nrow(curve))
      }
      
      
      CURVE_EXP <- c(start, curve$Cond, rep(end, times=PltLength))
      CURVE_SIG <- c(rep(start, times=PltLength), curve$Cond, rep(end, times=PltLength))
      
      timeSS0_EXP <- 0
      timeSS0_SIG <- seq(from=0, to=(timestep*PltLength-1), by=timestep)
      
      timeCurve_EXP <- (curve$Time-curve$Time[1]) + tail(timeSS0_EXP, n=1) + timestep
      timeCurve_SIG <- (curve$Time-curve$Time[1]) + tail(timeSS0_SIG, n=1) + timestep
      
      timeSS1_EXP <- seq(from=tail(timeCurve_EXP, n=1)+timestep, to=tail(timeCurve_EXP, n=1) + (timestep*PltLength), by=timestep)
      timeSS1_SIG <- seq(from=tail(timeCurve_SIG, n=1)+timestep, to=tail(timeCurve_SIG, n=1) + (timestep*PltLength), by=timestep)
      
      TIME_EXP <- c(timeSS0_EXP, timeCurve_EXP, timeSS1_EXP)
      TIME_SIG <- c(timeSS0_SIG, timeCurve_SIG, timeSS1_SIG)
      
      obs_EXP <- CURVE_EXP
      obs_SIG <- CURVE_SIG
      secondes_EXP <- TIME_EXP
      secondes_SIG <- TIME_SIG
      maxTime_EXP <- max(TIME_EXP)
      maxTime_SIG <- max(TIME_SIG)
      
      # Define light change based of direction of curve (otherwise it messes up later)
      if (SENS == "DOWN") {
        # If curve goes down
        LCD <- timestep*PltLength
        secondes_EXPD <- secondes_EXP
        secondes_SIGD <- secondes_SIG
      } else if (SENS == "UP") {
        # If curve goes up
        LCU <- timestep*PltLength
        secondes_EXPU <- secondes_EXP
        secondes_SIGU <- secondes_SIG
      }
      
      # Plot data
      if (plotall == TRUE){
        if (SENS == "DOWN") {
          # If curve goes down
          screen(1)
          plot(obs_SIG~secondes_SIG, main=paste("SIG | ", i, "|", SENS), xlim=c((timestep*PltLength)-500, (timestep*PltLength)+(nrow(curve)*timestep)+500))
          screen(2)
          plot(obs_EXP~secondes_EXP, main=paste("EXP | ", i, "|", SENS), xlim=c(0, (nrow(curve)*timestep)+500))
        } else if (SENS == "UP") {
          # If curve goes up
          screen(3)
          plot(obs_SIG~secondes_SIG, main=paste("SIG | ", i, "|", SENS), xlim=c((timestep*PltLength)-500, (timestep*PltLength)+(nrow(curve)*timestep)+500))
          screen(4)
          plot(obs_EXP~secondes_EXP, main=paste("EXP | ", i, "|", SENS), xlim=c(0, (nrow(curve)*timestep)+500))
        }
      }
      
      # fit curve to model
      for (it in 1:ITER)
      {
        # Create starting values (random for constants)
        startTau <- runif(1, 0, maxTime_EXP)
        startk <- runif(1, 0, maxTime_SIG)
        startl <- runif(1, 0, maxTime_SIG)
        startG <- end
        startg0 <- start
        #cat(it, "|", startTau,startG,startg0, "\n")
        
        # Put starting values in list
        start.val_EXP <- list(startTau, startG, startg0)
        start.val_SIG <- list(startk, startl, startG, startg0)
        
        # Find upper bounds
        upper_EXP <- list(maxTime_EXP, startG+0.02, startg0+0.02)
        upper_SIG <- list(maxTime_SIG, maxTime_SIG, startG+0.02, startg0+0.02)
        
        # Find lower bounds
        lower_EXP <- list(0, startG-0.02, startg0-0.02)
        lower_SIG <- list(0, 0, startG-0.02, startg0-0.02)
        
        #Optimisation
        opt_EXP <- nlminb(start.val_EXP, fit_EXP, lower = lower_EXP, upper = upper_EXP)
        opt_SIG <- nlminb(start.val_SIG, fit_SIG, lower = lower_SIG, upper = upper_SIG)
        
        # Extract result of fitting
        Tau <- as.numeric(opt_EXP$par[1])
        k <- as.numeric(opt_SIG$par[1])
        l <- as.numeric(opt_SIG$par[2])
        G_EXP <- as.numeric(opt_EXP$par[2])
        G_SIG <- as.numeric(opt_SIG$par[3])
        g0_EXP <- as.numeric(opt_EXP$par[3])
        g0_SIG <- as.numeric(opt_SIG$par[4])
        
        # Calculate prediction given the parameters (only for real data points)
        datasec_EXP <- head(secondes_EXP, n=-PltLength) # Delete ending plateau
        datasec_EXP <- tail(secondes_EXP, n=-1)         # Delete 1st point (starting plateau)
        datasec_SIG <- head(secondes_SIG, n=-PltLength) # Delete ending plateau
        datasec_SIG <- tail(secondes_SIG, n=-PltLength) # Delete starting plateau
        
        dataobs_EXP <- head(obs_EXP, n=-PltLength) # Delete ending plateau
        dataobs_EXP <- tail(obs_EXP, n=-1)         # Delete 1st point (starting plateau)
        dataobs_SIG <- head(obs_SIG, n=-PltLength) # Delete ending plateau
        dataobs_SIG <- tail(obs_SIG, n=-PltLength) # Delete starting plateau
        
        pred_EXP <- G_EXP+(g0_EXP-G_EXP)*exp(-(datasec_EXP/Tau)) 
        pred_SIG <- g0_SIG+(G_SIG-g0_SIG)*exp(-exp((l-datasec_SIG)/k))
        
        #plot(pred_EXP~dataobs_EXP)
        #plot(pred_SIG~dataobs_SIG)
        R2_EXP <- as.numeric(summary(lm(pred_EXP~dataobs_EXP))$r.squared)
        R2_SIG <- as.numeric(summary(lm(pred_SIG~dataobs_SIG))$r.squared)
        
        R2_EXP_all <- 1 - fit_EXP(opt_EXP$par)
        R2_SIG_all <- 1 - fit_SIG(opt_SIG$par)
        
        #plot(obs_EXP~secondes_EXP, main="EXP")
        #points(pred_EXP~datasec_EXP, type="l")
        #plot(obs_SIG~secondes_SIG, main="SIG")
        #points(pred_SIG~datasec_SIG, type="l")
        
        #cat(paste(1 - (sum((dataobs_SIG-pred_SIG)^2) / sum((dataobs_SIG-mean(dataobs_SIG))^2))), " | ")
        #cat(paste(1 - fit_SIG(opt_SIG$par)), "\n")
        
        ntempD <- nrow(allresD)
        ntempU <- nrow(allresU)
        
        if (SENS == "DOWN") {
          # If curve goes down
          allresD[ntempD+1,2:12] <- c(as.numeric(Tau), as.numeric(k), as.numeric(l), as.numeric(G_EXP), as.numeric(g0_EXP), as.numeric(G_SIG), as.numeric(g0_SIG), as.numeric(R2_EXP), as.numeric(R2_SIG), as.numeric(R2_EXP_all), as.numeric(R2_SIG_all))
          allresD[ntempD+1,1] <- i
        } else if (SENS == "UP") {
          # If curve goes up
          allresU[ntempU+1,2:12] <- c(as.numeric(Tau), as.numeric(k), as.numeric(l), as.numeric(G_EXP), as.numeric(g0_EXP), as.numeric(G_SIG), as.numeric(g0_SIG), as.numeric(R2_EXP), as.numeric(R2_SIG), as.numeric(R2_EXP_all), as.numeric(R2_SIG_all))
          allresU[ntempU+1,1] <- i
        }
      }
    }
    
    if (nrow(allresD[allresD$id == i,]) > 0){
      
      # Subset the current curve (necessary if more than 1 curve) 
      celuici_D <- allresD[allresD$id == i,]
      
      # Keep only best fits
      resquant_EXPD <- celuici_D[round(as.numeric(celuici_D$R2_EXP_ALL),digits=7)>=round(quantile(celuici_D$R2_EXP_ALL,probs=0.99, na.rm=T),digits=7),]
      resquant_SIGD <- celuici_D[round(celuici_D$R2_SIG_ALL,digits=7)>=round(quantile(celuici_D$R2_SIG_ALL,probs=0.99, na.rm=T),digits=7),]
      
      # Count number of good fit
      Nquant_EXPD <- nrow(resquant_EXPD)
      Nquant_SIGD <- nrow(resquant_SIGD)
      
      # Get which is best R2
      R2quant_EXPD_ALL <- round(quantile(celuici_D$R2_EXP_ALL,probs=0.95, na.rm=T),digits=7)
      R2quant_SIGD_ALL <- round(quantile(celuici_D$R2_SIG_ALL,probs=0.95, na.rm=T),digits=7)
      R2quant_EXPD <- round(quantile(celuici_D$R2_EXP,probs=0.95, na.rm=T),digits=7)
      R2quant_SIGD <- round(quantile(celuici_D$R2_SIG,probs=0.95, na.rm=T),digits=7)
      
      # Get mean parameters for time constants
      mkD <- mean(resquant_SIGD$k, na.rm=T) # At the end : D = down ; U = up
      sdkD <- sd(resquant_SIGD$k, na.rm=T)
      mlD <- mean(resquant_SIGD$l, na.rm=T)
      sdlD <- sd(resquant_SIGD$l, na.rm=T)
      mTauD <- mean(resquant_EXPD$Tau, na.rm=T)
      sdTauD <- sd(resquant_EXPD$Tau, na.rm=T)
      
      # Get mean parameters for steady-states
      mG_EXPD <- mean(resquant_EXPD$G_EXP, na.rm=T)
      sdG_EXPD <- sd(resquant_EXPD$G_EXP, na.rm=T)
      mg0_EXPD <- mean(resquant_EXPD$g0_EXP, na.rm=T)
      sdg0_EXPD <- sd(resquant_EXPD$g0_EXP, na.rm=T)
      mG_SIGD <- mean(resquant_SIGD$G_SIG, na.rm=T)
      sdG_SIGD <- sd(resquant_SIGD$G_SIG, na.rm=T)
      mg0_SIGD <- mean(resquant_SIGD$g0_SIG, na.rm=T)
      sdg0_SIGD <- sd(resquant_SIGD$g0_SIG, na.rm=T)
      
      # Other parameters
      pmax_EXPD <- -(mg0_EXPD - mG_EXPD)/mTauD
      pmax_SIGD <- 1/mkD*(mG_SIGD-mg0_SIGD)/exp(1)
      sdpmax_SIGD <- sqrt((sdkD*exp(-1)*(mG_SIGD-mg0_SIGD)/(mkD*mkD))**2+(sdG_SIGD*exp(-1)/mkD)**2+(sdg0_SIGD*exp(-1)/mkD)**2)
      int_EXPD <- mg0_EXPD
      int_SIGD <- (mG_SIGD-mg0_SIGD)/exp(1) * (1-mlD/mkD) + mg0_SIGD
      mltD <- mlD - LCD
      
      mtempD <- nrow(mresD)
      
      # Fill final data frame with results
      mresD[mtempD+1,2:30] <- c( Nquant_EXPD, Nquant_SIGD, mTauD, sdTauD, mkD, sdkD, mlD, sdlD, LCD, mltD, mg0_EXPD, sdg0_EXPD, mG_EXPD, sdG_EXPD, mg0_SIGD, sdg0_SIGD, mG_SIGD, sdG_SIGD, R2quant_EXPD, R2quant_SIGD, R2quant_EXPD_ALL, R2quant_SIGD_ALL, pmax_EXPD, pmax_SIGD, sdpmax_SIGD, int_EXPD, int_SIGD, A0d, A1d)
      mresD[mtempD+1,1] <- i
      
      pred_final_EXPD <- mG_EXPD+(mg0_EXPD-mG_EXPD)*exp(-(secondes_EXPD / mTauD))
      pred_final_SIGD <- mg0_SIGD+(mG_SIGD-mg0_SIGD)*exp(-exp((mlD-secondes_SIGD)/mkD))
      
      if (plotall == TRUE){
        screen(1, FALSE)
        screen(1, FALSE)
        points(pred_final_SIGD~secondes_SIGD, type="l", col="red", lwd=2)
        abline(int_SIGD,pmax_SIGD,col="blue")
        screen(2, FALSE)
        screen(2, FALSE)
        points(pred_final_EXPD~secondes_EXPD, type="l", col="red", lwd=2)
        abline(int_EXPD,pmax_EXPD,col="blue")
      }
    } 
    
    if (nrow(allresU[allresU$id == i,]) > 0){
      
      # Subset the current curve (necessary if more than 1 curve) 
      celuici_U <- allresU[allresU$id == i,]
      
      # Keep only best fits
      resquant_EXPU <- celuici_U[round(celuici_U$R2_EXP_ALL,digits=7)>=round(quantile(celuici_U$R2_EXP_ALL,probs=0.99, na.rm=T),digits=7),]
      resquant_SIGU <- celuici_U[round(celuici_U$R2_SIG_ALL,digits=7)>=round(quantile(celuici_U$R2_SIG_ALL,probs=0.99, na.rm=T),digits=7),]
      
      # Count number of good fit
      Nquant_EXPU <- nrow(resquant_EXPU)
      Nquant_SIGU <- nrow(resquant_SIGU)
      
      # Get which is best R2
      R2quant_EXPU_ALL <- round(quantile(celuici_U$R2_EXP_ALL,probs=0.95, na.rm=T),digits=7)
      R2quant_SIGU_ALL <- round(quantile(celuici_U$R2_SIG_ALL,probs=0.95, na.rm=T),digits=7)
      R2quant_EXPU <- round(quantile(celuici_U$R2_EXP,probs=0.95, na.rm=T),digits=7)
      R2quant_SIGU <- round(quantile(celuici_U$R2_SIG,probs=0.95, na.rm=T),digits=7)
      
      # Get mean parameters for time constants
      mkU <- mean(resquant_SIGU$k, na.rm=T)
      sdkU <- sd(resquant_SIGU$k, na.rm=T)
      mlU <- mean(resquant_SIGU$l, na.rm=T)
      sdlU <- sd(resquant_SIGU$l, na.rm=T)
      mTauU <- mean(resquant_EXPU$Tau, na.rm=T)
      sdTauU <- sd(resquant_EXPU$Tau, na.rm=T)
      
      # Get mean parameters for steady-states
      mG_EXPU <- mean(resquant_EXPU$G_EXP, na.rm=T)
      sdG_EXPU <- sd(resquant_EXPU$G_EXP, na.rm=T)
      mg0_EXPU <- mean(resquant_EXPU$g0_EXP, na.rm=T)
      sdg0_EXPU <- sd(resquant_EXPU$g0_EXP, na.rm=T)
      mG_SIGU <- mean(resquant_SIGU$G_SIG, na.rm=T)
      sdG_SIGU <- sd(resquant_SIGU$G_SIG, na.rm=T)
      mg0_SIGU <- mean(resquant_SIGU$g0_SIG, na.rm=T)
      sdg0_SIGU <- sd(resquant_SIGU$g0_SIG, na.rm=T)
      
      # Other parameters
      pmax_EXPU <- -(mg0_EXPU - mG_EXPU)/mTauU
      pmax_SIGU <- 1/mkU*(mG_SIGU-mg0_SIGU)/exp(1)
      sdpmax_SIGU <- sqrt((sdkU*exp(-1)*(mG_SIGU-mg0_SIGU)/(mkU*mkU))**2+(sdG_SIGU*exp(-1)/mkU)**2+(sdg0_SIGU*exp(-1)/mkU)**2)
      int_EXPU <- mg0_EXPU
      int_SIGU <- (mG_SIGU-mg0_SIGU)/exp(1) * (1-mlU/mkU) + mg0_SIGU
      mltU <- mlU - LCU
      
      mtempU <- nrow(mresU)
      
      # Fill final data frame with results
      mresU[mtempU+1,2:30] <- c(Nquant_EXPU, Nquant_SIGU, mTauU, sdTauU, mkU, sdkU, mlU, sdlU, LCU, mltU, mg0_EXPU, sdg0_EXPU, mG_EXPU, sdG_EXPU, mg0_SIGU, sdg0_SIGU, mG_SIGU, sdG_SIGU, R2quant_EXPU, R2quant_SIGU, R2quant_EXPU_ALL, R2quant_SIGU_ALL, pmax_EXPU, pmax_SIGU, sdpmax_SIGU, int_EXPU, int_SIGU, A0u, A1u)
      mresU[mtempU+1,1] <- i
      
      pred_final_EXPU <- mG_EXPU+(mg0_EXPU-mG_EXPU)*exp(-(secondes_EXPU / mTauU))
      pred_final_SIGU <- mg0_SIGU+(mG_SIGU-mg0_SIGU)*exp(-exp((mlU-secondes_SIGU)/mkU))
      
      if (plotall == TRUE){
        screen(3, FALSE)
        screen(3, FALSE)
        points(pred_final_SIGU~secondes_SIGU, type="l", col="red", lwd=2)
        abline(int_SIGU,pmax_SIGU,col="blue")
        screen(4, FALSE)
        screen(4, FALSE)
        points(pred_final_EXPU~secondes_EXPU, type="l", col="red", lwd=2)
        abline(int_EXPU,pmax_EXPU,col="blue")
      }
      
    }
    
    if (plotall == TRUE){
      close.screen(all=T)
    }
    
    if (savePlot == TRUE) {
      graphics.off()
      graphics.off()
    }
  }
  
  mresD$Sens <- "DOWN"
  mresU$Sens <- "UP"
  
  allresD$Sens <- "DOWN"
  allresU$Sens <- "UP"
  
  mres <- rbind(mresD, mresU)
  allres <- rbind(allresD, allresU)
  
  RES <- list(mres, allres)
  return(RES)
}
