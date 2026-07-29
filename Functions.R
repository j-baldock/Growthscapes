#' ---
#' title: Mechanistic Functions
#' ---
#' 


library(tidyverse)
library(ggpubr)
library(egg)
library(patchwork)

#| cache: false
source("Habitat.R")

#| cache: false
load("data/wt.growth.array.RData")


# Rescale a vector to a specific
fncRescale <- function(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE)) {
    (x - from[1]) / diff(from) * diff(to) + to[1]
}


fncGetBioEParms <- function(spp, pred.en.dens, prey.en.dens, oxy, pff, wt.nearest, 
  startweights = rep(initial.mass, numFish), pvals = rep(0.5, numFish), ration = rep(0.1, numFish)){
  
	N.sites <- nrow(wt.nearest) #here, sites are fish in each reach
  N.steps <- 1 #one time step
  Species <- spp
  SimMethod <- 1 #method that predicts growth
  Pred <- pred.en.dens #predator energy density
  Oxygen <- oxy #oxygen consumed
	PFF <- pff #percent indigestible prey
  stab.factor <- 0.5 #stability factor for other simulation methods
  epsilon <- 0.5 #also for other simulation methods
	endweights <- startweights*5
	TotalConsumption <- rep(100, N.sites)
	pvalues <- t(matrix(round(pvals, 5)))
	sitenames <- t(matrix(wt.nearest[,"pid"]))
	temperature <- t(wt.nearest[,"WT"])
  prey.energy.density <- t(matrix(rep(prey.en.dens, N.sites))) 
  ration <- t(matrix(round(ration, 5)))
  
  return(list(	
    "Species" = Species,
    "SimMethod" = SimMethod,
    "Wstart" = startweights,
    "Endweights" = endweights,
    "TotalConsumption" = TotalConsumption,
    "pp "= pvalues,
    "Temps" = temperature,
    "N.sites" = N.sites,
    "N.steps" = N.steps,
    "sitenames" = sitenames,
    "Pred" = Pred,
    "prey.energy.density" = prey.energy.density,
    "Oxygen" = Oxygen,
    "stab.factor" = stab.factor,
    "PFF" = PFF,
    "epsilon" = epsilon,
    "ration" = ration)
  )
  
}


# Constants, hard-coded - using RATION instead of p-values
fncReadConstants <- function() {
  
  # Get consumption constants
  # Cons <- data.frame(
  #   ConsEQ = 3,
  #   CA = 0.1317,
  #   CB = -0.1396,
  #   CQ = 3.0,
  #   CTO = 15.8,
  #   CTM = 17.5,
  #   CTL = 21.0,
  #   CK1 = 0.06,
  #   CK4 = 0.38)
  # 
  # # Get respiration constants
  # Resp <- data.frame(
  #   RespEQ = 1, 
  #   RA = 0.0009,
  #   RB = -0.1266,
  #   RQ = 0.0833,
  #   RTO = 0,
  #   RTM = 0,
  #   RTL = 0,
  #   RK1 = 1,
  #   RK4 = 0,
  #   ACT = 1,
  #   BACT = 0,
  #   SDA = 0.172)
  # 
  # # Get Excretion / Egestion Constants
  # Excr <- data.frame(
  #   ExcrEQ = 2, 
  #   FA = 0.212,  
  #   FB = -0.222, 
  #   FG = 0.631, 
  #   UA = 0.0314, 
  #   UB = 0.58, 
  #   UG = -0.299)
  
    Cons <- data.frame(
    ConsEQ = 4,
    CA = 0.628,
    CB = -0.3,
    CQ = 5,
    CTO = 20,
    CTM = 20,
    CTL = 24,
    CK1 = 0.33,
    CK4 = 0.2)
  
  # Get respiration constants
  Resp <- data.frame(
    RespEQ = 1, 
    RA = 0.00264,
    RB = -0.217,
    RQ = 0.06818,
    RTO = 0.0234,
    RTM = 0,
    RTL = 25,
    RK1 = 1,
    RK4 = 0.13,
    ACT = 9.7,
    BACT = 0.0405,
    SDA = 0.172)
  
  # Get Excretion / Egestion Constants
  Excr <- data.frame(
    ExcrEQ = 4, 
    FA = 0.212,  
    FB = -0.222, 
    FG = 0.631, 
    UA = 0.0314, 
    UB = 0.58, 
    UG = -0.299)
  
  # Return the Constants
  return(list("Consumption" = Cons, 
              "Respiration" = Resp, 
              "Excretion" = Excr))
}


# Consumption Equation 1
ConsumptionEQ1 <- function(W, TEMP, PP, PREY, CA, CB, CQ) {

	CMAX <- CA * (W ** CB)			#max specific feeding rate (g_prey/g_pred/d)
	CONS <<- (CMAX * PP * exp(CQ * TEMP))	#specific consumption rate (g_prey/g_pred/d) - grams prey consumed per gram of predator mass per day
	CONSj <<- CONS * PREY          		 #specific consumption rate (J/g_pred/d) - Joules consumed for each gram of predator for each day
	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
}

# Consumption Equation 2
ConsumptionEQ2 <- function(W, TEMP, PP, PREY, CA, CB, CTM, CTO, CQ) {

	Y <- log(CQ) * (CTM - CTO + 2)
	Z <- log(CQ) * (CTM - CTO)
	X <- (Z ^ 2 * (1 + (1 + 40 / Y) ^ .5) ^ 2) / 400
	V <- (CTM - TEMP) / (CTM - CTO)
	CMAX <- CA * (W ** CB)		
	CONS <- CMAX * PP * (V ** X) * exp(X * (1 - V))
	CONSj <- CONS*PREY             

	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
}

# Consumption Equation 3: Temperature Dependence for cool-cold water species
ConsumptionEQ3 <- function(W, TEMP, PP, PREY, CA, CB, CK1, CTO, CQ, CK4, CTL, CTM) {
  
	G1 <- (1 / (CTO - CQ)) * (log((0.98 * (1 - CK1)) / (CK1 * 0.02)))
	L1 <- exp(G1 * (TEMP - CQ))
	KA <- (CK1 * L1) / (1 + CK1 * (L1 - 1))
	G2 <- (1 / (CTL - CTM)) * (log((0.98 * (1 - CK4)) / (CK4 * 0.02)))
	L2 <- exp(G2 * (CTL - TEMP))
	KB <- (CK4 * L2) / (1 + CK4 * (L2 - 1))
	CMAX <- CA * (W ** CB)		#max specific feeding rate (g_prey/g_pred/d)
	CONS <- CMAX * PP * KA * KB		#specific consumption rate (g_prey/g_pred/d) - grams prey consumed per gram of predator mass per day
	CONSj <- CONS * PREY             #specific consumption rate (J/g_pred/d) - Joules consumed for each gram of predator for each day

	return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
	}

# Consumption Equation 4: equation 3 but give CONS based on RATION instead of PP
ConsumptionEQ4 <- function(W, TEMP, RATION, PREY, CA, CB, CK1, CTO, CQ, CK4, CTL, CTM) {
  
  G1 <- (1 / (CTO - CQ)) * (log((0.98 * (1 - CK1)) / (CK1 * 0.02)))
  L1 <- exp(G1 * (TEMP - CQ))
  KA <- (CK1 * L1) / (1 + CK1 * (L1 - 1))
  G2 <- (1 / (CTL - CTM)) * (log((0.98 * (1 - CK4)) / (CK4 * 0.02)))
  L2 <- exp(G2 * (CTL - TEMP))
  KB <- (CK4 * L2) / (1 + CK4 * (L2 - 1))
  CMAX <- CA * (W ** CB)		#max specific feeding rate (g_prey/g_pred/d)
  RATION[RATION >= CMAX] <- CMAX[RATION >= CMAX]
  CONS <- RATION * KA * KB  #specific consumption rate (g_prey/g_pred/d) - grams prey consumed per gram of predator mass per day
  CONSj <- CONS * PREY  #specific consumption rate (J/g_pred/d) - Joules consumed for each gram of predator for each day
  
  return(list("CMAX" = CMAX, "CONS" = CONS, "CONSj" = CONSj))
}

# Excretion Equation 1
ExcretionEQ1 <- function(CONS, CONSj, TEMP, PP, FA, UA) {

	EG <- FA * CONS				# egestion (fecal waste) in g_waste/g_pred/d
	U <- UA * (CONS - EG)	 			# excretion (nitrogenous waste) in g_waste/g_pred/d
	EGj <- FA * CONSj				# egestion in J/g/d
	Uj <- UA * (CONSj - EGj)			# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}	

# Excretion Equation 2
ExcretionEQ2 <- function(CONS, CONSj, TEMP, PP, FA, UA, FB, FG, UB, UG) {

	EG <- FA * TEMP ^ FB * exp(FG * PP) * CONS			# egestion (fecal waste) in g_waste/g_pred/d
	U <- UA * TEMP ^ UB * exp(UG * PP) * (CONS - EG)			# excretion (nitrogenous waste) in g_waste/g_pred/d
	EGj <- EG * CONSj / CONS					# egestion in J/g/d
	Uj <- U * CONSj / CONS					# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}

# Excretion Equation 3 (W/ correction for indigestible prey as per Stewart 1983)
ExcretionEQ3 <- function(CONS, CONSj, TEMP, PP, FA, UA, FB, FG, UB, UG, PFF) {

	#Note: In R, "F" means "FALSE", here we use EG as the variable name for egestion instead of F (as in the FishBioE 3.0 manual)
	#Note:  PFF = 0 assumes prey are entirely digestible, making this essentially the same as Equation 2
	PE <- FA * (TEMP ** FB) * exp(FG * PP)
	PF <- ((PE - 0.1) / 0.9) * (1 - PFF) + PFF
	EG <- PF * CONS					# egestion (fecal waste) in g_waste/g_pred/d
	U <- UA * (TEMP ** UB) * (exp(UG * PP)) * (CONS - EG)	# excretion (nitrogenous waste) in g_waste/g_pred/d
	EGj <- PF * CONSj				# egestion in J/g/d
	Uj <- UA * (TEMP ** UB) * (exp(UG * PP)) * (CONSj - EGj)	# excretion in J/g/d

	return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
	}	

# Excretion Equation 3 but using RATION instead of PP
ExcretionEQ4 <- function(CONS, CONSj, TEMP, RATION, CMAX, FA, UA, FB, FG, UB, UG, PFF) {
  
  #Note: In R, "F" means "FALSE", here we use EG as the variable name for egestion instead of F (as in the FishBioE 3.0 manual)
  #Note:  PFF = 0 assumes prey are entirely digestible, making this essentially the same as Equation 2
  RATION[RATION >= CMAX] <- CMAX[RATION >= CMAX]
  PP <- RATION / CMAX #calculating p-value based on ration and CMax inputs
  PE <- FA * (TEMP ** FB) * exp(FG * PP)
  PF <- ((PE - 0.1) / 0.9) * (1 - PFF) + PFF
  EG <- PF * CONS					# egestion (fecal waste) in g_waste/g_pred/d
  U <- UA * (TEMP ** UB) * (exp(UG * PP)) * (CONS - EG)	# excretion (nitrogenous waste) in g_waste/g_pred/d
  EGj <- PF * CONSj				# egestion in J/g/d
  Uj <- UA * (TEMP ** UB) * (exp(UG * PP)) * (CONSj - EGj)	# excretion in J/g/d
  
  return(list("EG" = EG, "EGj" = EGj, "U" = U, "Uj" = Uj))
}	

# Respiration Equation 1
RespirationEQ1 <- function(W, TEMP, CONS, EG, PREY, OXYGEN, RA, RB, ACT, SDA, RQ, RTO, RK1, RK4, RTL, BACT){

	VEL <- (RK1 * W ^ RK4) * (TEMP > RTL) + ACT * W ^ RK4 * exp(BACT * TEMP) * (1 - 1 * (TEMP > RTL))
	ACTIVITY <- exp(RTO * VEL)
	S <- SDA * (CONS - EG)					# proportion of assimilated energy lost to SDA in g/g/d (SDA is unitless)
	Sj <- S * PREY						# proportion of assimilated energy lost to SDA in J/g/d - Joules lost to digestion per gram of predator mass per day
	R <- RA * (W ** RB) * ACTIVITY * exp(RQ * TEMP)       	# energy lost to respiration (metabolism) in g/g/d
	Rj <- R * OXYGEN           				# energy lost to respiration (metabolism) in J/g/d - Joules per gram of predator mass per day
	return(list("R" = R, "Rj" = Rj, "S" = S, "Sj" = Sj))
	}

# Respiration Equation 2 (Temp dependent w/ ACT multiplier)
RespirationEQ2 <- function(W, TEMP, CONS, EG, PREY, OXYGEN, RA, RB, ACT, SDA, RTM, RTO, RQ) {

	V <- (RTM - TEMP) / (RTM - RTO)
	V[V < 0] <- 0.001 #AHF Added to stop errors when Water temps exceed RTM!
	Z <- (log(RQ)) * (RTM - RTO)
	Y <- (log(RQ)) * (RTM - RTO + 2)
	X <- ((Z ** 2) * (1 + (1 + 40 / Y) ** 0.5) ** 2) / 400
	S <- SDA * (CONS - EG)					# proportion of assimilated energy lost to SDA in g/g/d (SDA is unitless)
	Sj <- S * PREY						# proportion of assimilated energy lost to SDA in J/g/d - Joules lost to digestion per gram of predator mass per day
	R <- RA * (W ** RB) * ACT * ((V ** X) * (exp(X * (1 - V) )))     	# energy lost to respiration (metabolism) in g/g/d
	Rj <- R * OXYGEN           				# energy lost to respiration (metabolism) in J/g/d - Joules per gram of predator mass per day
	return(list("R" = R, "Rj" = Rj, "S" = S, "Sj" = Sj))
	}


CalculateGrowth <- function(Constants, Input) {
  
  pred <- Input$Pred
  Oxygen <- Input$Oxygen
  
  # Initialize Fish Weights
  W <- array(rep(0, (Input$N.sites * (Input$N.steps + 1))), c(Input$N.steps + 1, Input$N.sites))
  W[1,] <- as.numeric(Input$Wstart[1:ncol(W)])
  Growth <- array(rep(0, Input$N.sites * Input$N.steps), c(Input$N.steps, Input$N.sites))
  Growth_j <- Consumpt <- Consumpt_j <- Consumpt_cmax <- Prop_cmax <- Excret <- Excret_j <- Egest <- Egest_j <- Respirat <- Respirat_j <- S.resp <- Sj.resp <- Gg_WinBioE <- Gg_ELR <- Growth
  TotalC <- rep(0, Input$N.sites)
  
  ##Start Looping Through Time - for Known Consumption, solving for Weight
  t = 1
  for (t in 1:(Input$N.steps)) {
    
    ### Consumption 
    if(Constants$Consumption$ConsEQ == 1) {
      Cons <- with(Constants$Consumption, ConsumptionEQ1(W[t,], Input$Temps[t,], Input$pp[t,], Input$prey.energy.density[t,], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CQ))
    } else if(Constants$Consumption$ConsEQ == 2){
      Cons <- with(Constants$Consumption, ConsumptionEQ2(W[t,], Input$Temps[t,], Input$pp[t,], Input$prey.energy.density[t,], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CTM, Constants$Consumption$CTO, Constants$Consumption$CQ))
    } else if(Constants$Consumption$ConsEQ == 3){
      Cons <- with(Constants$Consumption, ConsumptionEQ3(W[t,], Input$Temps[t,], Input$pp[t,], Input$prey.energy.density[t,], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CK1, Constants$Consumption$CTO, Constants$Consumption$CQ, Constants$Consumption$CK4, Constants$Consumption$CTL, Constants$Consumption$CTM))
    } else if(Constants$Consumption$ConsEQ == 4){
      Cons <- with(Constants$Consumption, ConsumptionEQ4(W[t,], Input$Temps[t,], Input$ration[t,], Input$prey.energy.density[t,], Constants$Consumption$CA, Constants$Consumption$CB, Constants$Consumption$CK1, Constants$Consumption$CTO, Constants$Consumption$CQ, Constants$Consumption$CK4, Constants$Consumption$CTL, Constants$Consumption$CTM))
    }
    
    TotalC <- TotalC + Cons$CONS * W[t,]
    
    # store daily consumption 
    Consumpt[t,] <- as.numeric(Cons$CONS)
    Consumpt_j[t,] <- as.numeric(Cons$CONSj)
    Consumpt_cmax[t,] <- as.numeric(Cons$CMAX)
    Prop_cmax[t,] <- as.numeric(Cons$CONS) / as.numeric(Cons$CMAX)
    
    
    ### Excretion / Egestion
    if(Constants$Excretion$ExcrEQ == 1) {
      ExcEgest <- with(Constants$Excretion, ExcretionEQ1(Cons$CONS, Cons$CONSj, Input$Temps[t,], Input$pp[t,], Constants$Excretion$FA, Constants$Excretion$UA))
    } else if (Constants$Excretion$ExcrEQ == 2) {
      ExcEgest <- with(Constants$Excretion, ExcretionEQ2(Cons$CONS, Cons$CONSj, Input$Temps[t,], Input$pp[t,], Constants$Excretion$FA, Constants$Excretion$UA, Constants$Excretion$FB, Constants$Excretion$FG, Constants$Excretion$UB, Constants$Excretion$UG ))
    } else if (Constants$Excretion$ExcrEQ == 3) {
      ExcEgest <- with(Constants$Excretion, ExcretionEQ3(Cons$CONS, Cons$CONSj, Input$Temps[t,], Input$pp[t,], Constants$Excretion$FA, Constants$Excretion$UA, Constants$Excretion$FB, Constants$Excretion$FG, Constants$Excretion$UB, Constants$Excretion$UG, Input$PFF) )
    } else if (Constants$Excretion$ExcrEQ == 4) {
      ExcEgest <- with(Constants$Excretion, ExcretionEQ4(Cons$CONS, Cons$CONSj, Input$Temps[t,], Input$ration[t,], Cons$CMAX, Constants$Excretion$FA, Constants$Excretion$UA, Constants$Excretion$FB, Constants$Excretion$FG, Constants$Excretion$UB, Constants$Excretion$UG, Input$PFF) )
    }
    
    # store daily excretion and egestion
    Excret[t,] <- as.numeric(ExcEgest$U)
    Excret_j[t,] <- as.numeric(ExcEgest$Uj)
    Egest[t,] <- as.numeric(ExcEgest$EG)
    Egest_j[t,] <- as.numeric(ExcEgest$EGj)
    
    
    ### Respiration
    if(Constants$Respiration$RespEQ == 1) {
      Resp <- with(Constants$Respiration, RespirationEQ1(W[t,], Input$Temps[t,], Cons$CONS, ExcEgest$EG, Input$prey.energy.density[t,], Input$Oxygen, Constants$Respiration$RA, Constants$Respiration$RB, Constants$Respiration$ACT, Constants$Respiration$SDA, Constants$Respiration$RQ,Constants$Respiration$RTO, Constants$Respiration$RK1, Constants$Respiration$RK4, Constants$Respiration$RTL, Constants$Respiration$BACT))
    } else if (Constants$Respiration$RespEQ == 2) {
      Resp <- with(Constants$Respiration, RespirationEQ2(W[t,], Input$Temps[t,], Cons$CONS, ExcEgest$EG, Input$prey.energy.density[t,], Input$Oxygen, Constants$Respiration$RA, Constants$Respiration$RB, Constants$Respiration$ACT, Constants$Respiration$SDA, Constants$Respiration$RTM, Constants$Respiration$RTO, Constants$Respiration$RQ))
    }
    
    #store daily respiration results
    Respirat[t,] <- as.numeric(Resp$R)
    Respirat_j[t,] <- as.numeric(Resp$Rj)
    S.resp[t,] <- as.numeric(Resp$S)
    Sj.resp[t,] <- as.numeric(Resp$Sj)
    
    
    ### Now calculate Growth
    
    # specific growth in J/g/d - Joules allocated to growth for each gram of predator on each day
    Gj <- Cons$CONSj - Resp$Rj - ExcEgest$EGj - ExcEgest$Uj - Resp$Sj	
    # specific growth in g/g/d - grams allocated to growth for each gram of predator on each day
    G <- Cons$CONS - Resp$R - ExcEgest$EG - ExcEgest$U - Resp$S
    
    ## This original code by Fullerton is a little weird, seemingly inefficient. 
    ## Also the matrix naming scheme is a little wonky:
    ##   - Growth[,]: absolute growth in g/d
    ##   - Growth_j[,]: specific growth in J/g/d
    ##   - Gg_WinBioE[,]: specific growth in g/g/d
    # Growth_raw[t,] <- as.numeric(G) # g/g/d
    
    Growth[t,] <- as.numeric(Gj * W[t,]) / pred # g/d
    Growth_j[t,] <- as.numeric(Gj) # J/g/d
    
    # growth in g/g/d (DailyWeightIncrement divided by fishWeight)
    Gg_WinBioE[t,] <- as.numeric(Growth[t,] / W[t,])		# g/g/d
    
    # I would expect Growth_raw and Gg_WinBioE to be equal, as Gg_WinBioE is simply growth in g/g/d derived from growth in J/g/d. But they are not. Not entirely sure why this is. But like Fullerton's code, the FB4 code also derives growth in g/g/d from J/g/d. FB4 user manual also states that all calculations are done in J/g/d. Perhaps in the "G" equation (line 400), the individual components needs to be adjusted by some factor. Consumption is g_prey/g_pred/day, whereas the other terms are in g_waste/g_pred/day. Perhaps these two units are not directly comparable as I thought...maybe need to be adjust by respective energey densities?
    
    # Calculate absolute weight at time t+1
    W[t + 1,] <- W[t,] + Growth[t,]
    
    
  } # End of cycles through time
  
  # Return Results	
  return(list(
    "TotalC" = TotalC,
    "W" = W, 
    "Growth" = Growth, 
    "Gg_WinBioE" = Gg_WinBioE, 
    "Gg_ELR" = Gg_ELR,
    "Growth_j" = Growth_j,
    # "Growth_raw" = Growth_raw,
    "Consumption" = Consumpt,
    "Consumption_j" = Consumpt_j,
    "Consumption_max" = Consumpt_cmax,
    "Pvalue" = Prop_cmax,
    "Excretion" = Excret,
    "Excretion_j" = Excret_j,
    "Egestion" = Egest,
    "Egestion_j" = Egest_j,
    "Respiration" = Respirat, 
    "Respiration_j" = Respirat_j, 
    "S.resp" = S.resp,
    "Sj.resp" = Sj.resp
  ))
}


BioE <- function (Input, Constants) {

# for simulation method = 1 (we have p-vals, and want to solve for weights
if (Input$SimMethod == 1) {
  
	# Method 1: Calculate Growth from p-values and Temperatures
	Results <- CalculateGrowth(Constants, Input)
	W <- Results$W
	Growth <- Results$Growth
	
} else{

# We don't know p-values, but need to iteratively solve for them
# Method 2 or 3:  Calculate P-values from Total Growth or Consumption
# Need to assume p-values are constant with time
	
	# initialize first guess p-values of .5
	Input$pp <- array(rep(0.5, Input$N.sites * Input$N.step), 
		 c(Input$N.step, Input$N.sites))

	# set error at high value, iterate until it's small
	Error <- rep(99, Input$N.sites)
	iteration <- 0

### Interate until error is less than .1
	while(max(abs(Error)) > Input$epsilon) 
{
	iteration <- iteration + 1
		Results<- CalculateGrowth(Constants, Input)
		W <- Results$W
    TConsumption<- Results$TotalC

 # Find error (depending on which thing on which we're converging), and
 # and come up with new estimate for average p-value
		if (Input$SimMethod == 2)  {
		      Error <- (W[Input$N.step + 1,] - Input$Endweights)
  # Delta is the amount by which we'll change the p-value (prior to
  # scaling by the stability factor
		  	Delta <- (Input$Endweights) / (W[Input$N.step + 1,]) *
      	         Input$pp[1,] - Input$pp[1,]

		Pnew <- as.vector(Input$p[1,] + Input$stab.factor * Delta)
	# Guard against negative p-values
	 for (i in 1:length(Pnew)) {Pnew[i] = max(0, Pnew[i])}
		} else 
		{
	      Error <- (TConsumption - Input$TotalConsumption) / Input$TotalConsumption
	  	Delta <- Input$TotalConsumption/TConsumption * Input$p[1,] - Input$p[1,]
		Pnew <- as.vector(Input$p[1,] + Input$stab.factor * Delta)
	}
	# Update for user, show p-values and error, see if we're converging
		for (i in 1:Input$N.step) {Input$p[i,] <- as.numeric(Pnew)}
			print(paste("Pnew =",Pnew, "  Error =", Error))
	            print(paste("iteration = ", iteration))	
		} # end of while statement 
	} 

Results$pp = Input$pp
# We're done.  Let's return the results!
return(Results)
}


# PIDs: vector of fish IDs, or NA for all fish
# fish: data frame with columns pid, weight, ration, patch ("warm" or "cold")
# t: current day of year
fncGrowthFish <- function(PIDs = NA, fish, T_warm, T_cold) {
  
  # sequences must match the dimensions of the pre-calculated wt.growth array
  wt.seq <- seq(0.1, 25, 0.1)       # water temperature (250 values)
  ra.seq <- seq(0.001, 0.4, 0.001)  # ration (400 values)
  ma.seq <- seq(0.25, 1500, 0.25)         # fish mass (6000 values)
  
  if (any(is.na(PIDs))) {
    td <- fish[, c("pid", "weight", "ration", "patch")]  # all fish
  } else {
    td <- fish[fish$pid %in% PIDs, c("pid", "weight", "ration", "patch")]  # specific fish
  }
  td[, "weight"][td[, "weight"] < 1] <- 1  # minimum lookup weight is 1 g

  if (nrow(td) > 0) {
    n <- nrow(td)
    growth <- watemp <- vector(length = n)
    
    for (x in 1:n) {
      temp     <- if (td[x, "patch"] == "warm") T_warm else T_cold  # temperature at fish's current patch
      wt.idx   <- which.min(abs(wt.seq - temp))              # nearest water temperature index
      ra.idx   <- which.min(abs(ra.seq - td[x, "ration"]))   # nearest ration index
      ma.idx   <- which.min(abs(ma.seq - td[x, "weight"]))   # nearest mass index
      growth[x] <- wt.growth[wt.idx, ra.idx, ma.idx]         # growth rate in g/g/d
      watemp[x] <- wt.seq[wt.idx]
    }
    
    lookup <- data.frame(
      pid      = td[, "pid"],
      weight   = td[, "weight"],
      ration   = td[, "ration"],
      patch    = td[, "patch"],
      WT.actual = ifelse(td[, "patch"] == "warm", T_warm, T_cold),
      WT       = watemp,
      growth   = growth
    )
    return(lookup)
  } else {
    message("No fish IDs provided")
    return(NA)
  }
  
}


# fweight: fish mass in grams
# T_warm, T_cold: current water temperatures (°C) in each patch
# ration_warm, ration_cold: ration sizes available in each patch
fncGrowthPossible <- function(fweight, T_warm, T_cold, ration_warm, ration_cold) {
  
  # sequences must match the dimensions of the pre-calculated wt.growth array
  wt.seq <- seq(0.1, 25, 0.1)       # water temperature (250 values)
  ra.seq <- seq(0.001, 0.4, 0.001)  # ration (400 values)
  ma.seq <- seq(0.25, 1500, 0.25)         # fish mass (6000 values)
  
  fw     <- max(1, min(4500, round(fweight)))  # clamp weight to lookup range
  ma.idx <- which.min(abs(ma.seq - fw))

  patches <- c("warm", "cold")
  temps   <- c(T_warm, T_cold)
  rations <- c(ration_warm, ration_cold)
  growth  <- watemp <- vector(length = 2)
  
  for (x in 1:2) {
    wt.idx    <- which.min(abs(wt.seq - temps[x]))
    ra.idx    <- which.min(abs(ra.seq - rations[x]))
    growth[x] <- wt.growth[wt.idx, ra.idx, ma.idx]
    watemp[x] <- wt.seq[wt.idx]
  }
  
  lookup <- data.frame(
    patch     = patches,
    weight    = rep(fweight, 2),
    ration    = rations,
    WT.actual = c(T_warm, T_cold),
    WT        = watemp,
    growth    = growth
  )
  
  best_patch <- patches[which.max(growth)]
  
  return(list(lookup = lookup, best_patch = best_patch))
  
}


fncTempDepend <- function(temp) {
  # get bioenergetics constants
  constants <- fncReadConstants() 
  
  # Model 3 for cool and coldwater species
  G1 <- (1 / (constants$Consumption$CTO - constants$Consumption$CQ)) * (log((0.98 * (1 - constants$Consumption$CK1)) / (constants$Consumption$CK1 * 0.02)))
  L1 <- exp(G1 * (temp - constants$Consumption$CQ))
  KA <- (constants$Consumption$CK1 * L1) / (1 + constants$Consumption$CK1 * (L1 - 1))
  G2 <- (1 / (constants$Consumption$CTL - constants$Consumption$CTM)) * (log((0.98 * (1 - constants$Consumption$CK4)) / (constants$Consumption$CK4 * 0.02)))
  L2 <- exp(G2 * (constants$Consumption$CTL - temp))
  KB <- (constants$Consumption$CK4 * L2) / (1 + constants$Consumption$CK4 * (L2 - 1))
  return(KA * KB)
}



fncAllomCmax <- function(weight) {
  # get bioenergetics constants
  constants <- fncReadConstants() 
  # calculate allometric cmax = CA * (mass^CB)
  return(constants$Consumption$CA * (weight ^ constants$Consumption$CB))
} 




fncMoveSoftmax <- function(gwarm, gcold, tau = 0.001) {
  exp(gwarm / tau) / (exp(gwarm / tau) + exp(gcold / tau))
}




fncMoveCost <- function(gr = fish$growth) {
  return(abs(gr) * move_cost)
}


fncMoveCost_power <- function(w, c = 0.025, a = 1) {
  c / w^a
}


fncMoveCost_logistic <- function(w, c_max = 0.05, k = 3, w50 = 2) {
  c_max / (1 + exp(k * (w - w50)))
}


fncMoveCost_allometric <- function(w, c = 0.025, b = 0.75) {
  c * w^(b - 1)   # = c * w^(-0.25): slower decay than Option 1
}








fncDensityRationScalar <- function(fish_pop, k_warm = 50, k_cold = 50,
                                    A_warm = 1, A_cold = 1) {
  fish_pop %>%
    group_by(patch) %>%
    mutate(
      A_patch     = if_else(first(patch) == "warm", A_warm, A_cold),
      k           = if_else(first(patch) == "warm", k_warm, k_cold),
      eff_density = (fncEffDensity(weight, beta = 0) - 1) / A_patch,  # beta=0: pure scramble for summary use
      dd_scalar   = k / (k + eff_density)
    ) %>%
    ungroup() %>%
    select(-A_patch, -k, -eff_density)
}



# Dominance-weighted effective density for a vector of fish weights within a patch.
# For focal fish i, effective density = sum_j(alpha_ij) + 1 (for self), where:
#   alpha_ij = 1              if m_j >= m_i  (competitor j is dominant — full competitive cost)
#   alpha_ij = (m_j/m_i)^beta if m_j <  m_i  (focal fish dominates j — reduced cost)
# beta = 0 → pure scramble (all competitors equivalent, recovers standard hyperbolic scalar)
# beta = 1 → linear size-based dominance
# beta → Inf → pure contest (only larger fish impose any cost)
fncEffDensity <- function(weights, beta = 1) {
  n <- length(weights)
  if (n == 1L) return(1)
  ratio_mat      <- outer(weights, weights, FUN = function(mi, mj) ifelse(mj >= mi, 1, (mj / mi)^beta))
  diag(ratio_mat) <- 0   # exclude self from competitor sum
  rowSums(ratio_mat) + 1 # +1 to count focal fish in effective density
}




fncSurviveVari <- function(df, minprob = 0.96, b = 1, b_interact = 0.03, rescale = 0.01){
  # df:          data frame of fish table with only the survivors, e.g., fish[fish$survive == 1, c("weight", "growth")]
  # minprob:     smallest probability any fish can have of dying in any time step
  # b:           controls steepness of the base size-survival curve
  # b_interact:  controls steepness of weight buffering on negative growth; smaller values = growth
  #              effects persist to larger sizes. Default 0.1 means buffering is near-complete ~50g, i.e., no effect of growth >50g.
  # rescale:     controls the relative effect of growth rates on survival. As minprob declines, you need a larger rescale value to generate meaningful variation in survival across the range of observed growth rates
  
  # Weights
  w <- df$weight # weight of fish that are alive at this time step
  
  # Base size-survival curve: 0 for tiny fish, approaches 1 for large fish (controlled by b)
  w_scale <- 1 - 1 / exp(b * w)
  v <- minprob + (1 - minprob) * w_scale
  
  # Growth during this time step that reflects recent conditions (i.e., a hungry/stressed fish may behave in ways that make it more vulnerable to predation, etc.)
  g <- df$growth 
  g <- fncRescale(g, to = c(-rescale, rescale))
  
  # Weight x growth interaction: larger fish are buffered from negative growth penalties
  # (controlled independently from base curve via b_interact).
  # Small fish bear the full cost of negative growth; positive growth benefits are size-independent.
  w_buf    <- 1 - 1 / exp(b_interact * w)
  g_effect <- ifelse(g < 0, g * (1 - w_buf), g)
  
  # Amount of food eaten in previous time step (better survival if bellies full b/c hunkered down somewhere)
  # f <- df$pvals
  # f <- fncRescale(f, to = c(-0.001, 0.001))
  
  # Probability of survival
  prb.srv <- v + g_effect #+ f
  prb.srv[prb.srv > 1] <- 1 # set upper bound at 1
  
  # Sample from binomial distribution with probabilities of prb.srv to determine which fish survive this time step
  survivors <- rbinom(n = nrow(df), size = 1, prob = prb.srv)
  
  return(list(prb.srv, survivors))
  #return(survivors)
}





fncSurviveSize <- function(weight, minprob = 0.96, maxprob = 1, b = 1){
  # weight:      numeric vector of fish weights
  # minprob:     daily survival probability for the smallest fish (asymptote at w = 0)
  # maxprob:     daily survival probability asymptote for very large fish (default 1)
  # b:           controls steepness of the size-survival curve

  # Weights
  w <- weight # weight of fish that are alive at this time step
  
  # Base size-survival curve: 0 for tiny fish, approaches 1 for large fish (controlled by b)
  w_scale <- 1 - 1 / exp(b * w)
  v <- minprob + (maxprob - minprob) * w_scale
  
  # Probability of survival
  prb.srv <- v
  prb.srv[prb.srv > 1] <- 1 # set upper bound at 1
  
  # Sample from binomial distribution with probabilities of prb.srv to determine which fish survive this time step
  survivors <- rbinom(n = length(weight), size = 1, prob = prb.srv)
  
  return(list(prb.srv, survivors))
  #return(survivors)
}





fncSurviveFixed <- function(df, minprob = 0.96, b = 1, b_interact = 0.03, rescale = 0.01, g_range = c(-0.03, 0.03)){
  # df:          data frame of fish table with only the survivors, e.g., fish[fish$survive == 1, c("weight", "growth")]
  # minprob:     smallest probability any fish can have of dying in any time step
  # b:           controls steepness of the base size-survival curve
  # b_interact:  controls steepness of weight buffering on negative growth; smaller values = growth
  #              effects persist to larger sizes. Default 0.03 means buffering is near-complete ~50g, i.e., no effect of growth >50g.
  # rescale:     controls the relative effect of growth rates on survival. As minprob declines, you need a larger rescale value to generate meaningful variation in survival across the range of observed growth rates
  # g_range:     fixed [min, max] instantaneous growth rate (g/g/d) used to rescale growth effects.
  #              Unlike fncSurvive(), this is constant across days, so a fish growing slowly on a
  #              day when all fish grow slowly is not penalized relative to its peers.

  # Weights
  w <- df$weight # weight of fish that are alive at this time step

  # Base size-survival curve: 0 for tiny fish, approaches 1 for large fish (controlled by b)
  w_scale <- 1 - 1 / exp(b * w)
  v <- minprob + (1 - minprob) * w_scale

  # Growth during this time step that reflects recent conditions (i.e., a hungry/stressed fish may behave in ways that make it more vulnerable to predation, etc.)
  g <- df$growth
  # Clamp to fixed range before rescaling so extreme values don't extrapolate beyond [-rescale, rescale]
  g_clamped <- pmax(pmin(g, g_range[2]), g_range[1])
  g <- fncRescale(g_clamped, to = c(-rescale, rescale), from = g_range)

  # Weight x growth interaction: larger fish are buffered from negative growth penalties
  # (controlled independently from base curve via b_interact).
  # Small fish bear the full cost of negative growth; positive growth benefits are size-independent.
  w_buf    <- 1 - 1 / exp(b_interact * w)
  g_effect <- ifelse(g < 0, g * (1 - w_buf), g)

  # Probability of survival
  prb.srv <- v + g_effect
  prb.srv[prb.srv > 1] <- 1 # set upper bound at 1

  # Sample from binomial distribution with probabilities of prb.srv to determine which fish survive this time step
  survivors <- rbinom(n = nrow(df), size = 1, prob = prb.srv)

  return(list(prb.srv, survivors))
}



fncSurviveTemp <- function(temp, T1 = 30, T9 = 25.8) {
  # temp: numeric vector of water temperatures experienced by each fish (e.g., growth_df$WT.actual)
  # T1:   temperature at which daily survival = 0.1 (near-lethal)
  # T9:   temperature at which daily survival = 0.9 (onset of stress)
  # Returns a probability vector (length = length(temp)); multiply with other survival probabilities before sampling.
  T50 <- (T1 + T9) / 2
  k   <- 2 * log(9) / (T1 - T9)
  1 / (1 + exp(k * (temp - T50)))
}



fncSurviveStarve <- function(condition, K9 = 0.55, K1 = 0.45) {
  # condition: numeric vector of relative condition values (current_weight / peak_weight)
  # K9:        condition at which daily survival = 0.9 (onset of starvation effects)
  # K1:        condition at which daily survival = 0.1 (near-lethal)
  # Returns a probability vector; multiply with other survival probabilities before sampling.
  K50 <- (K1 + K9) / 2
  k   <- 2 * log(9) / (K9 - K1)  # note: sign flipped vs. temperature — survival rises with condition
  1 / (1 + exp(-k * (condition - K50)))
}




fncSurviveConsumption <- function(pcmax_dd, age_days,
                                  crit_pcmax_lo    = 0.30,
                                  crit_pcmax_hi    = 0.40,
                                  crit_period_days = 60L) {
  # pcmax_dd:         daily consumption as proportion of Cmax (pcmax_adjusted_dd)
  # age_days:         integer age in days
  # crit_pcmax_lo:    pcmax_dd at which ~1%  period-survival is achieved (near-lethal)
  # crit_pcmax_hi:    pcmax_dd at which ~99% period-survival is achieved (negligible cost)
  # crit_period_days: length of the critical period in days (default 60)
  # Returns a daily survival probability vector of the same length as pcmax_dd.

  out <- rep(1.0, length(pcmax_dd))
  in_crit <- !is.na(age_days) & (age_days <= crit_period_days) & !is.na(pcmax_dd)
  if (!any(in_crit)) return(out)

  # Sigmoid midpoint and steepness anchored to lo/hi thresholds.
  # At lo:  sigmoid ≈ 0.05  (5th percentile → maps to p_min)
  # At hi:  sigmoid ≈ 0.95  (95th percentile → maps to p_max)
  x0 <- (crit_pcmax_lo + crit_pcmax_hi) / 2
  k  <- -log(19) / (crit_pcmax_lo - x0)

  # Daily probability bounds derived from target period-level survival:
  #   p_min: 60-day survival ≈ 1%  → daily = 0.01^(1/T)
  #   p_max: 60-day survival ≈ 99% → daily = 0.99^(1/T)
  p_min <- 0.01 ^ (1 / crit_period_days)
  p_max <- 0.99 ^ (1 / crit_period_days)

  sigmoid      <- 1 / (1 + exp(-k * (pcmax_dd[in_crit] - x0)))
  out[in_crit] <- p_min + (p_max - p_min) * sigmoid
  out
}


fncSurviveAge <- function(age, age_thresh = 5, lambda = 0.05, k = 3, p_min = 0.99) {
  # age:        numeric vector of fish ages in years (pass d / 365 from the IBM loop)
  # age_thresh: age (years) at which senescent mortality begins (default 5)
  # lambda:     scale parameter; larger = faster approach to p_min (default 0.25)
  # k:          Weibull shape exponent; k > 1 gives a slow initial decline that
  #             accelerates with age — the mirror of the concave-down rise in fncSurvive().
  #             k = 1 reduces to a simple exponential (fast drop, decelerating). Default 2.
  # p_min:      asymptotic daily survival floor for very old fish (default 0.99)
  # Returns a probability vector: 1.0 at/below threshold, declining with increasing age above.
  # Default parameters: onset ~age 5, annual survival from aging alone approaches ~0 by age 8.
  excess <- pmax(0, age - age_thresh)
  p_min + (1 - p_min) * exp(-lambda * excess^k)
}



# length in mm
# weight in grams
fncLW <- function(values, input = c("lengths", "weights"), sigma = 0) {
  # values: vector of lengths (mm) or weights (g)
  # input:  type of data in 'values'; output will be the other type
  # sigma:  residual SD on the log10(weight) scale from the LW regression.
  #         Set sigma = 0 (default) for deterministic output.
  #         Biologically plausible values are typically 0.05–0.10.
  #         Noise is applied on the log10 scale (multiplicative on arithmetic scale).
  input <- match.arg(input)
  n <- length(values)
  if (input == "lengths") {
    # Length (mm) -> Weight (g): add noise on log10(W) scale
    noise <- rnorm(n, mean = 0, sd = sigma)
    10^(-5.023 + 3.024 * log10(values) + noise)
  } else {
    # Weight (g) -> Length (mm): noise propagated to log10(L) scale (sigma / b)
    noise <- rnorm(n, mean = 0, sd = sigma / 3.024)
    10^((log10(values) + 5.023) / 3.024 + noise)
  }
}




fncFecund <- function (lengths, sigma = 0) {
  n <- length(lengths)
  noise <- rnorm(n, mean = 0, sd = sigma)
  10^(log10(0.002) + 2.25 * log10(lengths) + noise)
}







fncFecundBromage <- function (weights, sigma = 0, survival = 1) {
  n <- length(weights)
  noise <- rnorm(n, mean = 0, sd = sigma)
  10^(3.461 + 0.504*log10(weights/1000) + noise) * survival
}








fncMaturitySize <- function(weight, K1 = 300, K9 = 500, w_min = 200) {
  # weight:    numeric vector of fish weights
  # K1:        weight at which daily probability of maturity = 0.1 (onset of maturity)
  # K9:        weight at which daily probability of maturity = 0.9 (totally ready)
  # w_min:     minimum weight considered (below this, probability is 0)
  # Returns a probability vector; multiply with other survival probabilities before sampling.
  K50 <- (K1 + K9) / 2
  k   <- 2 * log(9) / (K9 - K1)
  p   <- 1 / (1 + exp(-k * (weight - K50)))
  p[weight < w_min] <- 0
  p
}




fncMaturityCondition <- function(condition, K1 = 0.7, K9 = 0.9, c_min = 0.6) {
  # condition: numeric vector of relative condition values (current_weight / peak_weight)
  # K1:        condition at which daily maturity probability = 0.1 (onset of spawning)
  # K9:        condition at which daily maturity probability = 0.9 (spawning is almost guaranteed)
  # c_min:     minimum condition below which probability is exactly 0 (biologically implausible to spawn)
  # Returns a probability vector; multiply with other survival probabilities before sampling.
  K50 <- (K1 + K9) / 2
  k   <- 2 * log(9) / (K9 - K1)  # note: sign flipped vs. temperature — survival rises with condition
  p   <- 1 / (1 + exp(-k * (condition - K50)))
  p[condition < c_min] <- 0
  p
}




fncMaturityDate <- function(doy, peak_doy = 121, duration = 15, p_max = 1) {
  # doy:       numeric vector of day-of-year values (1–365)
  # peak_doy:  day of year at which daily spawning probability is highest
  # duration:  approximate spawning season width in days; defined as the span
  #            covering ~95% of the probability mass (i.e., ±2 SD from peak),
  #            so sigma = duration / 4
  # p_max:     maximum daily spawning probability, reached at peak_doy
  # Returns a probability vector scaled to [0, p_max]
  sigma <- duration / 4
  p_max * exp(-0.5 * ((doy - peak_doy) / sigma)^2)
}









