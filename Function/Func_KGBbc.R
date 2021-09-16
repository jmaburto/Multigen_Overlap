# =================================================================================
# This is the code for the multigenerational overlapping
# project with Christiaan Monden
# 2021/07/15
# updated: 2021/08/30
# Calculate the expected years of multigenerational overlapping from a cohort perspective
# for all countries available at the HFD and HMD
# =================================================================================

Fun_KGBbc <- function(CountryName, agelim){
  #### ---- Data ----
  # Death rates from the HMD
  Mx <- read.table(paste0("Data/", CountryName, "_cMx_1x1.txt"), header = TRUE, fill = TRUE, skip = 2)
  # cohort mean age at birth from the HFD
  cmab <- read.table(paste0("Data/", CountryName, "mabVH.txt"), header = TRUE, fill = TRUE, skip = 2)
  # standard deviation of the cohort mean age at birth from the HFD
  sdmab <- read.table(paste0("Data/", CountryName, "sdmabVH.txt"), header = TRUE, fill = TRUE, skip = 2)
  # cohort total fertility rate from the HFD
  ctfr <- read.table(paste0("Data/", CountryName, "tfrVH.txt"), header = TRUE, fill = TRUE, skip = 2)
  # cohort age-specific fertility rate from the HFD
  casfr <- read.table(paste0("Data/", CountryName, "asfrVH.txt"), header = TRUE, fill = TRUE, skip = 2)
  
  #### ---- Settings ----
  # settings
  alpha.f <- 15
  beta.f  <- 45
  alpha.m <- 15
  beta.m  <- 60
  maximum <- 100
  
  ## Mean age at childbirth
  cmab_db <- cmab %>% 
    as.data.frame() %>% 
    mutate(CMAB = as.numeric(as.character(CMAB)),
           Year = Cohort) %>% 
    filter(CMAB > 0)
  
  Ar <- as.vector(cmab_db$CMAB)
  
  ## The standard deviation of the mean age at childbirth
  sdmab_db <- sdmab %>% 
    as.data.frame() %>% 
    mutate(sdCMAB = as.numeric(as.character(sdCMAB))) %>% 
    filter(sdCMAB > 0)
  
  std <- as.vector(sdmab_db$sdCMAB)
  
  ## Sex ratio at birth
  SRB <- 1.055
  
  ## Cohort total fertility rate
  ctfr_db <- ctfr %>% 
    as.data.frame() %>% 
    mutate(CCF = as.numeric(as.character(CCF))) %>% 
    filter(CCF > 0)
  
  ccf <- as.vector(ctfr_db$CCF)
  
  ## Gross reproduction rate
  GRR <- ccf / (1 + SRB)
  
  ## birth cohort
  BC_fert <- as.vector(as.numeric(cmab_db$Year))
  BC_mort <- as.vector(as.numeric(unique(Mx$Year)))
  BC_matched <- BC_fert[BC_fert %in% BC_mort]
  
  MABbc <- Ar %>% 
    as.data.frame() %>% 
    mutate(Ar = round(Ar),
           bc = BC_fert)
  
  ## death rate for the focal (grandmother)
  wideMx <- Mx %>% 
    as.data.frame() %>% 
    select(Year, Age, Female) %>% 
    mutate(Age = ifelse(Age == "110+", "110", Age),
           Age = as.numeric(as.character(Age)),
           Female = as.numeric(as.character(Female)),
           Year = as.numeric(as.character(Year)),
           Female = ifelse(Year < 1909 & Age >= 105 & is.na(Female), 0, Female)) %>% 
    spread(key = Year, value = Female)
  
  # first birth cohort which does not contain NA
  #BC_NA0 <- which(colSums(is.na(wideMx)) == 0)[2] | which(names(wideMx) == "1900")
  BC_NA0 <- which(names(wideMx) == "1895")
  
  wideMx <- wideMx[, BC_NA0:ncol(wideMx)]
  
  # select Mx having nonNA values at ages 0-90
  AgeLim <- case_when(agelim == 90 ~ 20,
                      agelim == 84 ~ 26)
  selwideMx <- wideMx[, colSums(is.na(wideMx)) <= AgeLim]
  
  # the first column at the NA appears
  col_numNA <- colSums(is.na(selwideMx))
  firstbcNA <- which(col_numNA > 0)[1]
  
  #the last column at the NA appears
  lastbcNA <- ncol(selwideMx) - firstbcNA
  
  # assign the average of 5-year cohorts' mx to NA
  for(i in 0:lastbcNA){
    
    targetbc <- firstbcNA + i
    
    # mean values of 5-year cohorts
    mean5 <- apply(selwideMx[, c(targetbc-5, targetbc-4, targetbc-3, targetbc-2, targetbc-1)], MARGIN = 1, FUN = "mean")
    
    # age at NA
    ageNA <- which(is.na(selwideMx[, targetbc]))
    
    # assign the mean to NA
    selwideMx[ageNA, targetbc] <- mean5[ageNA]
  }
  
  selwideMx <- selwideMx[, names(selwideMx) %in% BC_matched]
  
  if(ncol(selwideMx) == 0){
    
    stop("too small number of columns in Mx")
  }
  
  ## survival function of the focal
  db_l.f <- matrix(NA, nrow = nrow(selwideMx), ncol = ncol(selwideMx))
  for(i in 1:ncol(selwideMx)){
    db_l.f[, i] <- Fun_MxtoLT(x = 1:nrow(selwideMx), Mx = selwideMx[, i], sex = "F", ax = NULL, outcome = "lx")
  }
  
  l.f_gm <- list()
  for(i in 1:ncol(db_l.f)){
    l.f_gm[[i]] <- db_l.f[, i]
  }
  
  ## life expectancy
  e0.f_gm <- c()
  for(i in 1:ncol(selwideMx)){
    e0.f_ele <- Fun_MxtoLT(x = 1:nrow(selwideMx), Mx = selwideMx[, i], sex = "F", ax = NULL, outcome = "ex")[1]
    e0.f_gm <- c(e0.f_gm, e0.f_ele)
  }
  
  e0.f_gm <- as.vector(as.numeric(e0.f_gm))
  
  ## Fertility
  db_casfr <- casfr %>% 
    mutate(Age = case_when(Age == "12-" ~ "12",
                           Age == "55+" ~ "55",
                           T ~ Age),
           Age = as.numeric(as.character(Age)),
           ASFR = as.numeric(as.character(ASFR)),
           ASFR = ifelse((Age <= 15 | Age >= 50) & is.na(ASFR), 0, ASFR),
           Cohort = as.numeric(as.character(Cohort))) %>% 
    filter(Cohort %in% BC_matched,
           Age %in% alpha.f:beta.f) %>% 
    spread(key = Cohort, value = ASFR) %>% 
    select(-Age)
  
  m.f_ele <- list(rep(0, maximum))
  m.f <- rep(m.f_ele, length(Ar))
  
  for (index in 1:length(Ar)) {
    for (age in alpha.f:beta.f) {
      #m.f[[index]][age] <- db_casfr[age-14, index]
      
      # m.f[[index]][age] <- GRR[index] / (beta.f-alpha.f) # fertility rate based on a uniform distribution assumption 
      
      m.f[[index]][age] <- GRR[index] * dnorm(age, mean = Ar[index], sd = std[index], log = FALSE) / 
        (pnorm(beta.f, mean = Ar[index], sd = std[index]) - pnorm(alpha.f, mean = Ar[index], sd = std[index])) # normal distribution with fixed sd
      
      # m.f[[index]][age] <- NRR.f[index] * dnorm(age, mean = Ar, sd = sigma.NRR.f[index], log = FALSE)/
      # (pnorm(beta.f, mean = Ar, sd = sigma.NRR.f[index])-pnorm(alpha.f, mean = Ar, sd = sigma.NRR.f[index]))/l.f[[index]][age] #truncated normal distribution
    }
  }
  names(m.f) <- paste0("bc", BC_matched)
  
  ## survival function for daughter and granddaughter
  wideMx_full <- Mx %>% 
    as.data.frame() %>% 
    select(Year, Age, Female) %>% 
    mutate(Age = ifelse(Age == "110+", "110", Age),
           Age = as.numeric(as.character(Age)),
           Female = as.numeric(as.character(Female)),
           Year = as.numeric(as.character(Year))) %>% 
    filter(Year >= BC_matched[1]) %>%
    spread(key = Year, value = Female)
  
  db_l.f_full <- matrix(NA, nrow = length(0:maximum), ncol = ncol(wideMx_full)-1)
  for(i in 2:ncol(wideMx_full)){
    db_l.f_full[, i-1] <- Fun_MxtoLT(x = 0:maximum, Mx = wideMx_full[1:(maximum + 1), i], sex = "F", ax = NULL, outcome = "lx")
  }
  
  l.f_full <- list()
  for(i in 1:ncol(db_l.f_full)){
    l.f_full[[i]] <- db_l.f_full[, i]
  }
  
  BC_mort_full <- as.numeric(as.character(names(wideMx_full)[-1]))
  names(l.f_full) <- paste0("bc", BC_mort_full)
  
  ## birth cohort of grandmother
  BC_mort_gm <- as.numeric(as.character(names(selwideMx)))
  
  ## Estimate the intrinsic growth rate based on l(x), m(x) and the characteristic function. 
  r.f <- c()
  
  for (index in 1:length(BC_mort_gm)) {
    char.func.f <- function(r){
      
      frac <- 0
      
      for (x in 1:100) { 
        frac <- frac + (l.f_gm[[index]][x] * m.f[[index]][x] * exp(-r * x))
      }
      
      frac <- frac -1
    }
    
    r.f[index] <- uniroot(char.func.f, c(-1, 1))$root
  }
  
  midpoint <- function(vec, type){
    
    n <- length(vec)
    
    if(n == 1){
      
      return(vec)
    } else {
      
      # outcome
      if(type == "linear"){
        
        linear <- (vec[1:(n-1)] + vec[2:n]) * 0.5
        return(linear)
      } else {
        
        #if(unique(vec) == 0){
        #  
        #  exponential <- rep(0, n-1)
        #} else {
        #  
        #  exponential <- vec[1:(n-1)] * (vec[2:n] / vec[1:(n-1)]) ^ 0.5
        #}
        db <- data.frame(a = vec[1:(n-1)],
                         b = vec[2:n])
        db <- db %>% 
          mutate(mid = ifelse(a == 0, 0, a * (b / a) ^ 0.5))
        
        return(db$mid)
      }
    }
  }
  
  ## Prospective Grandparent-grandchild concepts: Average years lived with any grandchild: Between Grandmother and daughter's daughter
  Fun_M2a <- function(a, lx_d, lx_gd, mx_gm, mx_d, type){
    
    store_outside <- c()
    
    for(x in alpha.f:a){
      
      store_inside <- c()
      
      for(y in alpha.f:(a-x)){
        if(a - x - alpha.f >= 0 && a - x - y + 1 > 0){
          
          inside <- lx_d[y+1] * mx_d[y+1] * lx_gd[a-x-y+1]
          
        } else {
          
          inside <- 0
        }
        
        store_inside <- c(store_inside, inside)
      }
      
      # midpoint approximation
      store_inside_mid <- midpoint(store_inside, type = type)
      outcome_inside <- sum(store_inside_mid)
      
      outside <- outcome_inside * mx_gm[x+1]
      store_outside <- c(store_outside, outside)
    }
    
    store_outside_mid <- midpoint(store_outside, type = type)
    outcome <- sum(store_outside_mid, na.rm = T)
    
    return(outcome)
  }
  Fun_EP1 <- function(lx_d, lx_gd, mx_gm, mx_d, e0, type){
    
    store <- c()
    for(a in (2 * alpha.f):e0){
      
      element <- 1 - exp(-Fun_M2a(a = a, lx_d = lx_d, lx_gd = lx_gd, 
                                  mx_gm = mx_gm, mx_d = mx_d, type = type))
      store <- c(store, element)
    }
    
    store_mid <- midpoint(store, type = type)
    EP1 <- sum(store_mid)
    return(EP1)
  }
  
  ## Retrospective Grandparent-grandchild concepts: Average years lived with any grandparent between age 10-19 years: Between Grandmother and daughter's daughter
  Fun_ER1 <- function(r.f, lx_gm, lx_d, mx_gm, mx_d, type){
    
    store_int3 <- c()
    for(a in 10:19){
      
      store_int2 <- c()
      for(y in alpha.f:beta.f){
        
        store_int1 <- c()
        for(x in alpha.f:beta.f){
          #if (a+x+y <= 100) {
          #
          #  ele <- exp(-r.f[i] * (x + y)) * lx_gm[x+y+a+1] * mx_gm[y+1] * lx_d[x+1] * mx_d[x+1]
          #}
          #else {
          #  
          #  ele <- 0
          #}
          
          ele <- exp(-r.f[i] * (x + y)) * lx_gm[x+y+a+1] * mx_gm[y+1] * lx_d[x+1] * mx_d[x+1]
          store_int1 <- c(store_int1, ele)
        }
        
        store_int1_mid <- midpoint(store_int1, type = type)
        outcome_int1 <- sum(store_int1_mid)
        
        store_int2 <- c(store_int2, outcome_int1)
      }
      
      store_int2_mid <- midpoint(store_int2, type = type)
      outcome_int2 <- sum(store_int2_mid)
      
      store_int3 <- c(store_int3, outcome_int2)
    }
    
    store_int3_mid <- midpoint(store_int3, type = type)
    outcome <- sum(store_int3_mid)
    
    return(outcome)
  }
  
  outcome_store_db <- c()
  for(i in 1:length(BC_mort_gm)){
    # birth cohort of grandmother (focal), daughter, granddaughter
    bc_gm <- BC_mort_gm[i]
    bc_d <- bc_gm + MABbc$Ar[which(MABbc$bc == bc_gm)]
    bc_gd <- bc_d + MABbc$Ar[which(MABbc$bc == bc_d)]
    
    name_bc_gm <- paste0("bc", bc_gm)
    name_bc_d <- paste0("bc", bc_d)
    name_bc_gd <- paste0("bc", bc_gd)
    
    mx_gm <- c(0, m.f[[name_bc_gm]])
    e0 <- e0.f_gm[i]
    lx_gm <- l.f_gm[[i]]
    lx_d <- l.f_full[[name_bc_d]]
    lx_gd <- l.f_full[[name_bc_gd]]
    mx_d <- c(0, m.f[[name_bc_d]])
    
    EP1 <- Fun_EP1(lx_d = lx_d, lx_gd = lx_gd, mx_gm = mx_gm, mx_d = mx_d, e0 = e0, type = "linear")
    ER1 <- Fun_ER1(r.f = r.f, lx_gm = lx_gm, lx_d = lx_d, mx_gm = mx_gm, mx_d = mx_d, type = "linear")
    ele_store <- c(bc_gm, bc_d, bc_gd, EP1, ER1, e0, bc_d - bc_gm, bc_gd - bc_d)
    outcome_store_db <- rbind(outcome_store_db, ele_store)
  }
  colnames(outcome_store_db) <- c("bc_gm", "bc_d", "bc_gd", "Cohort_EP1", "Cohort_ER1", "e0", "MAB_gm", "MAB_d")
  
  outcome_store_db <- outcome_store_db %>% 
    as.data.frame() %>% 
    mutate(diff = e0 - MAB_gm - MAB_d,
           Country = CountryName)
  
  return(outcome_store_db)
  
}