Fun_mxbc <- function(CountryName, agelim){
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
 
  return(m.f) 
}

DNK_bc_mx <- Fun_mxbc(CountryName = "DNK", agelim = 90)
DNK_bc_mx <- do.call(cbind.data.frame, DNK_bc_mx)
DNK_bc_mx <- DNK_bc_mx[15:45, ]

BC <- 1901:1970
Age <- 15:45
pdf("Graph/mx_bc_DNK.pdf")
image(x = BC, y = Age, z = t(DNK_bc_mx),
      xlab = "Birth cohort", ylab = "Age",
      cex = 1.5)
contour(x = BC, y = Age, z = t(DNK_bc_mx), add = TRUE, lwd = 3)
dev.off()