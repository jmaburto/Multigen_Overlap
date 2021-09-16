Fun_KGBper <- function(CountryName){
  #### ---- Data ----
  # period death rates from the HMD
  Mx <- read.table(paste0("Data/", CountryName, "_Mx_1x1.txt"), header = TRUE, fill = TRUE, skip = 2)
  # period life table from the HMD
  plt <- read.table(paste0("Data/", CountryName, "_fltper_1x1.txt"), header = TRUE, fill = TRUE, skip = 2)
  # period mean age at birth from the HFD
  pmab <- read.table(paste0("Data/", CountryName, "mabRR.txt"), header = TRUE, fill = TRUE, skip = 2)
  # standard deviation of the period mean age at birth from the HFD
  psdmab <- read.table(paste0("Data/", CountryName, "sdmabRR.txt"), header = TRUE, fill = TRUE, skip = 2)
  # period total fertility rate from the HFD
  ptfr <- read.table(paste0("Data/", CountryName, "tfrRR.txt"), header = TRUE, fill = TRUE, skip = 2)
  
  #### ---- Settings ----
  # settings
  alpha.f <- 15
  beta.f  <- 45
  alpha.m <- 15
  beta.m  <- 60
  maximum <- 100
  
  # survival function for the focal (grandmother)
  plt_wide <- plt %>% 
    as.data.frame() %>%
    select(Year, Age, lx) %>% 
    mutate(Age = ifelse(Age == "110+", "110", Age),
           Age = as.numeric(as.character(Age)),
           lx = lx / 100000) %>%
    spread(key = Year, value = lx) %>% 
    select(-Age)
  
  l.f <- list()
  for(i in 1:ncol(plt_wide)){
    l.f[[i]] <- plt_wide[, i]
  }
  
  # life expectancy
  e0.f_db <- plt %>% 
    as.data.frame() %>%
    filter(Age == "0") %>% 
    select(Year, ex)
  
  e0.f <- as.vector(as.numeric(e0.f_db$ex))
  
  # period for the mortality data
  Year_mort <- as.vector(as.numeric(e0.f_db$Year))
  
  # Mean age at childbirth
  pmab_db <- pmab %>% 
    as.data.frame() %>% 
    mutate(MAB = as.numeric(as.character(MAB)),
           Year = Year) %>% 
    filter(MAB > 0)
  
  Ar <- as.vector(pmab_db$MAB)
  
  # period for the fertility data
  Year_fert <- as.vector(pmab_db$Year)
  
  # The standard deviation of the mean age at childbirth
  psdmab_db <- psdmab %>% 
    as.data.frame() %>% 
    mutate(sdMAB = as.numeric(as.character(sdMAB))) %>% 
    filter(sdMAB > 0)
  
  std <- as.vector(psdmab_db$sdMAB)
  
  # Sex ratio at birth
  SRB <- 1.055
  
  # Period total fertility rate
  ptfr_db <- ptfr %>% 
    as.data.frame() %>% 
    mutate(TFR = as.numeric(as.character(TFR))) %>% 
    filter(TFR > 0)
  
  tfr <- as.vector(ptfr_db$TFR)
  
  # Gross reproduction rate
  GRR <- tfr / (1 + SRB)
  
  # Year that exist both in fertility and mortality data
  Year_matched <- Year_fert[Year_fert %in% Year_mort]
  
  Year_matched_order_fert <- which(Year_fert %in% Year_mort)
  Year_matched_order_mort <- which(Year_mort %in% Year_fert)
  
  # Select
  l.f_gm <- l.f[Year_matched_order_mort]
  e0.f_gm <- e0.f[Year_matched_order_mort]
  
  # Fertility
  m.f_ele <- list(rep(0, maximum))
  m.f <- rep(m.f_ele, length(Ar))
  
  for (index in 1:length(Ar)) {
    for (age in alpha.f:beta.f) {
      # m.f[[index]][age] <- GRR[index] / (beta.f-alpha.f) # fertility rate based on a uniform distribution assumption 
      m.f[[index]][age] <- GRR[index] * dnorm(age, mean = Ar[index], sd = std[index], log = FALSE) / 
        (pnorm(beta.f, mean = Ar[index], sd = std[index]) - pnorm(alpha.f, mean = Ar[index], sd = std[index])) # normal distribution with fixed sd
      # m.f[[index]][age] <- NRR.f[index] * dnorm(age, mean = Ar, sd = sigma.NRR.f[index], log = FALSE)/
      # (pnorm(beta.f, mean = Ar, sd = sigma.NRR.f[index])-pnorm(alpha.f, mean = Ar, sd = sigma.NRR.f[index]))/l.f[[index]][age] #truncated normal distribution
    }
  }
  
  # survival function for daughter and granddaughter
  wideMx <- Mx %>% 
    as.data.frame() %>% 
    select(Year, Age, Female) %>% 
    mutate(Age = ifelse(Age == "110+", "110", Age),
           Age = as.numeric(as.character(Age)),
           Female = as.numeric(as.character(Female)),
           Year = as.numeric(as.character(Year))) %>% 
    filter(Year >= Year_matched[1]) %>% 
    spread(key = Year, value = Female)
  
  db_l.f_full <- matrix(NA, nrow = length(0:maximum), ncol = ncol(wideMx)-1)
  for(i in 2:ncol(wideMx)){
    db_l.f_full[, i-1] <- Fun_MxtoLT(x = 0:maximum, Mx = wideMx[1:(maximum + 1), i], sex = "F", ax = NULL, outcome = "lx")
  }
  
  l.f_full <- list()
  for(i in 1:ncol(db_l.f_full)){
    l.f_full[[i]] <- db_l.f_full[, i]
  }
  
  # Years that exist both in fertility and mortality data
  Year_mort_full <- as.vector(unique(Mx[, "Year"]))
  Year_fert_full <- as.vector(pmab[, "Year"])
  
  Year_mort_full <- Year_mort_full[Year_mort_full >= Year_matched[1]]
  
  names(m.f) <- paste0("Year", Year_fert)
  names(l.f) <- paste0("Year", Year_mort)
  names(l.f_full) <- paste0("Year", Year_mort_full)
  MABp <- Ar %>% 
    as.data.frame() %>% 
    mutate(Ar = round(Ar),
           Year = Year_fert)
  
  #### ---- Prospective Grandparent-grandchild concepts: Average years lived with any grandchild ----
  # Panel A.1 Between Grandmother and daughter's daughter
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
  
  EP1_store_db <- c()
  for(i in 1:length(Year_matched)){
    # period of grandmother (focal), daughter, granddaughter
    Year_gm <- Year_matched[i]
    Year_d <- Year_gm + MABp$Ar[which(MABp$Year == Year_gm)]
    Year_gd <- Year_d + MABp$Ar[which(MABp$Year == Year_d)]
    
    name_Year_gm <- paste0("Year", Year_gm)
    name_Year_d <- paste0("Year", Year_d)
    name_Year_gd <- paste0("Year", Year_gd)
    
    mx_gm = c(0, m.f[[name_Year_gm]])
    e0 = e0.f_gm[i]
    lx_d = c(1, l.f_full[[name_Year_d]])
    lx_gd = c(1, l.f_full[[name_Year_gd]])
    mx_d = c(0, m.f[[name_Year_d]])
    
    EP1 <- Fun_EP1(lx_d = lx_d, lx_gd = lx_gd, mx_gm = mx_gm, mx_d = mx_d, e0 = e0, type = "linear")
    EP1_store <- c(Year_gm, Year_d, Year_gd, EP1, e0, Year_d - Year_gm, Year_gd - Year_d)
    EP1_store_db <- rbind(EP1_store_db, EP1_store)
  }
  colnames(EP1_store_db) <- c("Year_gm", "Year_d", "Year_gd", "Period_EP1", "e0", "MAB_gm", "MAB_d")
  
  #db_EP1_p <- data.frame(Year = seq(1900, 2010, 10),
  #                       Period_EP1 = EP1_SM)
  
  EP1_store_db <- EP1_store_db %>% 
    as.data.frame() %>% 
    filter(Year_gd != 0 & Year_gd <= 2020) %>% 
    #left_join(db_EP1_p, by = c("bc_gm" = "Year")) %>% 
    mutate(diff = e0 - MAB_gm - MAB_d,
           Country = CountryName)
  
  return(EP1_store_db)
}