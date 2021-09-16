source("Function/Func_MxtoLT.R")

CountryName <- "DNK"
# Death rates from the HMD
Mx <- read.table(paste0("Data/", CountryName, "_cMx_1x1.txt"), header = TRUE, fill = TRUE, skip = 2)
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

BC_NA0 <- which(names(wideMx) == "1895")
wideMx <- wideMx[, BC_NA0:ncol(wideMx)]
selwideMx <- wideMx[, colSums(is.na(wideMx)) <= 20]

## survival function of the focal
db_l.f <- matrix(NA, nrow = nrow(selwideMx), ncol = ncol(selwideMx))
for(i in 1:ncol(selwideMx)){
  db_l.f[, i] <- Fun_MxtoLT(x = 1:nrow(selwideMx), Mx = selwideMx[, i], sex = "F", ax = NULL, outcome = "lx")
}

colnames(db_l.f) <- colnames(selwideMx)

BC <- as.numeric(as.character(colnames(db_l.f)))
Age <- 0:110

pdf("Graph/lx_gm_DNK.pdf")
image(x = BC, y = Age, z = t(db_l.f),
      xlab = "Birth cohort of grandmother", ylab = "Age",
      cex = 1.5)
contour(x = BC, y = Age, z = t(db_l.f), add = TRUE, lwd = 3)
dev.off()

#filled.contour(x = BC, y = Age, z = t(db_l.f))

## life expectancy at age 55
eX.f_gm <- c()
for(i in 1:ncol(selwideMx)){
  eX.f_ele <- Fun_MxtoLT(x = 1:nrow(selwideMx), Mx = selwideMx[, i], sex = "F", ax = NULL, outcome = "ex")[56]
  eX.f_gm <- c(eX.f_gm, eX.f_ele)
}

eX.f_gm <- as.vector(as.numeric(eX.f_gm))
pdf("Graph/e55_gm_DNK.pdf")
plot(x = BC, y = eX.f_gm, type = "l", xlab = "Birth cohort of grandmother", ylab = "Life expectancy at age 55")
dev.off()

## qx
db_q.f <- matrix(NA, nrow = nrow(selwideMx), ncol = ncol(selwideMx))
for(i in 1:ncol(selwideMx)){
  db_q.f[, i] <- Fun_MxtoLT(x = 1:nrow(selwideMx), Mx = selwideMx[, i], sex = "F", ax = NULL, outcome = "qx")
}

colnames(db_q.f) <- colnames(selwideMx)

pdf("Graph/qx_gm_DNK.pdf")
image(x = BC, y = Age, z = t(db_q.f),
      xlab = "Birth cohort of grandmother", ylab = "Age",
      cex = 1.5)
contour(x = BC, y = Age, z = t(db_q.f), add = TRUE, lwd = 3)
dev.off()
