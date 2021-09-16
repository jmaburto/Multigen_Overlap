# Grandparent-grandchild exposure Song and Mare (2019 Demography)
# Estimate prospective and retrospective gp-gc concepts
# based on parameters from real historical life tables, GRR, r

# In this file, we estimate C1 (granddaughter-mother's mother; granddaughter-mother's father;
# granddaughter-father's mother; granddaughter-father's father;
# grandson-mother's mother; grandson-mother's father;
# grandson-father's mother; grandson-father's father)
# Last Modified by Xi Song, 9/24/2014



# C. Retrospective Grandchild-grandparent concepts
# C1. Average years lived with any grandparent

# Panel A.1 Between Granddaughter and maternal grandmother
# Assume a stable population

c1.gdmgm <- list(c1.gdmgm.level1=rep(0,100), c1.gdmgm.level2=rep(0,100), c1.gdmgm.level3=rep(0,100),
                 c1.gdmgm.level4=rep(0,100), c1.gdmgm.level5=rep(0,100), c1.gdmgm.level6=rep(0,100),
				 c1.gdmgm.level7=rep(0,100), c1.gdmgm.level8=rep(0,100), c1.gdmgm.level9=rep(0,100),
                 c1.gdmgm.level10=rep(0,100), c1.gdmgm.level11=rep(0,100), c1.gdmgm.level12=rep(0,100))
c1.gdmgm.a <- list(c1.gdmgm.a.level1=rep(0,100), c1.gdmgm.a.level2=rep(0,100), c1.gdmgm.a.level3=rep(0,100),
                   c1.gdmgm.a.level4=rep(0,100), c1.gdmgm.a.level5=rep(0,100), c1.gdmgm.a.level6=rep(0,100),
				   c1.gdmgm.a.level7=rep(0,100), c1.gdmgm.a.level8=rep(0,100), c1.gdmgm.a.level9=rep(0,100),
                   c1.gdmgm.a.level10=rep(0,100), c1.gdmgm.a.level11=rep(0,100), c1.gdmgm.a.level12=rep(0,100))

	for (index in 1 : 12) {
		for (a in 1 : maximum) {
			for (x in alpha.f : beta.f) {
				for (y in alpha.f : beta.f) {
					if (a+x+y <= 100) {
						c1.gdmgm.a[[index]][a] = c1.gdmgm.a[[index]][a] + exp(-r.f[index]*(x+y)) * l.f[[index]][y] * 
						  m.f[[index]][y] * l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][x+y+a]/l.f[[index]][y])
					} else {
						c1.gdmgm.a[[index]][a] = c1.gdmgm.a[[index]][a] 
					}
				}
			}
			if (a == 1) {
				c1.gdmgm[[index]][a] <- c1.gdmgm.a[[index]][a]
			} else {
				c1.gdmgm[[index]][a] <- c1.gdmgm[[index]][a-1] + c1.gdmgm.a[[index]][a] #one-sex estimation
			}					
		}
	}
	
	
c1.gdmgm.a; c1.gdmgm # output
# c1.gdmgm[[1]][e0.f[1]]; c1.gdmgm[[2]][e0.f[2]]; c1.gdmgm[[3]][e0.f[3]]; 
# c1.gdmgm[[4]][e0.f[4]]; c1.gdmgm[[5]][e0.f[5]]; c1.gdmgm[[6]][e0.f[6]] #expected #exposure at birth

c1.gdmgm[[1]][25]; c1.gdmgm[[2]][25]; c1.gdmgm[[3]][25]; 
c1.gdmgm[[4]][25]; c1.gdmgm[[5]][25]; c1.gdmgm[[6]][25];
c1.gdmgm[[7]][25]; c1.gdmgm[[8]][25]; c1.gdmgm[[9]][25]; 
c1.gdmgm[[10]][25]; c1.gdmgm[[11]][25]; c1.gdmgm[[12]][25]  #expected #exposure at age 25
	

# Panel A.2 Between Granddaughter and maternal grandfather
c1.gdmgf <- list(c1.gdmgf.level1=rep(0,100), c1.gdmgf.level2=rep(0,100), c1.gdmgf.level3=rep(0,100),
                 c1.gdmgf.level4=rep(0,100), c1.gdmgf.level5=rep(0,100), c1.gdmgf.level6=rep(0,100),
				 c1.gdmgf.level7=rep(0,100), c1.gdmgf.level8=rep(0,100), c1.gdmgf.level9=rep(0,100),
                 c1.gdmgf.level10=rep(0,100), c1.gdmgf.level11=rep(0,100), c1.gdmgf.level12=rep(0,100))
c1.gdmgf.a <- list(c1.gdmgf.a.level1=rep(0,100), c1.gdmgf.a.level2=rep(0,100), c1.gdmgf.a.level3=rep(0,100),
                   c1.gdmgf.a.level4=rep(0,100), c1.gdmgf.a.level5=rep(0,100), c1.gdmgf.a.level6=rep(0,100),
				   c1.gdmgf.a.level7=rep(0,100), c1.gdmgf.a.level8=rep(0,100), c1.gdmgf.a.level9=rep(0,100),
                   c1.gdmgf.a.level10=rep(0,100), c1.gdmgf.a.level11=rep(0,100), c1.gdmgf.a.level12=rep(0,100))

	for (index in 1 : 12) {
		for (a in 1 : maximum) {
			for (x in alpha.f : beta.f) {
				for (y in alpha.f : beta.f) {
					if (a+x+y <= 100) {
						c1.gdmgf.a[[index]][a] = c1.gdmgf.a[[index]][a] + exp(-r.f[index]*(x+y)) * l.m[[index]][y] * 
						  m.m[[index]][y] * l.f[[index]][x] * m.f[[index]][x] * (l.m[[index]][x+y+a]/l.m[[index]][y])
					} else {
						c1.gdmgf.a[[index]][a] = c1.gdmgf.a[[index]][a] 
					}
				}
			}
			if (a == 1) {
				c1.gdmgf[[index]][a] <- c1.gdmgf.a[[index]][a]
			} else {
				c1.gdmgf[[index]][a] <- c1.gdmgf[[index]][a-1] + c1.gdmgf.a[[index]][a] #one-sex estimation
			}					
		}
	}
	
	
c1.gdmgf.a; c1.gdmgf # output
# c1.gdmgf[[1]][e0.f[1]]; c1.gdmgf[[2]][e0.f[2]]; c1.gdmgf[[3]][e0.f[3]]; 
# c1.gdmgf[[4]][e0.f[4]]; c1.gdmgf[[5]][e0.f[5]]; c1.gdmgf[[6]][e0.f[6]] #expected #exposure at birth

c1.gdmgf[[1]][25]; c1.gdmgf[[2]][25]; c1.gdmgf[[3]][25]; 
c1.gdmgf[[4]][25]; c1.gdmgf[[5]][25]; c1.gdmgf[[6]][25];
c1.gdmgf[[7]][25]; c1.gdmgf[[8]][25]; c1.gdmgf[[9]][25]; 
c1.gdmgf[[10]][25]; c1.gdmgf[[11]][25]; c1.gdmgf[[12]][25] #expected #exposure at age 25



# Panel A.3 Between Granddaughter and paternal grandmother
c1.gdfgm <- list(c1.gdfgm.level1=rep(0,100), c1.gdfgm.level2=rep(0,100), c1.gdfgm.level3=rep(0,100),
                 c1.gdfgm.level4=rep(0,100), c1.gdfgm.level5=rep(0,100), c1.gdfgm.level6=rep(0,100),
				 c1.gdfgm.level7=rep(0,100), c1.gdfgm.level8=rep(0,100), c1.gdfgm.level9=rep(0,100),
                 c1.gdfgm.level10=rep(0,100), c1.gdfgm.level11=rep(0,100), c1.gdfgm.level12=rep(0,100))
c1.gdfgm.a <- list(c1.gdfgm.a.level1=rep(0,100), c1.gdfgm.a.level2=rep(0,100), c1.gdfgm.a.level3=rep(0,100),
                   c1.gdfgm.a.level4=rep(0,100), c1.gdfgm.a.level5=rep(0,100), c1.gdfgm.a.level6=rep(0,100),
				   c1.gdfgm.a.level7=rep(0,100), c1.gdfgm.a.level8=rep(0,100), c1.gdfgm.a.level9=rep(0,100),
                   c1.gdfgm.a.level10=rep(0,100), c1.gdfgm.a.level11=rep(0,100), c1.gdfgm.a.level12=rep(0,100))

	for (index in 1 : 12) {
		for (a in 1 : maximum) {
			for (x in alpha.f : beta.f) {
				for (y in alpha.f : beta.f) {
					if (a+x+y <= 100) {
						c1.gdfgm.a[[index]][a] = c1.gdfgm.a[[index]][a] + exp(-r.f[index]*(x+y)) * l.f[[index]][y] * 
						  m.f[[index]][y] * l.m[[index]][x] * m.m[[index]][x] * (l.f[[index]][x+y+a]/l.f[[index]][y])
					} else {
						c1.gdfgm.a[[index]][a] = c1.gdfgm.a[[index]][a] 
					}
				}
			}
			if (a == 1) {
				c1.gdfgm[[index]][a] <- c1.gdfgm.a[[index]][a]
			} else {
				c1.gdfgm[[index]][a] <- c1.gdfgm[[index]][a-1] +  c1.gdfgm.a[[index]][a] #one-sex estimation
			}					
		}
	}
	
	
c1.gdfgm.a; c1.gdfgm # output
# c1.gdfgm[[1]][e0.f[1]]; c1.gdfgm[[2]][e0.f[2]]; c1.gdfgm[[3]][e0.f[3]]; 
# c1.gdfgm[[4]][e0.f[4]]; c1.gdfgm[[5]][e0.f[5]]; c1.gdfgm[[6]][e0.f[6]] #expected #exposure at birth

c1.gdfgm[[1]][25]; c1.gdfgm[[2]][25]; c1.gdfgm[[3]][25]; 
c1.gdfgm[[4]][25]; c1.gdfgm[[5]][25]; c1.gdfgm[[6]][25];
c1.gdfgm[[7]][25]; c1.gdfgm[[8]][25]; c1.gdfgm[[9]][25]; 
c1.gdfgm[[10]][25]; c1.gdfgm[[11]][25]; c1.gdfgm[[12]][25] #expected #exposure at age 25


# Panel A.4 Between Granddaughter and paternal grandfather
c1.gdfgf <- list(c1.gdfgf.level1=rep(0,100), c1.gdfgf.level2=rep(0,100), c1.gdfgf.level3=rep(0,100),
                 c1.gdfgf.level4=rep(0,100), c1.gdfgf.level5=rep(0,100), c1.gdfgf.level6=rep(0,100),
				 c1.gdfgf.level7=rep(0,100), c1.gdfgf.level8=rep(0,100), c1.gdfgf.level9=rep(0,100),
                 c1.gdfgf.level10=rep(0,100), c1.gdfgf.level11=rep(0,100), c1.gdfgf.level12=rep(0,100))
c1.gdfgf.a <- list(c1.gdfgf.a.level1=rep(0,100), c1.gdfgf.a.level2=rep(0,100), c1.gdfgf.a.level3=rep(0,100),
                   c1.gdfgf.a.level4=rep(0,100), c1.gdfgf.a.level5=rep(0,100), c1.gdfgf.a.level6=rep(0,100),
				   c1.gdfgf.a.level7=rep(0,100), c1.gdfgf.a.level8=rep(0,100), c1.gdfgf.a.level9=rep(0,100),
                   c1.gdfgf.a.level10=rep(0,100), c1.gdfgf.a.level11=rep(0,100), c1.gdfgf.a.level12=rep(0,100))

	for (index in 1 : 12) {
		for (a in 1 : maximum) {
			for (x in alpha.f : beta.f) {
				for (y in alpha.f : beta.f) {
					if (a+x+y <= 100) {
						c1.gdfgf.a[[index]][a] = c1.gdfgf.a[[index]][a] + exp(-r.f[index]*(x+y)) * l.m[[index]][y] * 
						  m.m[[index]][y] * l.m[[index]][x] * m.m[[index]][x] * (l.m[[index]][x+y+a]/l.m[[index]][y])
					} else {
						c1.gdfgf.a[[index]][a] = c1.gdfgf.a[[index]][a] 
					}
				}
			}
			if (a == 1) {
				c1.gdfgf[[index]][a] <- c1.gdfgf.a[[index]][a]
			} else {
				c1.gdfgf[[index]][a] <- c1.gdfgf[[index]][a-1] +  c1.gdfgf.a[[index]][a] #one-sex estimation
			}					
		}
	}
	
	
c1.gdfgf.a; c1.gdfgf # output
# c1.gdfgf[[1]][e0.f[1]]; c1.gdfgf[[2]][e0.f[2]]; c1.gdfgf[[3]][e0.f[3]]; 
# c1.gdfgf[[4]][e0.f[4]]; c1.gdfgf[[5]][e0.f[5]]; c1.gdfgf[[6]][e0.f[6]] #expected #exposure at birth

c1.gdfgf[[1]][25]; c1.gdfgf[[2]][25]; c1.gdfgf[[3]][25]; 
c1.gdfgf[[4]][25]; c1.gdfgf[[5]][25]; c1.gdfgf[[6]][25];
c1.gdfgf[[7]][25]; c1.gdfgf[[8]][25]; c1.gdfgf[[9]][25]; 
c1.gdfgf[[10]][25]; c1.gdfgf[[11]][25]; c1.gdfgf[[12]][25] #expected #exposure at age 25


# Panel A Between Granddaughter and any grandparent
c1.gdgp <- list(c1.gdgp.level1=rep(0,100), c1.gdgp.level2=rep(0,100), c1.gdgp.level3=rep(0,100),
                c1.gdgp.level4=rep(0,100), c1.gdgp.level5=rep(0,100), c1.gdgp.level6=rep(0,100),
				c1.gdgp.level7=rep(0,100), c1.gdgp.level8=rep(0,100), c1.gdgp.level9=rep(0,100),
                c1.gdgp.level10=rep(0,100), c1.gdgp.level11=rep(0,100), c1.gdgp.level12=rep(0,100))
c1.gdgp.a <- list(c1.gdgp.a.level1=rep(0,100), c1.gdgp.a.level2=rep(0,100), c1.gdgp.a.level3=rep(0,100),
                  c1.gdgp.a.level4=rep(0,100), c1.gdgp.a.level5=rep(0,100), c1.gdgp.a.level6=rep(0,100),
				  c1.gdgp.a.level7=rep(0,100), c1.gdgp.a.level8=rep(0,100), c1.gdgp.a.level9=rep(0,100),
                  c1.gdgp.a.level10=rep(0,100), c1.gdgp.a.level11=rep(0,100), c1.gdgp.a.level12=rep(0,100))
				  
	for (index in 1 : 12) {
		for (a in 1 : maximum) {
			if (a == 1) {
				c1.gdgp[[index]][a] <- 1-(1-c1.gdmgm.a[[index]][a])*(1-c1.gdmgf.a[[index]][a])*(1-c1.gdfgm.a[[index]][a])*(1-c1.gdfgf.a[[index]][a])
			} else {
				c1.gdgp[[index]][a] <- c1.gdgp[[index]][a-1] +  
				1-(1-c1.gdmgm.a[[index]][a])*(1-c1.gdmgf.a[[index]][a])*(1-c1.gdfgm.a[[index]][a])*(1-c1.gdfgf.a[[index]][a]) #one-sex estimation
			}					
		}
	}


	c1.gdgp #output
	
	# c1.gdgp[[1]][e0.f[1]]; c1.gdgp[[2]][e0.f[2]]; c1.gdgp[[3]][e0.f[3]]; 
	# c1.gdgp[[4]][e0.f[4]]; c1.gdgp[[5]][e0.f[5]]; c1.gdgp[[6]][e0.f[6]] #expected #exposure at birth
	
	c1.gdgp[[1]][25]; c1.gdgp[[2]][25]; c1.gdgp[[3]][25]; 
	c1.gdgp[[4]][25]; c1.gdgp[[5]][25]; c1.gdgp[[6]][25];
	c1.gdgp[[7]][25]; c1.gdgp[[8]][25]; c1.gdgp[[9]][25]; 
	c1.gdgp[[10]][25]; c1.gdgp[[11]][25]; c1.gdgp[[12]][25]	#expected #exposure at age 25
	
	
	
##################################################################################################
##################################################################################################
# Panel B.1 Between Grandson and maternal grandmother

# Panel B.2 Between Grandson and maternal grandfather

# Panel B.3 Between Grandson and paternal grandmother

# Panel B.4 Between Grandson and paternal grandfather

# Panel B Between Grandson and any grandparent





