# Grandparent-grandchild exposure Song and Mare (2019 Demography)
# Estimate prospective and retrospective gp-gc concepts
# based on parameters from real historical life tables, GRR, r

# In this file, we estimate B1 (grandmother-daughter's daughter; grandmother-daughter's sons;
# grandmother-son's daughter; grandmother-son's sons
# Last modified by Xi Song, 6/24/2017


# B. Prospective Grandparent-grandchild concepts
# B1. Average years lived with any grandchild

# Panel A.1 Between Grandmother and daughter's daughter
b1.gmdgd <- list(b1.gmdgd.level1 = rep(0, 30), b1.gmdgd.level2 = rep(0, 30), b1.gmdgd.level3 = rep(0, 30),
                 b1.gmdgd.level4 = rep(0, 30), b1.gmdgd.level5 = rep(0, 30), b1.gmdgd.level6 = rep(0, 30),
		             b1.gmdgd.level7 = rep(0, 30), b1.gmdgd.level8 = rep(0, 30), b1.gmdgd.level9 = rep(0, 30),
		             b1.gmdgd.level10 = rep(0, 30), b1.gmdgd.level11 = rep(0, 30), b1.gmdgd.level12 = rep(0, 30))
		   
maximum <- 100
b1.gmdgd.a <- list(b1.gmdgd.a.level1 = rep(0, 30), b1.gmdgd.a.level2 = rep(0, 30), b1.gmdgd.a.level3 = rep(0, 30),
                   b1.gmdgd.a.level4 = rep(0, 30), b1.gmdgd.a.level5 = rep(0, 30), b1.gmdgd.a.level6 = rep(0, 30),
		               b1.gmdgd.a.level7 = rep(0, 30), b1.gmdgd.a.level8 = rep(0, 30), b1.gmdgd.a.level9 = rep(0, 30),
		               b1.gmdgd.a.level10 = rep(0, 30), b1.gmdgd.a.level11 = rep(0, 30), b1.gmdgd.a.level12 = rep(0, 30))


	for (index in 1:12) {
	
		for (cum in (2 * alpha.f + 1):maximum) {

			for (a in (2 * alpha.f):cum) {
				temp = 0
				for (x in alpha.f:a) {
					for (y in alpha.f:(a-x)) {
						if (a - x - alpha.f >= 0 && a - x - y > 0) {
							temp = temp + l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][a]/l.f[[index]][x]) * 
									l.f[[index]][y] * m.f[[index]][y] * l.f[[index]][a-x-y] 
						} else {
							temp = temp 
						}
					}
				}
				b1.gmdgd.a[[index]][a] <- min(temp, 1)
				b1.gmdgd[[index]][cum] <- b1.gmdgd[[index]][cum-1] +  b1.gmdgd.a[[index]][a]  
			}  
		}

	}
	
b1.gmdgd #output
b1.gmdgd[[1]][e0.f[1]]; b1.gmdgd[[2]][e0.f[2]]; b1.gmdgd[[3]][e0.f[3]]; 
b1.gmdgd[[4]][e0.f[4]]; b1.gmdgd[[5]][e0.f[5]]; b1.gmdgd[[6]][e0.f[6]];
b1.gmdgd[[7]][e0.f[7]]; b1.gmdgd[[8]][e0.f[8]]; b1.gmdgd[[9]][e0.f[9]]; 
b1.gmdgd[[10]][e0.f[10]]; b1.gmdgd[[11]][e0.f[11]]; b1.gmdgd[[12]][e0.f[12]]  #expected #exposure at birth


# Panel A.2 Between Grandmother and daughter's sons
b1.gmdgs <- list(b1.gmdgs.level1=rep(0,30), b1.gmdgs.level2=rep(0,30), b1.gmdgs.level3=rep(0,30),
                 b1.gmdgs.level4=rep(0,30), b1.gmdgs.level5=rep(0,30), b1.gmdgs.level6=rep(0,30),
				 b1.gmdgs.level7=rep(0,30), b1.gmdgs.level8=rep(0,30), b1.gmdgs.level9=rep(0,30),
                 b1.gmdgs.level10=rep(0,30), b1.gmdgs.level11=rep(0,30), b1.gmdgs.level12=rep(0,30))
		   
maximum <- 100
b1.gmdgs.a <- list(b1.gmdgs.a.level1=rep(0,30), b1.gmdgs.a.level2=rep(0,30), b1.gmdgs.a.level3=rep(0,30),
                   b1.gmdgs.a.level4=rep(0,30), b1.gmdgs.a.level5=rep(0,30), b1.gmdgs.a.level6=rep(0,30),
				   b1.gmdgs.a.level7=rep(0,30), b1.gmdgs.a.level8=rep(0,30), b1.gmdgs.a.level9=rep(0,30),
                   b1.gmdgs.a.level10=rep(0,30), b1.gmdgs.a.level11=rep(0,30), b1.gmdgs.a.level12=rep(0,30))


	for (index in 1 : 12) {
	
		for (cum in (2*alpha.f+1) : maximum) {

			for (a in (2*alpha.f) : cum) {
				temp = 0
				for (x in alpha.f : a) {
					for (y in alpha.f : (a-x)) {
						if (a-x-alpha.f >= 0 && a-x-y > 0) {
							temp = temp + l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][a]/l.f[[index]][x]) * 
									l.f[[index]][y] * m.f[[index]][y] * l.m[[index]][a-x-y] 
						} else {
							temp = temp 
						}
					}
				}
				b1.gmdgs.a[[index]][a] <- min(temp,1)
				b1.gmdgs[[index]][cum] <- b1.gmdgs[[index]][cum-1] +  b1.gmdgs.a[[index]][a]  
			}  
		}

	}
	
b1.gmdgs #output
b1.gmdgs[[1]][e0.f[1]]; b1.gmdgs[[2]][e0.f[2]]; b1.gmdgs[[3]][e0.f[3]]; 
b1.gmdgs[[4]][e0.f[4]]; b1.gmdgs[[5]][e0.f[5]]; b1.gmdgs[[6]][e0.f[6]];
b1.gmdgs[[7]][e0.f[7]]; b1.gmdgs[[8]][e0.f[8]]; b1.gmdgs[[9]][e0.f[9]]; 
b1.gmdgs[[10]][e0.f[10]]; b1.gmdgs[[11]][e0.f[11]]; b1.gmdgs[[12]][e0.f[12]]; #expected #exposure at birth


# Panel A.3 Between Grandmother and son's daughter
b1.gmsgd <- list(b1.gmsgd.level1=rep(0,30), b1.gmsgd.level2=rep(0,30), b1.gmsgd.level3=rep(0,30),
                 b1.gmsgd.level4=rep(0,30), b1.gmsgd.level5=rep(0,30), b1.gmsgd.level6=rep(0,30),
				 b1.gmsgd.level7=rep(0,30), b1.gmsgd.level8=rep(0,30), b1.gmsgd.level9=rep(0,30),
                 b1.gmsgd.level10=rep(0,30), b1.gmsgd.level11=rep(0,30), b1.gmsgd.level12=rep(0,30))
		   
maximum <- 100
b1.gmsgd.a <- list(b1.gmsgd.a.level1=rep(0,30), b1.gmsgd.a.level2=rep(0,30), b1.gmsgd.a.level3=rep(0,30),
                   b1.gmsgd.a.level4=rep(0,30), b1.gmsgd.a.level5=rep(0,30), b1.gmsgd.a.level6=rep(0,30),
				   b1.gmsgd.a.level7=rep(0,30), b1.gmsgd.a.level8=rep(0,30), b1.gmsgd.a.level9=rep(0,30),
                   b1.gmsgd.a.level10=rep(0,30), b1.gmsgd.a.level11=rep(0,30), b1.gmsgd.a.level12=rep(0,30))


	for (index in 1 : 12) {
	
		for (cum in (2*alpha.f+1) : maximum) {

			for (a in (2*alpha.f) : cum) {
				temp = 0
				for (x in alpha.f : a) {
					for (y in alpha.m : (a-x)) {
						if (a-x-alpha.m >= 0 && a-x-y > 0) {
							temp = temp + l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][a]/l.f[[index]][x]) * 
									l.m[[index]][y] * m.m[[index]][y] * l.f[[index]][a-x-y] 
						} else {
							temp = temp 
						}
					}
				}
				b1.gmsgd.a[[index]][a] <- min(temp,1)
				b1.gmsgd[[index]][cum] <- b1.gmsgd[[index]][cum-1] +  b1.gmsgd.a[[index]][a]  
			}  
		}

	}
	
b1.gmsgd #output
b1.gmsgd[[1]][e0.f[1]]; b1.gmsgd[[2]][e0.f[2]]; b1.gmsgd[[3]][e0.f[3]]; 
b1.gmsgd[[4]][e0.f[4]]; b1.gmsgd[[5]][e0.f[5]]; b1.gmsgd[[6]][e0.f[6]];
b1.gmsgd[[7]][e0.f[7]]; b1.gmsgd[[8]][e0.f[8]]; b1.gmsgd[[9]][e0.f[9]]; 
b1.gmsgd[[10]][e0.f[10]]; b1.gmsgd[[11]][e0.f[11]]; b1.gmsgd[[12]][e0.f[12]]; #expected #exposure at birth



# Panel A.4 Between Grandmother and son's sons
b1.gmsgs <- list(b1.gmsgs.level1=rep(0,30), b1.gmsgs.level2=rep(0,30), b1.gmsgs.level3=rep(0,30),
                 b1.gmsgs.level4=rep(0,30), b1.gmsgs.level5=rep(0,30), b1.gmsgs.level6=rep(0,30),
				 b1.gmsgs.level7=rep(0,30), b1.gmsgs.level8=rep(0,30), b1.gmsgs.level9=rep(0,30),
                 b1.gmsgs.level10=rep(0,30), b1.gmsgs.level11=rep(0,30), b1.gmsgs.level12=rep(0,30))
		   
maximum <- 100
b1.gmsgs.a <- list(b1.gmsgs.a.level1=rep(0,30), b1.gmsgs.a.level2=rep(0,30), b1.gmsgs.a.level3=rep(0,30),
                   b1.gmsgs.a.level4=rep(0,30), b1.gmsgs.a.level5=rep(0,30), b1.gmsgs.a.level6=rep(0,30),
				   b1.gmsgs.a.level7=rep(0,30), b1.gmsgs.a.level8=rep(0,30), b1.gmsgs.a.level9=rep(0,30),
                   b1.gmsgs.a.level10=rep(0,30), b1.gmsgs.a.level11=rep(0,30), b1.gmsgs.a.level12=rep(0,30))


	for (index in 1 : 12) {
	
		for (cum in (2*alpha.f+1) : maximum) {

			for (a in (2*alpha.f) : cum) {
				temp = 0
				for (x in alpha.f : a) {
					for (y in alpha.m : (a-x)) {
						if (a-x-alpha.m >= 0 && a-x-y > 0) {
							temp = temp + l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][a]/l.f[[index]][x]) * 
									l.m[[index]][y] * m.m[[index]][y] * l.m[[index]][a-x-y] 
						} else {
							temp = temp 
						}
					}
				}
				b1.gmsgs.a[[index]][a] <- min(temp,1)
				b1.gmsgs[[index]][cum] <- b1.gmsgs[[index]][cum-1] +  b1.gmsgs.a[[index]][a]  
			}  
		}

	}
	
b1.gmsgs #output
b1.gmsgs[[1]][e0.f[1]]; b1.gmsgs[[2]][e0.f[2]]; b1.gmsgs[[3]][e0.f[3]]; 
b1.gmsgs[[4]][e0.f[4]]; b1.gmsgs[[5]][e0.f[5]]; b1.gmsgs[[6]][e0.f[6]]; 
b1.gmsgs[[7]][e0.f[7]]; b1.gmsgs[[8]][e0.f[8]]; b1.gmsgs[[9]][e0.f[9]]; 
b1.gmsgs[[10]][e0.f[10]]; b1.gmsgs[[11]][e0.f[11]]; b1.gmsgs[[12]][e0.f[12]] #expected #exposure at birth



# Panel A Between Grandmother and grandchild (any granddaughter or grandson)

b1.gmgc <- list(b1.gmgc.level1=rep(0,30), b1.gmgc.level2=rep(0,30), b1.gmgc.level3=rep(0,30),
                b1.gmgc.level4=rep(0,30), b1.gmgc.level5=rep(0,30), b1.gmgc.level6=rep(0,30),
				b1.gmgc.level7=rep(0,30), b1.gmgc.level8=rep(0,30), b1.gmgc.level9=rep(0,30),
                b1.gmgc.level10=rep(0,30), b1.gmgc.level11=rep(0,30), b1.gmgc.level12=rep(0,30))
		   
maximum <- 100
b1.gmgc.a <- list(b1.gmgc.a.level1=rep(0,30), b1.gmgc.a.level2=rep(0,30), b1.gmgc.a.level3=rep(0,30),
                  b1.gmgc.a.level4=rep(0,30), b1.gmgc.a.level5=rep(0,30), b1.gmgc.a.level6=rep(0,30),
				  b1.gmgc.a.level7=rep(0,30), b1.gmgc.a.level8=rep(0,30), b1.gmgc.a.level9=rep(0,30),
                  b1.gmgc.a.level10=rep(0,30), b1.gmgc.a.level11=rep(0,30), b1.gmgc.a.level12=rep(0,30))


	for (index in 1 : 12) {
	
		for (cum in (2*alpha.f+1) : maximum) {

			for (a in (2*alpha.f) : cum) {
				temp1 = 0; temp2 = 0; temp3 = 0; temp4 = 0
				for (x in alpha.f : a) {
					for (y in alpha.f : (a-x)) {
						if (a-x-alpha.f >= 0 && a-x-y > 0) {
							temp1 = temp1 + l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][a]/l.f[[index]][x]) * 
									l.f[[index]][y] * m.f[[index]][y] * l.f[[index]][a-x-y] 
						} else {
							temp1 = temp1 
						}
					}
				} # daughter's daughter
				
				for (x in alpha.f : a) {
					for (y in alpha.f : (a-x)) {
						if (a-x-alpha.f >= 0 && a-x-y > 0) {
							temp2 = temp2 + l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][a]/l.f[[index]][x]) * 
									l.f[[index]][y] * m.f[[index]][y] * l.m[[index]][a-x-y] 
						} else {
							temp2 = temp2 
						}
					}
				} # daughter's son
				
				for (x in alpha.f : a) {
					for (y in alpha.m : (a-x)) {
						if (a-x-alpha.m >= 0 && a-x-y > 0) {
							temp3 = temp3 + l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][a]/l.f[[index]][x]) * 
									l.m[[index]][y] * m.m[[index]][y] * l.f[[index]][a-x-y] 
						} else {
							temp3 = temp3 
						}
					}
				} # son's daughter
				
				for (x in alpha.f : a) {
					for (y in alpha.m : (a-x)) {
						if (a-x-alpha.m >= 0 && a-x-y > 0) {
							temp4 = temp4 + l.f[[index]][x] * m.f[[index]][x] * (l.f[[index]][a]/l.f[[index]][x]) * 
									l.m[[index]][y] * m.m[[index]][y] * l.m[[index]][a-x-y] 
						} else {
							temp4 = temp4 
						}
					}
				} # son's son
				
				temp = temp1 + temp2 + temp3 + temp4
				
				b1.gmgc.a[[index]][a] <- min(temp,1)
				b1.gmgc[[index]][cum] <- b1.gmgc[[index]][cum-1] +  b1.gmgc.a[[index]][a]  
			}  
		}

	}
	

b1.gmgc #output
b1.gmgc[[1]][e0.f[1]]; b1.gmgc[[2]][e0.f[2]]; b1.gmgc[[3]][e0.f[3]]; 
b1.gmgc[[4]][e0.f[4]]; b1.gmgc[[5]][e0.f[5]]; b1.gmgc[[6]][e0.f[6]];
b1.gmgc[[7]][e0.f[7]]; b1.gmgc[[8]][e0.f[8]]; b1.gmgc[[9]][e0.f[9]]; 
b1.gmgc[[10]][e0.f[10]]; b1.gmgc[[11]][e0.f[11]]; b1.gmgc[[12]][e0.f[12]] #expected #exposure at birth


