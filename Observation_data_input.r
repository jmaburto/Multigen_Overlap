# Grandparent-grandchild exposure Song and Mare (2019 Demography)
# Estimate prospective and retrospective gp-gc concepts
# based on real historical life table, GRR, etc. 

# In this file, we allow GRR, r, e0, and life table to vary across years
# This file provides observational estimates of GP-GC exposures. 
# Last modified by Xi Song, 12/6/2016 corrected an error in 1900 mortality age 45 males
 
# Input parameters

rm(list=ls())


alpha.f <- 15; beta.f <- 45
alpha.m <- 15; beta.m <- 60
maximum <- 100

l.f <- list(l.level1 = rep(0, 100), l.level2 = rep(0, 100), l.level3 = rep(0, 100), l.level4 = rep(0, 100), 
            l.level5 = rep(0, 100), l.level6 = rep(0, 100), l.level7 = rep(0, 100), l.level8 = rep(0, 100),
            l.level9 = rep(0, 100), l.level10 = rep(0, 100), l.level11 = rep(0, 100), l.level12 = rep(0, 100)) # mortality female
l.m <- list(l.level1 = rep(0, 100), l.level2 = rep(0, 100), l.level3 = rep(0, 100), l.level4 = rep(0, 100), 
            l.level5 = rep(0, 100), l.level6 = rep(0, 100), l.level7 = rep(0, 100), l.level8 = rep(0, 100),
			      l.level9 = rep(0, 100), l.level10 = rep(0, 100), l.level11 = rep(0, 100), l.level12 = rep(0, 100)) # male 

m.f <- list(m.level1 = rep(0, 100), m.level2 = rep(0, 100), m.level3 = rep(0, 100), m.level4 = rep(0, 100), 
            m.level5 = rep(0, 100), m.level6 = rep(0, 100), m.level7 = rep(0, 100), m.level8 = rep(0, 100), 
		        m.level9 = rep(0, 100), m.level10 = rep(0, 100), m.level11 = rep(0, 100), m.level12 = rep(0, 100)) # fertility
m.m <- list(m.level1 = rep(0, 100), m.level2 = rep(0, 100), m.level3 = rep(0, 100), m.level4 = rep(0, 100), 
            m.level5 = rep(0, 100), m.level6 = rep(0, 100), m.level7 = rep(0, 100), m.level8 = rep(0, 100), 
		        m.level9 = rep(0, 100), m.level10 = rep(0, 100), m.level11 = rep(0, 100), m.level12 = rep(0, 100))	  

		  

  ###########################################################################
	######################### Male ##########################################
	# 1890 (whites and blacks) l is the L in the life table with radix = 1
	#l.0.m  <- 1;             l.m[[0]][1]  <- 0.83666; l.m[[0]][5]  <- 0.77125;
	#l.m[[0]][10] <- 0.75194; l.m[[0]][15] <- 0.74011; l.m[[0]][20] <- 0.72189; 
	#l.m[[0]][25] <- 0.69519; l.m[[0]][30] <- 0.66657; l.m[[0]][35] <- 0.63551;
	#l.m[[0]][40] <- 0.60169; l.m[[0]][45] <- 0.56537; l.m[[0]][50] <- 0.52444; 
	#l.m[[0]][55] <- 0.47746; l.m[[0]][60] <- 0.41931; l.m[[0]][65] <- 0.35179; 
	#l.m[[0]][70] <- 0.27306; l.m[[0]][75] <- 0.19016; l.m[[0]][80] <- 0.10964; 
	#l.m[[0]][85] <- 0.00000; l.m[[0]][90] <- 0.00000; l.m[[0]][95] <- 0.00000; l.m[[0]][100] <- 0.00000
	
	# 1900 (social security administration)
  l.0.m  <- 1;             l.m[[1]][1]  <- 0.85404; l.m[[1]][5]  <- 0.78591;
  l.m[[1]][10] <- 0.76776; l.m[[1]][15] <- 0.75662; l.m[[1]][20] <- 0.73898; 
  l.m[[1]][25] <- 0.71267; l.m[[1]][30] <- 0.68573; l.m[[1]][35] <- 0.65573;
  l.m[[1]][40] <- 0.62338; l.m[[1]][45] <- 0.58866; l.m[[1]][50] <- 0.54895; 
  l.m[[1]][55] <- 0.50165; l.m[[1]][60] <- 0.44295; l.m[[1]][65] <- 0.37322; 
  l.m[[1]][70] <- 0.29047; l.m[[1]][75] <- 0.19924; l.m[[1]][80] <- 0.11386; 
  l.m[[1]][85] <- 0.04729; l.m[[1]][90] <- 0.01323; l.m[[1]][95] <- 0.00218; l.m[[1]][100] <- 0.00018
	
	# 1900
#	l.0.m  <- 1;             l.m[[1]][1]  <- 0.86426; l.m[[1]][5]  <- 0.80548;
#	l.m[[1]][10] <- 0.78775; l.m[[1]][15] <- 0.77681; l.m[[1]][20] <- 0.75984; 
#	l.m[[1]][25] <- 0.73472; l.m[[1]][30] <- 0.70747; l.m[[1]][35] <- 0.67752;
#	l.m[[1]][40] <- 0.64447; l.m[[1]][45] <- 0.60849; l.m[[1]][50] <- 0.56736 ; 
#	l.m[[1]][55] <- 0.51939; l.m[[1]][60] <- 0.45896; l.m[[1]][65] <- 0.38736; 
#	l.m[[1]][70] <- 0.30217; l.m[[1]][75] <- 0.21076; l.m[[1]][80] <- 0.12084; 
#	l.m[[1]][85] <- 0.05179; l.m[[1]][90] <- 0.01508; l.m[[1]][95] <- 0.00262; l.m[[1]][100] <- 0.00022

	# 1910
	l.0.m  <- 1;             l.m[[2]][1]  <- 0.87505; l.m[[2]][5]  <- 0.82718;
	l.m[[2]][10] <- 0.81249; l.m[[2]][15] <- 0.80261; l.m[[2]][20] <- 0.78792; 
	l.m[[2]][25] <- 0.76675; l.m[[2]][30] <- 0.74378; l.m[[2]][35] <- 0.71614;
	l.m[[2]][40] <- 0.68297; l.m[[2]][45] <- 0.64518; l.m[[2]][50] <- 0.60118; 
	l.m[[2]][55] <- 0.54970; l.m[[2]][60] <- 0.48343; l.m[[2]][65] <- 0.40264; 
	l.m[[2]][70] <- 0.31023; l.m[[2]][75] <- 0.21213; l.m[[2]][80] <- 0.11942; 
	l.m[[2]][85] <- 0.05059; l.m[[2]][90] <- 0.01502; l.m[[2]][95] <- 0.00289; l.m[[2]][100] <- 0.00033
	
	# 1920
	l.0.m  <- 1;             l.m[[3]][1]  <- 0.91745;  l.m[[3]][5]  <- 0.88505;
	l.m[[3]][10] <- 0.87184; l.m[[3]][15] <- 0.86156;  l.m[[3]][20] <- 0.84440; 
	l.m[[3]][25] <- 0.82252; l.m[[3]][30] <- 0.79890;  l.m[[3]][35] <- 0.77514;
	l.m[[3]][40] <- 0.74432; l.m[[3]][45] <- 0.71244;  l.m[[3]][50] <- 0.67553; 
	l.m[[3]][55] <- 0.62965; l.m[[3]][60] <- 0.56917;  l.m[[3]][65] <- 0.49218; 
	l.m[[3]][70] <- 0.39668; l.m[[3]][75] <- 0.28316;  l.m[[3]][80] <- 0.17128; 
	l.m[[3]][85] <- 0.07920; l.m[[3]][90] <- 0.02527;  l.m[[3]][95] <- 0.00556; l.m[[3]][100] <- 0.00062 

	
	# 1930
	l.0.m  <- 1;             l.m[[4]][1]  <- 0.93440;  l.m[[4]][5]  <- 0.91294;
	l.m[[4]][10] <- 0.90346; l.m[[4]][15] <- 0.89561;  l.m[[4]][20] <- 0.88220; 
	l.m[[4]][25] <- 0.86359; l.m[[4]][30] <- 0.84346;  l.m[[4]][35] <- 0.82075;
	l.m[[4]][40] <- 0.79357; l.m[[4]][45] <- 0.75882;  l.m[[4]][50] <- 0.71518; 
	l.m[[4]][55] <- 0.65981; l.m[[4]][60] <- 0.58909;  l.m[[4]][65] <- 0.50154; 
	l.m[[4]][70] <- 0.39516; l.m[[4]][75] <- 0.27718;  l.m[[4]][80] <- 0.16172; 
	l.m[[4]][85] <- 0.07107; l.m[[4]][90] <- 0.02283;  l.m[[4]][95] <- 0.00451; l.m[[4]][100] <- 0.00040

	
	# 1940
	l.0.m  <- 1;             l.m[[5]][1]  <- 0.94762;  l.m[[5]][5]  <- 0.93624;
	l.m[[5]][10] <- 0.93054; l.m[[5]][15] <- 0.92508;  l.m[[5]][20] <- 0.91617; 
	l.m[[5]][25] <- 0.90385; l.m[[5]][30] <- 0.89009;  l.m[[5]][35] <- 0.87371;
	l.m[[5]][40] <- 0.85246; l.m[[5]][45] <- 0.82336;  l.m[[5]][50] <- 0.78254; 
	l.m[[5]][55] <- 0.72627; l.m[[5]][60] <- 0.65142;  l.m[[5]][65] <- 0.55776; 
	l.m[[5]][70] <- 0.44588; l.m[[5]][75] <- 0.31864;  l.m[[5]][80] <- 0.18995; 
	l.m[[5]][85] <- 0.08693; l.m[[5]][90] <- 0.02787;  l.m[[5]][95] <- 0.00586; l.m[[5]][100] <- 0.00078
	
	# 1950
	l.0.m  <- 1;             l.m[[6]][1]  <- 0.96661; l.m[[6]][5]  <- 0.96077;
	l.m[[6]][10] <- 0.95726; l.m[[6]][15] <- 0.95366; l.m[[6]][20] <- 0.94695; 
	l.m[[6]][25] <- 0.93791; l.m[[6]][30] <- 0.92861; l.m[[6]][35] <- 0.91760;
	l.m[[6]][40] <- 0.90207; l.m[[6]][45] <- 0.87819; l.m[[6]][50] <- 0.84158; 
	l.m[[6]][55] <- 0.78781; l.m[[6]][60] <- 0.71246; l.m[[6]][65] <- 0.61566; 
	l.m[[6]][70] <- 0.49950; l.m[[6]][75] <- 0.36756; l.m[[6]][80] <- 0.25237; 
	l.m[[6]][85] <- 0.11750; l.m[[6]][90] <- 0.04197; l.m[[6]][95] <- 0.00955; l.m[[6]][100] <- 0.00121
	
	# 1960
	l.0.m  <- 1;             l.m[[7]][1]  <- 0.97087; l.m[[7]][5]  <- 0.96643;
	l.m[[7]][10] <- 0.96375; l.m[[7]][15] <- 0.96107; l.m[[7]][20] <- 0.95491; 
	l.m[[7]][25] <- 0.94631; l.m[[7]][30] <- 0.93826; l.m[[7]][35] <- 0.92889;
	l.m[[7]][40] <- 0.91572; l.m[[7]][45] <- 0.89492; l.m[[7]][50] <- 0.86199; 
	l.m[[7]][55] <- 0.81039; l.m[[7]][60] <- 0.73887; l.m[[7]][65] <- 0.64177; 
	l.m[[7]][70] <- 0.52244; l.m[[7]][75] <- 0.38950; l.m[[7]][80] <- 0.25300; 
	l.m[[7]][85] <- 0.12845; l.m[[7]][90] <- 0.04609; l.m[[7]][95] <- 0.00970; l.m[[7]][100] <- 0.00117
	
	#1970
	l.0.m  <- 1;             l.m[[8]][1]  <- 0.97755; l.m[[8]][5]  <- 0.97395;
	l.m[[8]][10] <- 0.97151; l.m[[8]][15] <- 0.96904; l.m[[8]][20] <- 0.96126; 
	l.m[[8]][25] <- 0.95040; l.m[[8]][30] <- 0.94072; l.m[[8]][35] <- 0.92997;
	l.m[[8]][40] <- 0.91541; l.m[[8]][45] <- 0.89369; l.m[[8]][50] <- 0.86070; 
	l.m[[8]][55] <- 0.81139; l.m[[8]][60] <- 0.73958; l.m[[8]][65] <- 0.64318; 
	l.m[[8]][70] <- 0.52296; l.m[[8]][75] <- 0.38797; l.m[[8]][80] <- 0.24921; 
	l.m[[8]][85] <- 0.13168; l.m[[8]][90] <- 0.05107; l.m[[8]][95] <- 0.01326; l.m[[8]][100] <- 0.00222
	
	# 1980
	l.0.m  <- 1;             l.m[[9]][1]  <- 0.98607; l.m[[9]][5] <- 0.98333;
	l.m[[9]][10] <- 0.98160; l.m[[9]][15] <- 0.97972; l.m[[9]][20] <- 0.97316; 
	l.m[[9]][25] <- 0.96361; l.m[[9]][30] <- 0.95430; l.m[[9]][35] <- 0.94501;
	l.m[[9]][40] <- 0.93345; l.m[[9]][45] <- 0.91649; l.m[[9]][50] <- 0.89007; 
	l.m[[9]][55] <- 0.84936; l.m[[9]][60] <- 0.79012; l.m[[9]][65] <- 0.70646; 
	l.m[[9]][70] <- 0.59681; l.m[[9]][75] <- 0.46272; l.m[[9]][80] <- 0.31810; 
	l.m[[9]][85] <- 0.18020; l.m[[9]][90] <- 0.07732; l.m[[9]][95] <- 0.02279; l.m[[9]][100] <- 0.00423
	
	# 1990
	l.0.m  <- 1;              l.m[[10]][1]  <- 0.98961; l.m[[10]][5]  <- 0.98754;
	l.m[[10]][10] <- 0.98627; l.m[[10]][15] <- 0.98464; l.m[[10]][20] <- 0.97854; 
	l.m[[10]][25] <- 0.97049; l.m[[10]][30] <- 0.96166; l.m[[10]][35] <- 0.95091;
	l.m[[10]][40] <- 0.93761; l.m[[10]][45] <- 0.92139; l.m[[10]][50] <- 0.89865; 
	l.m[[10]][55] <- 0.86492; l.m[[10]][60] <- 0.81378; l.m[[10]][65] <- 0.73971; 
	l.m[[10]][70] <- 0.64107; l.m[[10]][75] <- 0.51385; l.m[[10]][80] <- 0.36749; 
	l.m[[10]][85] <- 0.21815; l.m[[10]][90] <- 0.09878; l.m[[10]][95] <- 0.02927; l.m[[10]][100] <- 0.00529
	
	# 2000
	l.0.m  <- 1;              l.m[[11]][1]  <- 0.99239; l.m[[11]][5]  <- 0.99095;
	l.m[[11]][10] <- 0.99008; l.m[[11]][15] <- 0.98890; l.m[[11]][20] <- 0.98426; 
	l.m[[11]][25] <- 0.97747; l.m[[11]][30] <- 0.97114; l.m[[11]][35] <- 0.96385;
	l.m[[11]][40] <- 0.95389; l.m[[11]][45] <- 0.93940; l.m[[11]][50] <- 0.91818; 
	l.m[[11]][55] <- 0.88897; l.m[[11]][60] <- 0.84551; l.m[[11]][65] <- 0.78241; 
	l.m[[11]][70] <- 0.69491; l.m[[11]][75] <- 0.57688; l.m[[11]][80] <- 0.42769; 
	l.m[[11]][85] <- 0.26527; l.m[[11]][90] <- 0.12473; l.m[[11]][95] <- 0.03855; l.m[[11]][100] <- 0.00645

	# 2010
	l.0.m  <- 1;              l.m[[12]][1]  <- 0.99333; l.m[[12]][5]  <- 0.99215;
	l.m[[12]][10] <- 0.99151; l.m[[12]][15] <- 0.99071; l.m[[12]][20] <- 0.98727; 
	l.m[[12]][25] <- 0.98105; l.m[[12]][30] <- 0.97441; l.m[[12]][35] <- 0.96724;
	l.m[[12]][40] <- 0.95880; l.m[[12]][45] <- 0.94699; l.m[[12]][50] <- 0.92822; 
	l.m[[12]][55] <- 0.90010; l.m[[12]][60] <- 0.85984; l.m[[12]][65] <- 0.80663; 
	l.m[[12]][70] <- 0.73371; l.m[[12]][75] <- 0.63519; l.m[[12]][80] <- 0.50405; 
	l.m[[12]][85] <- 0.34247; l.m[[12]][90] <- 0.17493; l.m[[12]][95] <- 0.05666; l.m[[12]][100] <- 0.00963

	
	
    ########################### Female #########################################
	# 1890
	#l.0.f  <- 1;             l.f[[0]][1]  <- 0.84235; l.f[[0]][5]  <- 0.77454;
	#l.f[[0]][10] <- 0.75447; l.f[[0]][15] <- 0.74208; l.f[[0]][20] <- 0.72230; 
	#l.f[[0]][25] <- 0.69532; l.f[[0]][30] <- 0.66579; l.f[[0]][35] <- 0.63473;
	#l.f[[0]][40] <- 0.60307; l.f[[0]][45] <- 0.56977; l.f[[0]][50] <- 0.53307; 
	#l.f[[0]][55] <- 0.48954; l.f[[0]][60] <- 0.43553; l.f[[0]][65] <- 0.37266; 
	#l.f[[0]][70] <- 0.29756; l.f[[0]][75] <- 0.21438; l.f[[0]][80] <- 0.13158; 
	#l.f[[0]][85] <- 0.00000; l.f[[0]][90] <- 0.00000; l.f[[0]][95] <- 0.00000; l.f[[0]][100] <- 0.00000
	      
	
	#1900
	l.0.f  <- 1;             l.f[[1]][1]  <- 0.88733; l.f[[1]][5]  <- 0.83119;
	l.f[[1]][10] <- 0.81390; l.f[[1]][15] <- 0.80307; l.f[[1]][20] <- 0.78555; 
	l.f[[1]][25] <- 0.76119; l.f[[1]][30] <- 0.73394; l.f[[1]][35] <- 0.70463;
	l.f[[1]][40] <- 0.67407; l.f[[1]][45] <- 0.64121; l.f[[1]][50] <- 0.60415; 
	l.f[[1]][55] <- 0.55908; l.f[[1]][60] <- 0.50155; l.f[[1]][65] <- 0.43246; 
	l.f[[1]][70] <- 0.34721; l.f[[1]][75] <- 0.24994; l.f[[1]][80] <- 0.15129; 
	l.f[[1]][85] <- 0.07063; l.f[[1]][90] <- 0.02306; l.f[[1]][95] <- 0.00452; l.f[[1]][100] <- 0.00043
	
	# 1910
	l.0.f  <- 1;             l.f[[2]][1]  <- 0.89623; l.f[[2]][5]  <- 0.85117;
	l.f[[2]][10] <- 0.83728; l.f[[2]][15] <- 0.82813; l.f[[2]][20] <- 0.81418; 
	l.f[[2]][25] <- 0.79481; l.f[[2]][30] <- 0.77247; l.f[[2]][35] <- 0.74719;
	l.f[[2]][40] <- 0.71894; l.f[[2]][45] <- 0.68755; l.f[[2]][50] <- 0.65001; 
	l.f[[2]][55] <- 0.60392; l.f[[2]][60] <- 0.54226; l.f[[2]][65] <- 0.46438; 
	l.f[[2]][70] <- 0.36916; l.f[[2]][75] <- 0.26155; l.f[[2]][80] <- 0.15682; 
	l.f[[2]][85] <- 0.07051; l.f[[2]][90] <- 0.02269; l.f[[2]][95] <- 0.00441; l.f[[2]][100] <- 0.00049
	
	# 1920
	l.0.f  <- 1;             l.f[[3]][1]  <- 0.93383;  l.f[[3]][5]  <- 0.90380;
	l.f[[3]][10] <- 0.89186; l.f[[3]][15] <- 0.88247;  l.f[[3]][20] <- 0.86556; 
	l.f[[3]][25] <- 0.84135; l.f[[3]][30] <- 0.81463;  l.f[[3]][35] <- 0.78713;
	l.f[[3]][40] <- 0.75907; l.f[[3]][45] <- 0.72954;  l.f[[3]][50] <- 0.69452; 
	l.f[[3]][55] <- 0.65099; l.f[[3]][60] <- 0.59438;  l.f[[3]][65] <- 0.52126; 
	l.f[[3]][70] <- 0.42741; l.f[[3]][75] <- 0.31344;  l.f[[3]][80] <- 0.19613; 
	l.f[[3]][85] <- 0.09515; l.f[[3]][90] <- 0.03314;  l.f[[3]][95] <- 0.00728; l.f[[3]][100] <- 0.00072 
	
	# 1930
	l.0.f  <- 1;             l.f[[4]][1]  <- 0.94728;  l.f[[4]][5]  <- 0.92789;
	l.f[[4]][10] <- 0.92008; l.f[[4]][15] <- 0.91364;  l.f[[4]][20] <- 0.90116; 
	l.f[[4]][25] <- 0.88328; l.f[[4]][30] <- 0.86398;  l.f[[4]][35] <- 0.84304;
	l.f[[4]][40] <- 0.81927; l.f[[4]][45] <- 0.79041;  l.f[[4]][50] <- 0.75456; 
	l.f[[4]][55] <- 0.70832; l.f[[4]][60] <- 0.64795;  l.f[[4]][65] <- 0.56924; 
	l.f[[4]][70] <- 0.46774; l.f[[4]][75] <- 0.34600;  l.f[[4]][80] <- 0.21578; 
	l.f[[4]][85] <- 0.10322; l.f[[4]][90] <- 0.03656;  l.f[[4]][95] <- 0.00807; l.f[[4]][100] <- 0.00082
	
	
	# 1940
	l.0.f  <- 1;             l.f[[5]][1]  <- 0.95848;  l.f[[5]][5]  <- 0.94858;
	l.f[[5]][10] <- 0.94402; l.f[[5]][15] <- 0.94000;  l.f[[5]][20] <- 0.93293; 
	l.f[[5]][25] <- 0.92322; l.f[[5]][30] <- 0.91182;  l.f[[5]][35] <- 0.89810;
	l.f[[5]][40] <- 0.88092; l.f[[5]][45] <- 0.85856;  l.f[[5]][50] <- 0.82828; 
	l.f[[5]][55] <- 0.78708; l.f[[5]][60] <- 0.73093;  l.f[[5]][65] <- 0.65523; 
	l.f[[5]][70] <- 0.55449; l.f[[5]][75] <- 0.42425;  l.f[[5]][80] <- 0.27524; 
	l.f[[5]][85] <- 0.13972; l.f[[5]][90] <- 0.05044;  l.f[[5]][95] <- 0.01195; l.f[[5]][100] <- 0.00179
	
	# 1950
	l.0.f  <- 1;             l.f[[6]][1]  <- 0.97406; l.f[[6]][5]  <- 0.96908;
	l.f[[6]][10] <- 0.96652; l.f[[6]][15] <- 0.96431; l.f[[6]][20] <- 0.96066; 
	l.f[[6]][25] <- 0.95583; l.f[[6]][30] <- 0.94993; l.f[[6]][35] <- 0.94206;
	l.f[[6]][40] <- 0.93101; l.f[[6]][45] <- 0.91469; l.f[[6]][50] <- 0.89075; 
	l.f[[6]][55] <- 0.85694; l.f[[6]][60] <- 0.80890; l.f[[6]][65] <- 0.74119; 
	l.f[[6]][70] <- 0.64873; l.f[[6]][75] <- 0.52111; l.f[[6]][80] <- 0.36486; 
	l.f[[6]][85] <- 0.20668; l.f[[6]][90] <- 0.08548; l.f[[6]][95] <- 0.02207; l.f[[6]][100] <- 0.00298
	
	# 1960
	l.0.f  <- 1;             l.f[[7]][1]  <- 0.97744; l.f[[7]][5]  <- 0.97371;
	l.f[[7]][10] <- 0.97173; l.f[[7]][15] <- 0.97016; l.f[[7]][20] <- 0.96756; 
	l.f[[7]][25] <- 0.96418; l.f[[7]][30] <- 0.95996; l.f[[7]][35] <- 0.95409;
	l.f[[7]][40] <- 0.94560; l.f[[7]][45] <- 0.93265; l.f[[7]][50] <- 0.91327; 
	l.f[[7]][55] <- 0.88451; l.f[[7]][60] <- 0.84430; l.f[[7]][65] <- 0.78462; 
	l.f[[7]][70] <- 0.70100; l.f[[7]][75] <- 0.58394; l.f[[7]][80] <- 0.43063; 
	l.f[[7]][85] <- 0.25269; l.f[[7]][90] <- 0.10056; l.f[[7]][95] <- 0.02193; l.f[[7]][100] <- 0.00264

	#1970
	l.0.f  <- 1;             l.f[[8]][1]  <- 0.98254; l.f[[8]][5]  <- 0.97955;
	l.f[[8]][10] <- 0.97784; l.f[[8]][15] <- 0.97636; l.f[[8]][20] <- 0.97331; 
	l.f[[8]][25] <- 0.96966; l.f[[8]][30] <- 0.96544; l.f[[8]][35] <- 0.95966;
	l.f[[8]][40] <- 0.95097; l.f[[8]][45] <- 0.93793; l.f[[8]][50] <- 0.91852; 
	l.f[[8]][55] <- 0.89066; l.f[[8]][60] <- 0.85139; l.f[[8]][65] <- 0.79698; 
	l.f[[8]][70] <- 0.71955; l.f[[8]][75] <- 0.61107; l.f[[8]][80] <- 0.46445; 
	l.f[[8]][85] <- 0.29538; l.f[[8]][90] <- 0.14160; l.f[[8]][95] <- 0.04565; l.f[[8]][100] <- 0.00954
	
	# 1980
	l.0.f  <- 1;             l.f[[9]][1]  <- 0.98880; l.f[[9]][5]  <- 0.98666;
	l.f[[9]][10] <- 0.98544; l.f[[9]][15] <- 0.98432; l.f[[9]][20] <- 0.98184; 
	l.f[[9]][25] <- 0.97883; l.f[[9]][30] <- 0.97551; l.f[[9]][35] <- 0.97140;
	l.f[[9]][40] <- 0.96531; l.f[[9]][45] <- 0.95570; l.f[[9]][50] <- 0.94060; 
	l.f[[9]][55] <- 0.91760; l.f[[9]][60] <- 0.88414; l.f[[9]][65] <- 0.83520; 
	l.f[[9]][70] <- 0.76720; l.f[[9]][75] <- 0.67186; l.f[[9]][80] <- 0.54372; 
	l.f[[9]][85] <- 0.37772; l.f[[9]][90] <- 0.20578; l.f[[9]][95] <- 0.07862; l.f[[9]][100] <- 0.01927

	
	# 1990
	l.0.f  <- 1;              l.f[[10]][1]  <- 0.99172; l.f[[10]][5]  <- 0.99006;
	l.f[[10]][10] <- 0.98911; l.f[[10]][15] <- 0.98814; l.f[[10]][20] <- 0.98597; 
	l.f[[10]][25] <- 0.98325; l.f[[10]][30] <- 0.98013; l.f[[10]][35] <- 0.97596;
	l.f[[10]][40] <- 0.97033; l.f[[10]][45] <- 0.96222; l.f[[10]][50] <- 0.94932; 
	l.f[[10]][55] <- 0.92881; l.f[[10]][60] <- 0.89742; l.f[[10]][65] <- 0.85075; 
	l.f[[10]][70] <- 0.78522; l.f[[10]][75] <- 0.69287; l.f[[10]][80] <- 0.56986; 
	l.f[[10]][85] <- 0.41115; l.f[[10]][90] <- 0.23666; l.f[[10]][95] <- 0.09346; l.f[[10]][100] <- 0.02251

	                                        
	# 2000                                  
	l.0.f  <- 1;              l.f[[11]][1]  <- 0.99375; l.f[[11]][5]  <- 0.99261;
	l.f[[11]][10] <- 0.99190; l.f[[11]][15] <- 0.99111; l.f[[11]][20] <- 0.98915; 
	l.f[[11]][25] <- 0.98682; l.f[[11]][30] <- 0.98418; l.f[[11]][35] <- 0.98052;
	l.f[[11]][40] <- 0.97493; l.f[[11]][45] <- 0.96648; l.f[[11]][50] <- 0.95425; 
	l.f[[11]][55] <- 0.93609; l.f[[11]][60] <- 0.90767; l.f[[11]][65] <- 0.86433; 
	l.f[[11]][70] <- 0.80219; l.f[[11]][75] <- 0.71311; l.f[[11]][80] <- 0.58455; 
	l.f[[11]][85] <- 0.41830; l.f[[11]][90] <- 0.23936; l.f[[11]][95] <- 0.09560; l.f[[11]][100] <- 0.02183

	                                        
	# 2010                                  
	l.0.f  <- 1;              l.f[[12]][1]  <- 0.99445; l.f[[12]][5]  <- 0.99351;
	l.f[[12]][10] <- 0.99301; l.f[[12]][15] <- 0.99241; l.f[[12]][20] <- 0.99102; 
	l.f[[12]][25] <- 0.98880; l.f[[12]][30] <- 0.98604; l.f[[12]][35] <- 0.98247;
	l.f[[12]][40] <- 0.97745; l.f[[12]][45] <- 0.96996; l.f[[12]][50] <- 0.95798; 
	l.f[[12]][55] <- 0.94018; l.f[[12]][60] <- 0.91575; l.f[[12]][65] <- 0.88040; 
	l.f[[12]][70] <- 0.82760; l.f[[12]][75] <- 0.75037; l.f[[12]][80] <- 0.63820; 
	l.f[[12]][85] <- 0.48344; l.f[[12]][90] <- 0.29178; l.f[[12]][95] <- 0.12005; l.f[[12]][100] <- 0.02758

	                                             
	
	#linear interpolation
	for (index in 1 : 12) {
		for (x in 2 : 4) {
			l.f[[index]][x] <- l.f[[index]][1] + (l.f[[index]][5] - l.f[[index]][1])/4 * (x-1)
			l.m[[index]][x] <- l.m[[index]][1] + (l.m[[index]][5] - l.m[[index]][1])/4 * (x-1)			
		}

		for (x in 6 : (maximum-1)) {
			l.f[[index]][x] = l.f[[index]][5*(x-x%%5)/5] + 
							(l.f[[index]][5*((x-x%%5)/5+1)]-l.f[[index]][5*(x-x%%5)/5])/5 * (x%%5)
			l.m[[index]][x] = l.m[[index]][5*(x-x%%5)/5] + 
							(l.m[[index]][5*((x-x%%5)/5+1)]-l.m[[index]][5*(x-x%%5)/5])/5 * (x%%5)
		}
	}
	
	
	  # Based on census data GRR for women                                  
	  Ar <- c(29.2, 29.2, 28.6, 28.1, 27.3, 26.7, 26.4, 26.0, 26.0, 26.6, 27.4, 28.3)
	  std <- c(6.8, 6.8, 6.9, 6.9, 6.6, 6.3, 6.1, 5.9, 5.7, 6.1, 6.2, 6.3)
	  GRR <- c(1937, 1745, 1415, 1121, 1123, 1508, 1782, 1210, 897, 1015, 1003, 942)/1000
	  
	  # Reported in national vital statistics report volume 63, number 7. November 2014. 
	  e0.m <- c(47.88, 49.86, 55.50, 57.71, 61.60, 65.47, 66.80, 67.04, 70.11, 71.83, 74.13, 76.20) 
	  e0.f <- c(50.70, 53.24, 57.40, 60.90, 65.89, 70.96, 73.24, 74.64, 77.62, 78.81, 79.47, 81.04)
	  e25.m <- c(38.38, 38.59, 41.11, 40.79, 42.51, 44.36, 45.19, 45.07, 47.37, 48.67, 50.57, 52.44) 
	  e25.f <- c(39.92, 40.69, 41.86, 43.11, 45.87, 48.99, 50.79, 51.80, 54.16, 55.03, 55.42, 56.87)
	  # To estimate e0, we use e0 = sum_i=0 to 100(L_i)/l_0; L_i = (l_i+l_i+1)/2
	  # To estimate e25, we use e25 = sum_i=25 to 100(L_i)/l_25;
	  
	  # e0.m <- rep(0,12); e0.f <- rep(0,12)
	  # e25.m <- rep(0,12); e25.f <- rep(0,12)  
	  # for (period in 1 : 12) {
	  #		e0.m[period] = e0.m[period]+ (1+l.m[[period]][1])/2
	  #	    e0.f[period] = e0.f[period]+ (1+l.f[[period]][1])/2
	  #	for (i in 1 : 99) {
	  #		e0.m[period] = e0.m[period]+ (l.m[[period]][i]+l.m[[period]][i+1])/2
	  #		e0.f[period] = e0.f[period]+ (l.f[[period]][i]+l.f[[period]][i+1])/2
	  #	}
	  #		e0.m[period] = e0.m[period] + l.m[[period]][100]/2
	  #		e0.f[period] = e0.f[period] + l.f[[period]][100]/2
	  # }
	  
	  # for (period in 1 : 12) {
	  #	for (i in 25 : 99) {
	  #		e25.m[period] = e25.m[period]+ (l.m[[period]][i]+l.m[[period]][i+1])/2
	  #		e25.f[period] = e25.f[period]+ (l.f[[period]][i]+l.f[[period]][i+1])/2
	  #	}
	  #		e25.m[period] = e25.m[period] + l.m[[period]][100]/2
	  #		e25.f[period] = e25.f[period] + l.f[[period]][100]/2
			
	  #		e25.m[period] = e25.m[period]/l.m[[period]][25]
	  #		e25.f[period] = e25.f[period]/l.f[[period]][25]
	  # }
	  
	  NRR.f <- GRR * c(l.f[[1]][Ar[1]], l.f[[2]][Ar[2]], l.f[[3]][Ar[3]], l.f[[4]][Ar[4]], l.f[[5]][Ar[5]], l.f[[6]][Ar[6]],
				   l.f[[7]][Ar[7]], l.f[[8]][Ar[8]], l.f[[9]][Ar[9]], l.f[[10]][Ar[10]], l.f[[11]][Ar[11]], l.f[[12]][Ar[12]])
	  NRR.m <- GRR * c(l.m[[1]][Ar[1]], l.m[[2]][Ar[2]], l.m[[3]][Ar[3]], l.m[[4]][Ar[4]], l.m[[5]][Ar[5]], l.m[[6]][Ar[6]],
				   l.m[[7]][Ar[7]], l.m[[8]][Ar[8]], l.m[[9]][Ar[9]], l.m[[10]][Ar[10]], l.m[[11]][Ar[11]], l.m[[12]][Ar[12]])
	  #sigma.NRR.f <- sqrt(2*(abs(r.f)*Ar-log(NRR.f))/(abs(r.f)^2))
	  #sigma.NRR.m <- sqrt(2*(abs(r.m)*Ar-log(NRR.m))/(abs(r.m)^2))
	  

	for (index in 1 : 12) {
		for (age in alpha.f : beta.f) {
			# m.f[[index]][age] <- GRR[index] / (beta.f-alpha.f) # fertility rate based on a uniform distribution assumption 
			 m.f[[index]][age] <- GRR[index] * dnorm(age, mean = Ar[index], sd = std[index], log = FALSE)/
			 (pnorm(beta.f, mean = Ar[index], sd = std[index])-pnorm(alpha.f, mean = Ar[index], sd = std[index])) # normal distribution with fixed sd
			# m.f[[index]][age] <- NRR.f[index] * dnorm(age, mean = Ar, sd = sigma.NRR.f[index], log = FALSE)/
			# (pnorm(beta.f, mean = Ar, sd = sigma.NRR.f[index])-pnorm(alpha.f, mean = Ar, sd = sigma.NRR.f[index]))/l.f[[index]][age] #truncated normal distribution
		}
		for (age in alpha.m : beta.m) {
			# m.m[[index]][age] <- GRR[index] / (beta.m-alpha.m) # uniform distribution 
			 m.m[[index]][age] <- GRR[index] * dnorm(age, mean = Ar[index], sd = std[index], log = FALSE)/
			 (pnorm(beta.m, mean = Ar[index], sd = std[index])-pnorm(alpha.m, mean = Ar[index], sd = std[index])) # normal distribution with fixed sd
			# m.m[[index]][age] <- NRR.m[index] * dnorm(age, mean = Ar, sd = sigma.NRR.m[index], log = FALSE)/
			# (pnorm(beta.m, mean = Ar, sd = sigma.NRR.m[index])-pnorm(alpha.m, mean = Ar, sd = sigma.NRR.m[index]))/l.m[[index]][age] #truncated normal distribution
		}
	}
	# I choose 6 for the sd of fertility rate because the range of sd for model fertility schedules 
	# in Coale and Trussell (1974) is from 5 to 7.5. 
  # By fixing sd, we control away both fertility level effect (in GRR) and 
	# fertility timing effect (in mean childbearing age and variance in timing). 	
	
	# Estimate the intrinsic growth rate based on l(x), m(x) and the characteristic function. 
	r.f <- c(); r.m <- c()
	
	for (index in 1 : 12) {
		char.func.f <- function(r) { 
						frac = 0 
						for (x in 1:100) { 
							frac = frac + l.f[[index]][x]*m.f[[index]][x]*exp(-r*x)
						}
						frac = frac -1
				   }
		char.func.m <- function(r) { 
						frac = 0 
						for (x in 1:100) { 
							frac = frac + l.m[[index]][x]*m.m[[index]][x]*exp(-r*x)
						}
						frac = frac -1
				   }
		r.f[index] <- uniroot(char.func.f, c(-1,1))$root
		r.m[index] <- uniroot(char.func.m, c(-1,1))$root
	}
	
	obs.r <- c(21.0, 21.0, 15.0, 16.2, 7.3, 14.5, 18.5, 13.3, 11.5, 9.8, 13.2, 9.7)/1000 #observed annual r for the population