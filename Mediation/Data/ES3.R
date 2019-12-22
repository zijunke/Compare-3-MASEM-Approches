# cohen d to correlation
# cohen d: the standardized between-group difference 
#		   in mean change scores
d2r = function(d,n1,n2){
	a = ((n1+n2)^2)/n1/n2
	r = d/sqrt((d^2) + a)
	return(r)
}

# group difference in mean change score and its CI to correlation
# group difference in mean change score : 
#		the standardized between-group difference 
#		in mean change scores
mdCI2r = function(mdd,CIup,CIlow,alpha,n1,n2){
	sp	= (CIup-CIlow)/2/qt(1-alpha/2,n1+n2-2)/sqrt(1/n1+1/n2)
	r	= d2r(mdd/sp,n1,n2)
	return(r)
}


# mean change scores and their SDs to correlation
# mean change score: Post-Pre or change score or improvement
# SDs: standard deviations or CI
msd2r = function(md1,md2,SD1,SD2,n1,n2){
	a = n1*n2/(n1+n2)
	mdd = md2-md1
	SDT = (SD1^2)*(n1-1)+(SD2^2)*(n2-1)+a*(mdd^2)
	SDT = sqrt(SDT/(n1+n2-1))
	r = mdd/SDT*sqrt(a/(n1+n2-1))
	return(r)
}



# F statistic to correlation
# F statistic: Time by group interaction
# 2*2 design; 2*3 design does not apply
F2r = function(F,n1,n2){
	d = sqrt(F*(n1+n2)/n1/n2)
	return(d2r(d,n1,n2))
}

# Paired sample cohen ds to correlation
# Cohen d: The two within-group 
#			standardized mean change scores
pd2r = function(d1,d2,md1,md2,n1,n2){
	SD1 = md1/d1
	SD2 = md2/d2
	a = n1*n2/(n1+n2)
	mdd =md2-md1
	SDT = (SD1^2)*(n1-1)+(SD2^2)*(n2-1)+a*(mdd^2)
	SDT = sqrt(SDT/(n1+n2-1))
	r = mdd/SDT*sqrt(a/(n1+n2-1))
	return(r)
}



# calculate means and sds for a total sample
# from two subgroups
gm2tm = function(m1,m2,n1,n2){
	return((m1*n1+m2*n2)/(n1+n2))
}

gsd2tsd = function(sd1,sd2,n1,n2){
	sdt = sqrt(((sd1^2)*(n1-1)+(sd2^2)*(n2-1))/(n1+n2-1))
	return(sdt)
}

#-----------------------------------------
# chen
#-------------------------------------------
#Duarte2016
F2r(3.59,29,19)
-F2r(0.49,29,19)

t1depr.m = gm2tm(2.83,3.37,35,36) # time 1 depression total mean 
t1depr.sd = gsd2tsd(2.29,4.45,35,36) # time 1 depression total sd
t1depr.m/42
t1depr.sd/42

#Forkman2014
F2r(0.089,64,66)
-F2r(0.409,64,66)

t1depr.m = gm2tm(10.27,10.21,64,66) # time 1 depression total mean 
t1depr.sd = gsd2tsd(3.69,3.55,64,66) # time 1 depression total sd
t1depr.m/52
t1depr.sd/52

#CladderMicus2018 #ANCOVA
F2r(35.702,63,52)
-F2r(9.26,63,52)

# time 1 depression total mean 
 # time 1 depression total sd
7.75/52
5.08/52

#Batink,2013 # Table 2
msd2r(3.1,14.8,9,16.8,66,64)
msd2r(-0.5,-3.2,4.3,4.7,66,64) #Hamilton Depression Rating Scale (HDRS).
msd2r(-3.1,-7.7,8.4,8.8,66,64) #Inventory of Depressive Symptoms
(50+48)/130

# time 1 depression total mean 
 # time 1 depression total sd
t1depr.m = gm2tm(10.3,10.2,64,66) # time 1 depression total mean 
t1depr.sd = gsd2tsd(3.7,3.6,64,66) # time 1 depression total sd
t1depr.m/52
t1depr.sd/52
 
t1depr.m = gm2tm(22.4,22.5,64,66) # time 1 depression total mean 
t1depr.sd = gsd2tsd(10.7,8.7,64,66) # time 1 depression total sd
t1depr.m/84
t1depr.sd/84

# Gayner2012 # 2*3 interactions
F2r(5.9,78,39)
-F2r(0.8,78,39)

gm2tm(42.9,45.5,78,39)
gsd2tsd(7.1,6.7,78,39)

t1depr.m = gm2tm(8.5,8.6,78,39) # time 1 depression total mean 
t1depr.sd = gsd2tsd(3.8,3.1,78,39) # time 1 depression total sd
t1depr.m/21
t1depr.sd/21

#Foley2010
F2r(18.51,55,60)
-F2r(18.78,55,60)
gm2tm(0.764,0.783,55,60)

t1depr.m = gm2tm(16.02,14.38,55,60) # time 1 depression total mean 
t1depr.sd = gsd2tsd(7.28,8.12,55,60) # time 1 depression total sd
t1depr.m/52
t1depr.sd/52

#Raedt2011
F2r(4.31,44,26)
-F2r(8.51,44,26)

gm2tm(45.2,45,45,26)
gsd2tsd(9.8,8.6,45,26)
(33+19)/(33+12+26)

t1depr.m = gm2tm(17.64,7.88,44,26) # time 1 depression total mean 
t1depr.sd = gsd2tsd(10.34,8.17,44,26) # time 1 depression total sd
t1depr.m/63
t1depr.sd/63

#Eisendrath2015
msd2r(-5.2,-12.2,3.9,5.3,17,19)
msd2r(-8.4,-11.2,6.7,4.9,17,19)
gm2tm(35.3,37,23,20)
gsd2tsd(10,10.8,23,20)
gm2tm(.78,.65,23,20)

t1depr.m = gm2tm(20.1,19.3,23,20) # time 1 depression total mean 
t1depr.sd = gsd2tsd(3.3,3.7,23,20) # time 1 depression total sd
t1depr.m/52
t1depr.sd/52

t1depr.m = gm2tm(21.8,14.1,23,20) # time 1 depression total mean 
t1depr.sd = gsd2tsd(5.5,4.0,23,20) # time 1 depression total sd
t1depr.m/48
t1depr.sd/48

#Frank2015 # Table 3
msd2r(-0.31,0.43,0.63,0.87,18,18) #SCS mindfulness
msd2r(-0.15,6.04,2.97,3.71,18,18) # FFM observe
msd2r(-0.47,-1.39,3.73,3.35,18,18)

t1depr.m = gm2tm(3.56,3.17,18,18) # time 1 depression total mean 
t1depr.sd = gsd2tsd(4.06,3.56,18,18) # time 1 depression total sd
t1depr.m/24
t1depr.sd/24

#Bränström2010
-F2r(2.91,32,39)
msd2r(-0.69,-1.74,2.78,3.01,39,32) # Table 2 not used because the precision is lower

t1depr.m = gm2tm(6.41,7.18,32,39) # time 1 depression total mean 
t1depr.sd = gsd2tsd(4.46,3.55,32,39) # time 1 depression total sd
t1depr.m/21
t1depr.sd/21

#Armstrong2016 #ANCOVA
F2r(3.27,17,17)
-F2r(0.04,17,17)

gm2tm(29.4,29.7,17,17)
gsd2tsd(8.7,8.4,17,17)
(16+15)/34

t1depr.m = gm2tm(6.6,11.1,17,17) # time 1 depression total mean 
t1depr.sd = gsd2tsd(5.6,4.9,17,17) # time 1 depression total sd
t1depr.m/27
t1depr.sd/27

#-----------------------------------------
# Liu
#-----------------------------------------
#Labelle2010
F2r(12.81,46,31)
-F2r(5.25,46,31)
a = 0.35
b = -0.19
c = -0.23
cp = c-a*b
rMY = a*c+b*(1-a^2)

t1depr.m = gm2tm(10.65,7.55,46,31) # time 1 depression total mean 
t1depr.sd = gsd2tsd(7.06,4.98,46,31) # time 1 depression total sd
t1depr.m/30
t1depr.sd/30

#McManus2012
F2r(3.66,36,38) # 2*3 design # not used
-F2r(1.37,36,38) # 2*3 design # not used
pd2r(-0.84/sqrt(38),3.62/sqrt(36),36.9-38.5,38.9-33.9,38,36) #post-pre
pd2r(1.42/sqrt(38),0.32/sqrt(36),36.9-38.5,38.9-33.9,38,36) #1yfu-pre

gm2tm(41.28,43.92,36,38)
gsd2tsd(11.98,10.98,36,38)
(27+31)/74

t1depr.m = gm2tm(19.89,19.74,36,38) # time 1 depression total mean 
t1depr.sd = gsd2tsd(13.96,10.15,36,38) # time 1 depression total sd
t1depr.m/63
t1depr.sd/63

#James2018 # ANCOVA
# post
F2r(16.5,28,32)
-F2r(0.5,28,32)

eta2 = 0.22
d = 2*sqrt(eta2/(1-eta2))
d2r(d,28,32)
eta2 = 0.01
d = 2*sqrt(eta2/(1-eta2))
d2r(d,28,32)

# 10 week follow up
F2r(7.6,28,32)
-F2r(0.5,28,32)

eta2 = 0.12
d = 2*sqrt(eta2/(1-eta2))
d2r(d,28,32)

(23+26)/60

t1depr.m = gm2tm(16.1,14.0,28,32) # time 1 depression total mean 
t1depr.sd = gsd2tsd(12.1,10.3,28,32) # time 1 depression total sd
t1depr.m/42
t1depr.sd/42

#Hazlett-Stevens2017
F2r(8.9,25,43)
-F2r(16.4,25,43)

51/(51+16+1)=0.75

t1depr.m = gm2tm(8.7,7.7,25,43) # time 1 depression total mean 
t1depr.sd = gsd2tsd(9.1,8.0,25,43) # time 1 depression total sd
t1depr.m/42
t1depr.sd/42

#Helmes2017
# post-pre
d1 = 1.11/sqrt(26)
d2 = 4.39/sqrt(26)
pd2r(d1,d2,1.00,6.54,26,26)

# followup-pre
d1 = .90/sqrt(26)
d2 = 5.59/sqrt(26)
pd2r(d1,d2,.72,10.30,26,26)

34/52

#Key2017
F2r(5.21,18,18)
-F2r(13.52,18,18)
d2r(0.77,18,18)
d2r(-1.25,18,18)
msd2r(-0.19,6.91,8.66,9.72,18,18)
msd2r(2.65,-2.17,3.39,4.27,18,18)
# Self-compassion
F2r(4.82,18,18)

c(gm2tm(40.53,46.06,18,18),
gsd2tsd(12.42,15.25,18,18))
17/36

t1depr.m = gm2tm(15.33,21.24,18,18) # time 1 depression total mean 
t1depr.sd = gsd2tsd(9.38,13.72,18,18) # time 1 depression total sd
t1depr.m/63
t1depr.sd/63


# Moore2016
# Paper-and-pencil
#F2r(1.2,32,35) # not used because the F test does not seem to be the usual one
#-F2r(0,32,35)	# the dfs of F test have decimals
d2r(0.2,32,35)
d2r(-0.1,32,35)

c(gm2tm(69.8,72,32,35),
gsd2tsd(4.9,5.5,32,35))
gsd2tsd(.72,.81,32,35)

t1depr.m = gm2tm(9.4,9.5,32,35) # time 1 depression total mean 
t1depr.sd = gsd2tsd(3.8,3.6,32,35) # time 1 depression total sd
t1depr.m/32
(t1depr.m-8)/32
t1depr.sd/32

#Hoge2015
d2r(1.81*sqrt(1/19+1/19),19,19)
msd2r(5.53,17.2,11.0,25.8,19,19)

#Kingston2015
#F2r(8.88,6,7) not used because it is only a subscale
gm2tm(49.75,50.38,8,8)
gsd2tsd(10/89,13.17,8,8)
(10)/16

t1depr.m = gm2tm(8.0,8.5,7,6) # time 1 depression total mean 
t1depr.sd = gsd2tsd(6.2,6.4,7,6) # time 1 depression total sd
t1depr.m/60
t1depr.sd/60

t1depr.m = gm2tm(11.9,11.9,8,8) # time 1 depression total mean 
t1depr.sd = gsd2tsd(2.3,1.0,8,8) # time 1 depression total sd
t1depr.m/21
t1depr.sd/21

#Hoelzel2016 # Table 2
#F2r(7.44,22,23) Not used because Age and days were included as covariates
d1 = -0.712/sqrt(22)
d2 = 3.292/sqrt(23)
pd2r(d1,d2,4.00-4.08,4.17-3.63,22,23)
c(gm2tm(35,28.9,23,23),
gsd2tsd(10,8,23,23))
(12+17)/46


#Jing2014 Table 6.3
msd2r(2.76,6.77,12.07,14.58,57,62)
msd2r(-2.85,-6.72,4.82,6.38,55,64)

msd2r(1.35,6.98,13.09,16.84,51,60)
msd2r(-1.39,-4.70,8.52,6.32,46,56)

t1depr.m = gm2tm(16.91,17.55,70,71) # time 1 depression total mean 
t1depr.sd = gsd2tsd(8.97,8.92,70,71) # time 1 depression total sd
t1depr.m/60
t1depr.sd/60

# Manicavasagar2012
d1 = -4.6/sqrt(26)
d2 = -4.8/sqrt(19)
pd2r(d1,d2,23.62-36.23,21.21-32.42,26,19)

c(gm2tm(45,47,26,19),
gsd2tsd(12.9,13.8,26,19))

t1depr.m = gm2tm(36.23,32.42,26,19) # time 1 depression total mean 
t1depr.sd = gsd2tsd(11.11,9.01,26,19) # time 1 depression total sd
t1depr.m/63
t1depr.sd/63

#Keune2011
F2r(14.60,37,40)
-F2r(4.10,37,40)
msd2r(0.73,-3.5,1.55*sqrt(37),1.41*sqrt(40),37,40) # not used because precision is lower
c(gm2tm(45.24,48.93,27,30),
gsd2tsd(10.50,9.68,27,30))
(27+30)/77

t1depr.m = gm2tm(12.70,9.05,37,40) # time 1 depression total mean 
t1depr.sd = gsd2tsd(9.19,8.60,37,40) # time 1 depression total sd
t1depr.m/63
t1depr.sd/63

#Jasbi2018
-F2r(53.22,24,24) #Table 4
pd2r(-0.26,-2.62,17.82-18.25,13.20-18.33,24,24) #Table 3 & 5 # Not used because precision is lower

c(gm2tm(50.03,52.91,24,24),
gsd2tsd(2.45,2.91,24,24))

t1depr.m = gm2tm(18.33,18.25,24,24) # time 1 depression total mean 
t1depr.sd = gsd2tsd(1.88,1.75,24,24) # time 1 depression total sd
t1depr.m/42
t1depr.sd/42

# Greenberg2017
-F2r(22.5,12,13) #BDI
-F2r(4.77,16,9) #HAM

eta2 = 0.505
d = 2*sqrt(eta2/(1-eta2))
d2r(d,12,13)

eta2 = 0.178
d = 2*sqrt(eta2/(1-eta2))
d2r(d,16,9)

c(gm2tm(39.77,36.89,22,18),
gsd2tsd(10.6,15.83,22,18))
gm2tm(0.59,.66,22,18)

#BDI
t1depr.m = gm2tm(23.51,22.18,12,13) # time 1 depression total mean 
t1depr.sd = gsd2tsd(7.09,9.39,12,13) # time 1 depression total sd
t1depr.m/63
t1depr.sd/63

#HAM
t1depr.m = gm2tm(22.80,24.00,16,9) # time 1 depression total mean 
t1depr.sd = gsd2tsd(8.62,6.34,16,9) # time 1 depression total sd
t1depr.m/52
t1depr.sd/52

#Johns2016 Tabel 5
pd2r(-0.94,-1.05,7.80-12.53,6.27-11.35,36,35)
pd2r(-.66,-.98,8.51-12.53,6.55-11.35,36,35)

gm2tm(0.943,.861,35,36)

t1depr.m = gm2tm(11.35,12.53,35,36) # time 1 depression total mean 
t1depr.sd = gsd2tsd(5.57,4.90,35,36) # time 1 depression total sd
t1depr.m/24
t1depr.sd/24

#-----------------------------------------
# Liang
#---------------------------------------------
#Wong2018 Table 2
# sp = (3.16-0.33)/2/qt(.975,137)/sqrt(1/69+1/70)
# d2r(1.74/sp,69,70) # Observe
# sp = (1.26+1.04)/2/qt(.975,137)/sqrt(1/69+1/70)
# d2r(0.19/sp,69,70) # Describe
# sp = (1.28+1.22)/2/qt(.975,137)/sqrt(1/69+1/70)
# d2r(0.03/sp,69,70) # Act with Awareness
# sp = (0.36+1.98)/2/qt(.975,137)/sqrt(1/69+1/70)
# d2r(-0.81/sp,69,70) # Nonjudege
# sp = (1.67+0.23)/2/qt(.975,137)/sqrt(1/69+1/70)
# d2r(0.72/sp,69,70) # Nonreact

-F2r(4.78,69,70)
mdCI2r(-0.77,-0.07,-1.47,.05,69,70)

t1depr.m = gm2tm(6.96,6.59,98,99) # time 1 depression total mean 
t1depr.sd = gsd2tsd(2.95,2.45,98,99) # time 1 depression total sd
t1depr.m/15
t1depr.sd/15

#Vollestad2011
# Not used because Table 6 are direct results
F2r(19.1,31,34)
-F2r(25.0,31,34)

# time 1 depression total mean 
 # time 1 depression total sd
16.9/63
8.6/63

#VanSon2013
#-F2r(5.37,70,69) # 2*3 design
#-F2r(8.38,70,69) # 2*3 design
d2r(-0.59,70,69) 
d2r(-0.71,70,69)

c(gm2tm(56,57,70,69),
gsd2tsd(13,13,70,69))
(37+32)/139

t1depr.m = gm2tm(7.9,8.9,70,69) # time 1 depression total mean 
t1depr.sd = gsd2tsd(3.8,3.9,70,69) # time 1 depression total sd
t1depr.m/21
t1depr.sd/21

t1depr.m = gm2tm(25.3,26.6,70,69) # time 1 depression total mean 
t1depr.sd = gsd2tsd(5.8,6.3,70,69) # time 1 depression total sd
t1depr.m/32
t1depr.sd/32

# Sevinc2018 # Table 2
# SCS mindfulness
pd2r(2.98/sqrt(19),2.18/sqrt(18),3.64-3.25,3.65-3.42,19,18)
c(gm2tm(42.1,35.3,22,28),
gsd2tsd(11.5,9.3,22,28))
32/50

#Song2015 # ANCOVA
F2r(5.03,23,21) #
-F2r(10.99,23,21)
c(gm2tm(19.6,19.5,21,23),
gsd2tsd(1.7,2,21,23))
36/45

t1depr.m = gm2tm(8.3,8.6,21,23) # time 1 depression total mean 
t1depr.sd = gsd2tsd(5.1,8.9,21,23) # time 1 depression total sd
t1depr.m/42
t1depr.sd/42

# Tovote2014
-F2r(9.71,31,31)
-F2r(17.52,31,31)
pd2r(-0.65/sqrt(31),-4.75/sqrt(31),23.5-24.3,17.1-23.6,31,31)
pd2r(-0.71/sqrt(31),-6.33/sqrt(31),7.1-7.5,4.7-8.9,31,31)

t1depr.m = gm2tm(23.6,24.3,31,31) # time 1 depression total mean 
t1depr.sd = gsd2tsd(7.7,8.0,31,31) # time 1 depression total sd
t1depr.m/63
t1depr.sd/63

t1depr.m = gm2tm(8.9,7.5,31,31) # time 1 depression total mean 
t1depr.sd = gsd2tsd(3.5,2.8,31,31) # time 1 depression total sd
t1depr.m/28
t1depr.sd/28


# Zhang2015
-F2r(4.469,30,30)
c(gm2tm(78.57,77.63,30,30),
gsd2tsd(2.94,3.01,30,30))
25/60

t1depr.m = gm2tm(14.20,15.97,30,30) # time 1 depression total mean 
t1depr.sd = gsd2tsd(2.5,4.17,30,30) # time 1 depression total sd
t1depr.m/30
t1depr.sd/30

# Polusny2015 # Table 2
mdCI2r(9.14,14.37,3.91,.05,58,58)
mdCI2r(-1.17,2.56,-.22,.05,58,58)

mdCI2r(9.73,15.04,4.42,.05,58,58)
mdCI2r(-1.34,2.75,-.07,.05,58,58)
# time 1 depression total mean 
# time 1 depression total sd
15/27
5.3/27

# Nathan2017 Table2
mdCI2r(-4.69,-2.43,-6.96,.05,32,30)
t1depr.m = gm2tm(12.06,11.30,30,30) # time 1 depression total mean 
sd1 = (13.68-10.44)/2/1.96*sqrt(32)
sd2 = (12.95-9.65)/2/1.96*sqrt(30)
t1depr.sd = gsd2tsd(sd1,sd2,32,30) # time 1 depression total sd
t1depr.m/27
t1depr.sd/27

#Reich2017 # Removed
# Figure 2
#0.33*0.03
#0.17*(-0.76)
#0.33*0.03*0.17*(-0.76)

#Schellekens2017
# pooled post and follow-ups and adjusted for baseline
mdCI2r(13.94,22.11,5.77,.05,21,18)
d2r(0.84,21,18)
 
mdCI2r(-5.19,-1.43,-8.94,.05,21,18)
d2r(-.69,21,18)

c(gm2tm(60.6,57.0,31,32),
gsd2tsd(6.8,8.5,31,32))
(18+15)/63

t1depr.m = gm2tm(11.97,12.94,31,32) # time 1 depression total mean 
t1depr.sd = gsd2tsd(6.97,8.05,31,32) # time 1 depression total sd
t1depr.m/21
t1depr.sd/21

# VanDijk2017
d2r(0.35,67,58)
c(gm2tm(23.3,23.7,84,83),
gsd2tsd(1.77,1.91,84,83))
(71+60)/167























