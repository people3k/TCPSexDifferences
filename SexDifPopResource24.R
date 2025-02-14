####Analysis script for ``The effects of resources and population expansion on
#female-male protein consumption among hunter-gatherers#
###Set working directory

###Load packages
library(splines)
library(reshape2)
library(ggplot2)
library(nlme)
library(rcarbon)
library(cowplot)
library(rcarbon)
library(effects)
library(BayesFactor)

###########Analysis of Isotope data from Texas===========================================

####Load data for Central Texas and Texas Coastal Plain HG Sex Difference

d2 <- read.csv("data/TexasHGIsotopesMaleFemale.csv")

##Calculate median nitrogen by site for the entire data set
medianN1<-aggregate((d2$delta13ccol), list(d2$site), FUN=median)
medianN1

#Partition the data by biogeographic setting and sex identification. In this case, we partition to Riverine ecosystems between 7800 and 400 cal BP

MM<-subset(d2, sex=="M" & Setting2=="Riverine" & PeriodID>0 & delta15n>0 & delta13ccol<0)
FF<-subset(d2, sex=="F" & Setting2=="Riverine" & PeriodID>0 & delta15n>0 & delta13ccol<0)
d3<-rbind(MM,FF)

##Calculate median nitrogen by site for the Riverine Zone
medianN1<-aggregate((d3$delta13ccol), list(d3$site), FUN=median)
medianN1


# Perform Bayesian t-test on delta 15N by sex identification in Riverine ecosystems
bayesian_ttest <- ttestBF(x = MM$delta15n, y = FF$delta15n, rscale=.707)
# Show the result
print(bayesian_ttest)
plot(bayesian_ttest)

###run Mann-Whittney U for delta 13C col and 15N col by sex for the Riverine ecosystems
wilcox.test(delta15n ~ sex, data=d3)

#wilcox.test(delta13ccol ~ sex, data=d3)
#kruskal.test(delta13ccol~sex, data = d3)

#Calculate median nitrogen values by sex
meansC1<-aggregate((d3$delta15n), list(d3$sex), FUN=median)
meansC1

#Create violin plot of the data sex vs. delta 15N. Note the Bayes Factors and Mann-Whittney U
pn1 <- ggplot(d3, aes(factor(sex), (delta15n), fill=sex))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  scale_y_continuous(limits=c(5,17.5), breaks=c(6,8,10,12,14,16))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Sex ID (Female, Male)", y="Delta 15N collagen", title = "B. Riverine HG Delta 15 N Collagen by Sex")+
  #facet_wrap(~factor(Setting2))+
  annotate("text", x =1.5, y = 15.5, label = "W=1276.5**", size = 6)+
  annotate("text", x =1.5, y = 16.5, label = "BF=45.81", size = 6)+
annotate("text", x =.95, y = 15.5, label = "female n=60", size = 6)+
  annotate("text", x =.95, y = 14.5, label = "Med=10.45", size = 6)+
annotate("text", x =2.05, y = 15.5, label = "male n=75", size = 6)+
  annotate("text", x =2.05, y = 14.5, label = "Med=11.10", size = 6)
pn1

###Inland: Run the same analysis for inland mortuary sites==============================================================================================
#Partition the data by inland cemeteries
MMi<-subset(d2, sex=="M" & Setting2=="Inland" & PeriodID>0 & delta15n>0)
FFi<-subset(d2, sex=="F" & Setting2=="Inland" & PeriodID>0 & delta15n>0)
d3i<-rbind(MMi,FFi)

##Calculate median nitrogen by site for the inland Zone
medianN1<-aggregate((d3i$delta13ccol), list(d3i$site), FUN=median)
medianN1

# Perform Bayesian t-test
bayesian_ttesti <- ttestBF(x = MMi$delta15n, y = FFi$delta15n, rscale=.707)

# Show the result
print(bayesian_ttesti)
plot(bayesian_ttesti)

###run Mann-Whitney U for 15N col by sex for the inland ecosystems
wilcox.test(delta15n ~ sex, data=d3i)

#Calculate median nitrogen values by sex
meansC1<-aggregate((d3i$delta15n), list(d3i$sex), FUN=median)
meansC1

pn2 <- ggplot(d3i, aes(factor(sex), (delta15n), fill=sex))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  scale_y_continuous(limits=c(5,17.5), breaks=c(6,8,10,12,14, 16))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Sex ID (Female, Male)", y="Delta 15N collagen", title = "C. Inland HG Delta 15 N Collagen by Sex")+
  #facet_wrap(~factor(Setting2))+
  annotate("text", x =1.5, y = 14, label = "W=21**", size = 6)+
  annotate("text", x =1.5, y = 15.5, label = "BF=0.96", size = 6)+
  annotate("text", x =.95, y = 15.5, label = "female n=10", size = 6)+
  annotate("text", x =.95, y = 14.5, label = "Med=8.74", size = 6)+
  annotate("text", x =2.25, y = 15.5, label = "male n=13", size = 6)+
  annotate("text", x =2.25, y = 14.5, label = "Med=9.5", size = 6)
#annotate("text", x =19.5, y = 34.5, label = "Tenure at USU", size = 5, angle=90)
pn2

###Run the same analysis for the Coast===================

MMc<-subset(d2, sex=="M" & Setting2=="Coast" & PeriodID>0 & delta15n>0)
FFc<-subset(d2, sex=="F" & Setting2=="Coast" & PeriodID>0 & delta15n>0)
d3c<-rbind(MMc,FFc)

##Calculate median nitrogen by site for the coastal Zone
medianN1<-aggregate((d3c$delta13ccol), list(d3c$site), FUN=median)
medianN1

# Perform Bayesian t-test
bayesian_ttestc <- ttestBF(x = MMc$delta15n, y = FFc$delta15n, rscale=.707)

# Show the result
print(bayesian_ttestc)
plot(bayesian_ttestc)

###run Mann-Whitney U for delta 15N col by sex in the coastal ecosystems

wilcox.test(delta15n ~ sex, data=d3c)

#Calculate median nitrogen values by sex
meansC1<-aggregate((d3c$delta15n), list(d3c$sex), FUN=median)
meansC1

##Creat violin plot
pn3 <- ggplot(d3c, aes(factor(sex), (delta15n), fill=sex))+
  geom_violin()+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  stat_boxplot(geom ='errorbar')+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  theme_bw() +
  scale_y_continuous(limits=c(5,17.5), breaks=c(6,8,10,12,14,16))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Sex ID (Female, Male)", y="Delta 15N collagen", title = "A. Coastal HG Delta 15 N Collagen by Sex")+
  #facet_wrap(~factor(Setting2))+
  annotate("text", x =1.5, y = 14, label = "W=700", size = 6)+
  annotate("text", x =1.5, y = 15.5, label = "BF=0.51", size = 6)+
  annotate("text", x =.95, y = 15.5, label = "female n=38", size = 6)+
  annotate("text", x =.95, y = 14.5, label = "Med=11.22", size = 6)+
  annotate("text", x =2.05, y = 15.5, label = "male n=43", size = 6)+
  annotate("text", x =2.05, y = 14.5, label = "Med=11.36", size = 6)
#annotate("text", x =19.5, y = 34.5, label = "Tenure at USU", size = 5, angle=90)
pn3


##Stack the plots of the coastal, riverine and inland ecosystem analytical units
Fig5Rev<-plot_grid(pn3,pn1,pn2, ncol=3, align="hv", axis = "rl")
Fig5Rev

###Export the stacked graphs as a pdf

pdf("graphics/TexasGeoSerxDifAaotated.pdf", width=18.55, height=14.55)
Fig5Rev
dev.off()


##Changes in Sex differences by demographic state/Time========================================
#We have enough data to analyze the riverine ecosystems over time and assess changes in sex differences
#relative to changes in estimated population density and per capita growth

#subset the riverine data with males and females by time periods 1,2,3,4,and 5.
MM2<-subset(MM, PeriodID=="5")
FF2<-subset(FF, PeriodID=="5")
MM2<-subset(MM, PeriodID>3)
FF2<-subset(FF, PeriodID>3)
PerFiveRiver<-rbind(MM2, FF2)
# Perform Bayesian t-test on delta 15N by sex identification in Riverine ecosystems
bayesian_ttest <- ttestBF(x = MM2$delta15n, y = FF2$delta15n, rscale=.707)

# Show the result
print(bayesian_ttest)
#plot(bayesian_ttest)

##Run Mann-Whittney U test for each period
wilcox.test(delta15n ~ sex, data=PerFiveRiver)

#Calculate median nitrogen values by sex and cultural historical period
meansC1<-aggregate((PerFiveRiver$delta15n), list(PerFiveRiver$sex), FUN=median)
meansC1

##Remove the one individual with a sex identification between 7200 and 5801 cal BP.
MMs<-subset(d2, sex=="M" & Setting2=="Riverine" & PeriodID>0.1 & delta15n>0 & delta13ccol<0)
FFs<-subset(d2, sex=="F" & Setting2=="Riverine" & PeriodID>0.1 & delta15n>0 & delta13ccol<0)
d3s<-rbind(MMs,FFs)

#Calculate median nitrogen values by cultural historical periods.
meansC1<-aggregate((d3s$delta15n), list(d3s$sex), FUN=median)
meansC1

#Plot nitrogen isotope values by sex and time periods
ptime1 <- ggplot(d3s, aes(factor(PeriodID), (delta15n), fill=sex))+
  geom_violin(trim = FALSE, position = position_dodge(width = 0.75)) +  # Separate violins for each gender
  stat_summary(fun = median, geom = "point", size = 4, color = "black", position = position_dodge(width = 0.75)) +  # Median points
  stat_boxplot(geom ='errorbar', position = position_dodge(width = 0.75))+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  scale_y_continuous(limits=c(6,16), breaks=c(6,8,10,12,14,16))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Time period", y="Delta 15N collagen", title = "B. Riverine Delta 15N Collagen by Sex and Period")+
 # facet_wrap(~factor(FloodID))+
  geom_hline(yintercept = 10.45)+
  annotate("text", x =1, y = 16, label = "BF=1.35", size = 6)+
  annotate("text", x =1, y = 15.5, label = "W=10*", size = 6)+
  annotate("text", x =1, y = 7.5, label = "female n=8", size = 6)+
  annotate("text", x =1, y = 7, label = "Med=11.04", size = 6)+
  annotate("text", x =1, y = 6.5, label = "male n=6", size = 6)+
  annotate("text", x =1, y = 6, label = "Med=11.73", size = 6)+
#2
  annotate("text", x =2, y = 16, label = "BF=0.74", size = 6)+
  annotate("text", x =2, y = 15.5, label = "W=25", size = 6)+
  annotate("text", x =2, y = 7.5, label = "female n=7", size = 6)+
  annotate("text", x =2, y = 7, label = "Med=10.29", size = 6)+
  annotate("text", x =2, y = 6.5, label = "male n=10", size = 6)+
  annotate("text", x =2, y = 6, label = "Med=10.61", size = 6)+
#3
    annotate("text", x =3, y = 16, label = "BF=0.66", size = 6)+
  annotate("text", x =3, y = 15.5, label = "W=445*", size = 6)+
  annotate("text", x =3, y = 7.5, label = "female n=30", size = 6)+
  annotate("text", x =3, y = 7, label = "Med=10.45", size = 6)+
  annotate("text", x =3, y = 6.5, label = "male n=47", size = 6)+
  annotate("text", x =3, y = 6, label = "Med=11.00", size = 6)+
#4
  annotate("text", x =4, y = 16, label = "BF=1.58", size = 6)+
  annotate("text", x =4, y = 15.5, label = "W=2.5*", size = 6)+
  annotate("text", x =4, y = 7.5, label = "female n=3", size = 6)+
  annotate("text", x =4, y = 7, label = "Med=10.5", size = 6)+
  annotate("text", x =4, y = 6.5, label = "male n=7", size = 6)+
  annotate("text", x =4, y = 6, label = "Med=11.00", size = 6)+
#5
  annotate("text", x =5, y = 16, label = "BF=5.41", size = 6)+
  annotate("text", x =5, y = 15.5, label = "W=.001**", size = 6)+
  annotate("text", x =5, y = 7.5, label = "female n=12", size = 6)+
  annotate("text", x =5, y = 7, label = "Med=10.4", size = 6)+
  annotate("text", x =5, y = 6.5, label = "male n=4", size = 6)+
  annotate("text", x =5, y = 6, label = "Med=11.32", size = 6)
ptime1

###Plot delta 13 collagen by sex=========================
##Calculate median for all individuals
median(d3s$delta13ccol)

##subset the data by time period
pp_a<-subset(d3s, PeriodID=="4")
pp_b<-subset(d3s, PeriodID=="5")
colperiod<-rbind(pp_a, pp_b)

# Perform Bayesian t-test on delta 15N by sex identification in Riverine ecosystems
bayesian_ttest <- ttestBF(x = pp_a$delta13ccol, y = pp_b$delta13ccol, rscale=.707)

# Show the result
print(bayesian_ttest)
#plot(bayesian_ttest)

##Run Mann-Whittney U test for each period
wilcox.test(delta13ccol ~ PeriodID, data=colperiod)

#Plot all individuals, including unidentified sex, by time period and pre vs. post territoriality (phase)
ptimec2 <- ggplot(d3s, aes(factor(PeriodID), (delta13ccol), fill=Phase))+
  geom_violin(trim = FALSE, position = position_dodge(width = 0.75)) +  # Separate violins for each gender
  stat_summary(fun = median, geom = "point", size = 4, color = "black", position = position_dodge(width = 0.75)) +  # Median points
  stat_boxplot(geom ='errorbar', position = position_dodge(width = 0.75))+
  #facet_wrap( ~ factor(CompID1))+
  scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  scale_y_continuous(limits=c(-22,-7), breaks=c(-20,-18,-16,-14, -12,-10,-8))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Time period", y="Delta 13C collagen", title = "A. Riverine Delta 13C Collagen by Territory Regime and Period")+
  #facet_wrap(~factor(FloodID))
  geom_hline(yintercept = -19)+
  annotate("text", x =1.5, y = -14, label = "BF=7.63", size = 6)+
  annotate("text", x =1.5, y = -14.5, label = "W=51.5**", size = 6)+
  annotate("text", x =2.5, y = -14, label = "BF=7.43", size = 6)+
  annotate("text", x =2.5, y = -14.5, label = "W=1006**", size = 6)+
  annotate("text", x =3.5, y = -14, label = "BF=0.33", size = 6)+
  annotate("text", x =3.5, y = -14.5, label = "W=327", size = 6)+
  annotate("text", x =4.5, y = -14, label = "BF=1.23", size = 6)+
  annotate("text", x =4.5, y = -14.5, label = "W=112.5*", size = 6)
ptimec2

#Calculate median carbon collagen values by cultural historical periods.
meansC1<-aggregate((d3s$delta13ccol), list(d3s$sex), FUN=median)
meansC1

#Plot individuals with sex ids by time period and pre vs. post territoriality (phase)

ptimec1 <- ggplot(d3s, aes(factor(PeriodID), (delta13ccol), fill=sex))+
  geom_violin(trim = FALSE, position = position_dodge(width = 0.75)) +  # Separate violins for each gender
  stat_summary(fun = median, geom = "point", size = 4, color = "black", position = position_dodge(width = 0.75)) +  # Median points
  stat_boxplot(geom ='errorbar', position = position_dodge(width = 0.75))+
  #facet_wrap( ~ factor(CompID1))+
  scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  scale_y_continuous(limits=c(-22,-7), breaks=c(-20,-18,-16,-14, -12,-10,-8))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Time period", y="Delta 13C collagen", title = "B. Riverine Delta 13C Collagen by Sex and Period")+
  #facet_wrap(~factor(FloodID))
  geom_hline(yintercept = -19.1)
ptimec1

##Stack the plots of riverine carbon collage
Fig5Rev<-plot_grid(ptimec2,ptimec1, ncol=1, align="hv", axis = "rl")
Fig5Rev

###Export the stacked graphs as a pdf

pdf("graphics/RiverineCollagen.pdf", width=18.55, height=14.55)
Fig5Rev
dev.off()

##Coastal ecosystems over time===========================================================

#For statistical analysis subset to periods 4 and 5 because period 3 does not have males and females identified
#Very few burials on the coast prior to period 3.
##Only one female during period 3, thus we cannot run a statistical comparison of males and females during this period
d3c<-subset(d3c, MidPoint<2101)

MMc2<-subset(MMc, PeriodID=="5")
FFc2<-subset(FFc, PeriodID=="5")
CoastSex<-rbind(MMc2, FFc2)

#Calculate medians by sex and cultural historical time period on the coast
meansC2<-aggregate((CoastSex$delta15n), list(CoastSex$sex), FUN=median)
meansC2
# Perform Bayesian t-test on delta 15N by sex identification in Riverine ecosystems
bayesian_ttest <- ttestBF(x = MMc2$delta15n, y = FFc2$delta15n, rscale=.707)

# Show the result
print(bayesian_ttest)
#plot(bayesian_ttest)

##Run Mann-Whittney U test for each period
wilcox.test(delta15n ~ sex, data=CoastSex)

###Plot delta 15N over time
ptime2 <- ggplot(d3c, aes(factor(PeriodID), (delta15n), fill=sex))+
  geom_violin(trim = FALSE, position = position_dodge(width = 0.75)) +  # Separate violins for each gender
  stat_summary(fun = median, geom = "point", size = 4, color = "black", position = position_dodge(width = 0.75)) +  # Median points
  stat_boxplot(geom ='errorbar', position = position_dodge(width = 0.75))+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  scale_y_continuous(limits=c(5,16), breaks=c(6,8,10,12,14,16))+
  scale_x_discrete(limits=c("1","2","3","4","5"))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Time period", y="Delta 15N collagen", title = "B. Coastal Delta 15N Collagen by Sex and Period")+
 # facet_wrap(~factor(sex))+
  geom_hline(yintercept = 11.22)+
  annotate("text", x =4, y = 16, label = "BF=0.51", size = 6)+
  annotate("text", x =4, y = 15.5, label = "W=30.5", size = 6)+
  annotate("text", x =4, y = 7.5, label = "female n=13", size = 6)+
  annotate("text", x =4, y = 7, label = "Med=11.38", size = 6)+
  annotate("text", x =4, y = 6.5, label = "male n=7", size = 6)+
  annotate("text", x =4, y = 6, label = "Med=12.06", size = 6)+
  #
  annotate("text", x =5, y = 16, label = "BF=0.51", size = 6)+
  annotate("text", x =5, y = 15.5, label = "W=.289.5", size = 6)+
  annotate("text", x =5, y = 7.5, label = "female n=24", size = 6)+
  annotate("text", x =5, y = 7, label = "Med=11.15", size = 6)+
  annotate("text", x =5, y = 6.5, label = "male n=28", size = 6)+
  annotate("text", x =5, y = 6, label = "Med=11.20", size = 6)
ptime2

###Plot collagen values for coastal burials
ptime2c <- ggplot(d3c, aes(factor(PeriodID), (delta13ccol), fill=sex))+
  geom_violin(trim = FALSE, position = position_dodge(width = 0.75)) +  # Separate violins for each gender
  stat_summary(fun = median, geom = "point", size = 4, color = "black", position = position_dodge(width = 0.75)) +  # Median points
  stat_boxplot(geom ='errorbar', position = position_dodge(width = 0.75))+
  #facet_wrap( ~ factor(CompID1))+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  #scale_y_continuous(limits=c(5,16), breaks=c(6,8,10,12,14,16))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Time period", y="Delta 13C collagen", title = "B. Coastal Delta 13C Collagen by Sex and Period")
  # facet_wrap(~factor(sex))+
  #geom_hline(yintercept = 10.25)
ptime2c


###Texas Radiocarbon data====================================================


###Plot mean KDE against the per capita growth rate in the North
us30pc<- read.csv("data/TCPPerCap.csv")

##Restrict the KDE to the relevant time period that corresponds with isotope data.
us30pc<-subset(us30pc, calBP<5800 & calBP>399)

pcUSdom <- ggplot(us30pc,aes(x=(MKDE*100), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=factor(PeriodID2)), size=3.5) +
  geom_path(aes(),linewidth=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "A. TCP KDE Per Capita Growth vs. Density")
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
pcUSdom

uscpt <- ggplot(us30pc,aes(x=(calBP), y=(MKDE*100))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=MKDE*100, color=factor(PeriodID2)), size=3.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  #scale_x_reverse()+
  scale_x_reverse(breaks=c(5500, 4500, 3500, 2500,1500, 500), limits=c(5800,400))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"), legend.position = "none")+
  labs(x = "Years cal BP", y="Mean KDE", title = "A. TCP KDE vs. Time")+
geom_vline(xintercept = 1400)+
  geom_vline(xintercept = 2950)+
  geom_vline(xintercept = 4300)+
  geom_vline(xintercept = 1200)+
annotate("text", x =5500, y = 1.25, label = "Period 1", size = 6)+
annotate("text", x =3800, y = 1.25, label = "Period 2", size = 6)+
annotate("text", x =2200, y = 1.25, label = "Period 3", size = 6)+
annotate("text", x =1300, y = 1.25, label = "Period 4", size = 6)+
annotate("text", x =800, y = 1.25, label = "Period 5", size = 6)
uscpt


###Create Stacked Plots of population dynamics and changes in isotopes over time.

Fig5<-plot_grid(uscpt, pcUSdom, ncol=1, align="hv", axis = "rl")
Fig5

pdf("graphics/Test2.pdf", width=15.55, height=17.55)
Fig5
dev.off()

#Riverine ecosystems plot
Fig5<-plot_grid(uscpt, ptime1, ncol=1, align="hv", axis = "rl")
Fig5

pdf("graphics/RiverineSexRev2.pdf", width=15.55, height=17.55)
Fig5
dev.off()

##Coastal ecosystems plot
Figcoast<-plot_grid(uscpt, ptime2, ncol=1, align="hv", axis = "rl")
Figcoast

pdf("graphics/CoastSexRev2.pdf", width=15.55, height=17.55)
Figcoast
dev.off()

####Texas Coastal Plain Radiocarbon Population Expansion Estimates===========================================================
box<- read.csv("data/FinalRCDTexas3.csv")
box2<- subset(box, Region=="TCP")

counties<-map_data("state")

ArchGlobeMap<-ggplot() +
  geom_polygon(data = counties, mapping = aes(x = long, y = lat, group = group),
               fill = "grey", color = "white") +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)+
  #    fill = "white", color="black") +
  theme_bw()+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=20, face = "bold"))+
  labs(x = "Longitude", y="Latitude", title = "US ArchaeoGlobe Regions and Radiocarbon")+
  geom_point(data=box2, aes(Longitude, Latitude, color=factor(Region)),
             inherit.aes = FALSE, alpha = 0.5, size = 2, shape=23)
ArchGlobeMap

#3Calibrate the dates
cptcal <- calibrate(x = box2$Age,  errors = box2$Error, calCurves = "intcal20",  normalised = FALSE)
boxbins <- binPrep(sites = box2$Trinomial, ages = box2$Age, h = 100)

####construct SPD from 8200 to 200 cal BP
spd.CTx <- spd(cptcal, bins=boxbins, runm=200, timeRange=c(8200,200))
plot(spd.CTx, runm=200, xlim=c(8200,200), type="simple")

#save SPD and calBP for later
PrDens<-spd.CTx$grid$PrDens
calBP<-spd.CTx$grid$calBP

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))


####make KDEs for the TCP
US.randates = sampleDates(cptcal, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(8200,200),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "data/TCPKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum
dd2c<- read.csv("data/TCPKDE50bin.csv") %>%
  dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
sums<-cbind(MKDE, out50)
write.table(sums, file = "data/TCPSumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("data/TCPSumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
write.table(pcgrowth, file = "data/TCPPerCap.csv", sep = ",", col.names=NA)


###Test BayesFactor comparison of means. This helps understand the how the Bayseian t-test operates

##Simulate data
set.seed(123)
group1 <- rnorm(30, mean = 50, sd = 10)  # group 1 data
group2 <- rnorm(30, mean = 59, sd = 10)  # group 2 data
test<-cbind(group1,group2)

#Plot histograms for group1 and group2

m <- ggplot(test, aes(x=(group1))) + geom_histogram(binwidth = 1) + geom_vline(xintercept = 50)+
  scale_x_continuous(limits=c(40, 65))

z <- ggplot(test, aes(x=(group2))) + geom_histogram(binwidth = 1)+ geom_vline(xintercept = 55) +
  scale_x_continuous(limits=c(40, 65))

#Stack the histograms
Fig5Rev<-plot_grid(m,z, ncol=1, align="hv", axis = "rl")
Fig5Rev

##Bayesian t-test of simulated data
bayesian_ttest <- ttestBF(x = group1, y = group2)
# Show the result
print(bayesian_ttest)

##Change the rscale parameter of the Bayesian t-test
bayesian_ttest_custom <- ttestBF(x = group1, y = group2, rscale = .707)

# Print the result with custom prior
print(bayesian_ttest_custom)

