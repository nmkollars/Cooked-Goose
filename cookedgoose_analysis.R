##########
#Project title: "Sequential disturbances alter the outcome of intergentoyptic interactions in a clonal plant" 
#Project nickname: "Cooked Goose" 
#Authors: Kollars, N. M.; DuBois, K.; Stachowicz, J. J. 
#Script written by: Nicole M. Kollars (nmkollars@ucdavis.edu)
##########

#Nicole's working directory
setwd("C:/Users/nmkol/OneDrive - UC Davis/Cooked_Goose")

#Needed packages
library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(lme4)
library(AER)
library(cowplot)
library(rcompanion)
library(scales)
 
#Loading data 
dat<-read.csv("cookedgoose_dat.csv")

##########
#Variable descriptions
##########

#pair:	pair id
#trt:	treatment indentification
#temp: temperature treatment; hold or cold
#clip: clipping treatment; clipped or not
#color:	color of tag used to identify genotype; "UK" means that the identity of genotype for ramet was unknown at breakdown
#genotype: the genotype identification according to naming by Abbott et al. 2018 and Hughes et al. 2009
#start.len: length of the longest leaf at initial planting
#start.rhiz.diam: diameter of the rhizome at initial planting
#start.rhiz.len: length of rhizome at initial planting
#no.shoots.pre.heat: # shoots post acclimation period but before heating treatment began
#no.shoots.post.heat:	# shoots after the heating treatment
#length1.T1: length # 1 immediately after heating treatment
#length2.T1: length #2 immediately after heating treatment
#length3.T1: length #3 immediately after heating treatment
#length4.T1: length #4 immediately after heating treatment
#lengthtot.T1: total length grown after heating treatment
#days.T1: total days between punching and length measurements
#growth.rate.T1: total length divided by total # days calculates growth rate immediately after the heating treatment
#dead.T1: identifies whether genotype was dead ("1") or alive ("0") immediately after heating treatment; assessed at the same time as growth measurements
#dead.TR: identifies whether genotype was dead ("1") or alive ("0") at the time of transplant; assessed at the same time as growth measurements
#dead.preC: identifies whether genotype was dead ("1") or alive ("0") before implementing the clipping treatment; assessed at the same time as growth measurements
#no.shoots.preC1: number of shoots at the plot level before the first clipping; note that if we cannot distinguish # shoots per genotype, the shoots are listed under the "blue" genotype. If we can distinguish the genotypes by knowing which genotypes had previously died, the # of shoots is listed next to the living genotype and the # of shoots for the dead genotype are indicated by 0's
#no.shoots.postC1: number of shoots at the plot level after the first clipping; note that if we cannot distinguish # shoots per genotype, the shoots are listed under the "blue" genotype. If we can distinguish the genotypes by knowing which genotypes had previously died, the # of shoots is listed next to the living genotype and the # of shoots for the dead genotype are indicated by 0's
#dead.postC: identifies whether genotype was dead ("1") or alive ("0") post clipping treatment; assessed at the same time as breakdown growth measurements
#status.bd: the status of the genotype at the breakdown; "NA" means wasn't measured, "dead" means it either died before breakdown or after the 2nd clipping, "measured" means the shoot was identified and measured during breakdown. Note an assumption: if there is only one identifiable genet in the pot at breakdown, assume that the other genet died and isn't missing (aka "NA"), even if there was no evidence of death prior to clipping 
#no.shoots.bd:	#. shoots at the breakdown, "0" if genet is dead but other genet still present, "NA" if both genets are dead
#width.bd: sheath width of terminal shoot (or largest shoot) at breakdown
#date.measured.bd: date that the leaf extensions on the shoots were measured
#days.bd:	# of days between punching and measurement at breakdown
#length1.bd: length #1 at breakdown
#length2.bd: length #2 at breakdown
#length3.bd: length #3 at breakdown
#length4.bd: length #4 at breakdown
#length5.bd: length #5 at breakdown
#length6.bd: length #6 at breakdown
#length7.bd: length #7 at breakdown
#lengthtot.bd: total length grown at breakdown
#growthrate.bd: total length divided by total # days calculates growth rate at breakdwon
#grazing: evidence for amphipod grazing on the leaves at breakdown
#laby: evidence for the presence of labyrinthyla on the leaves at breakdown
#ss.np:	# of side shoots that were not punched to be assessed for new growth at the breakdown; these shoots would not have been included in the "new growth" biomass
#dessicated: evidence that the leaves were dessicated or dead at breakdown; this likely resulted from being emmerged in the sun too long during punching
#flowering: # of flowering shoots during breakdown
#gen.id: id number for the tissue sample set aside for genetic analysis if needed
#old	biomass: of the "old growth" tissue for the entire genotype
#new	biomass: of the "new growth" tissue for the entire genotype
#rhiz: biomass of the rhizome for the entire genotype
#roots:	biomass of the roots for the entire genotype
#rhizbrokeT1:	this indicates whether the rhizome broke ("1") or not ("0") while measuring growth at T1 
#rhizbroke.TR: this indicates whether the rhizome broke ("1") or not ("0") during transplanting
#squashed: this indicates whether the shoot was squashed ("1") or not ("0"), likely near the tag, prior to transplant
#SS.lost: this indicates whether a side shoot was lost ("1") or not ("0") at any point during the experiment; some note of this is sometimes repeated (sometimes not) in the notes
#replicate:	"1" for yes, "0" for no - no's include plants died before trt application, unknowns, contradictory, accidentally clipped samples, "2" for pots dead at breakdown, "3"  for only reps that had all four treatments present at breakdown
#notes: any observation on the individual genotype that may be useful in analysis 

##########
#preparing the data 
##########

#removes pots that are not considered replicates 
#see supplementary Table 1 for reasons for exclusion 
dat<-dat[dat$replicate>"0",]

#removing fill-in NAs at the end of the dataset 
dat<-dat[1:310,]

#removes replicates where both shoots died by breakdown
alive<-dat[dat$replicate!="2",]

###########
#Aggregate (pot-level) density and biomass
###########

#DENSITY (shoot counts)
dat1 <- alive %>%
  select_at(vars(pair,trt, temp,clip, color, no.shoots.pre.heat, no.shoots.post.heat,
                 no.shoots.preC1,no.shoots.postC1, no.shoots.bd)) 
dat1$no.shoots.pre.heat<-as.numeric(dat1$no.shoots.pre.heat)

#The NAs for "R" for shoot counts for preC1 and post C1 are just placeholders 
#because we only counted shoots at the pot level. #So I will convert these NAs to zeros,
dat1$no.shoots.preC1[is.na(dat1$no.shoots.preC1)]<-0
dat1$no.shoots.postC1[is.na(dat1$no.shoots.postC1)]<-0

#create B- and R- specific shoot count variables for each time step
#no.shoots.pre.heat
dat1.B.sc.in <- dat1 %>%
  filter(color == "B") %>%
  rename(no.shoots.B.in = no.shoots.pre.heat) %>%
  select_at(vars(pair,temp,clip, trt, no.shoots.B.in))
dat1.R.sc.in <- dat1 %>%
  filter(color == "R") %>%
  rename(no.shoots.R.in = no.shoots.pre.heat)%>%
  select_at(vars(pair,temp,clip,trt, no.shoots.R.in))
#no.shoots.post.heat
dat1.B.sc.t1 <- dat1 %>%
  filter(color == "B") %>%
  rename(no.shoots.B.t1 = no.shoots.post.heat) %>%
  select_at(vars(pair,temp,clip, trt, no.shoots.B.t1))
dat1.R.sc.t1 <- dat1 %>%
  filter(color == "R") %>%
  rename(no.shoots.R.t1 = no.shoots.post.heat)%>%
  select_at(vars(pair,temp,clip, trt, no.shoots.R.t1))
#no.shoots.preC1
dat1.B.sc.tr <- dat1 %>%
  filter(color == "B") %>%
  rename(no.shoots.B.tr = no.shoots.preC1) %>%
  select_at(vars(pair, trt, temp, clip, no.shoots.B.tr))
dat1.R.sc.tr <- dat1 %>%
  filter(color == "R") %>%
  rename(no.shoots.R.tr = no.shoots.preC1)%>%
  select_at(vars(pair, trt,temp,clip, no.shoots.R.tr))
#no.shoots.postC1
dat1.B.sc.pc <- dat1 %>%
  filter(color == "B") %>%
  rename(no.shoots.B.pc = no.shoots.postC1) %>%
  select_at(vars(pair, trt,temp,clip, no.shoots.B.pc))
dat1.R.sc.pc <- dat1 %>%
  filter(color == "R") %>%
  rename(no.shoots.R.pc = no.shoots.postC1)%>%
  select_at(vars(pair, trt,temp,clip, no.shoots.R.pc))
#no.shoots.bd
dat1.B.sc.bd <- dat1 %>%
  filter(color == "B") %>%
  rename(no.shoots.B.bd = no.shoots.bd) %>%
  select_at(vars(pair, trt,temp,clip, no.shoots.B.bd))
dat1.R.sc.bd <- dat1 %>%
  filter(color == "R") %>%
  rename(no.shoots.R.bd = no.shoots.bd)%>%
  select_at(vars(pair, trt,temp,clip, no.shoots.R.bd))
#combine all the new variables into one datasheet - there is probably
#an easier way to do this..but brute force won the day
dat1.sc.in<-full_join(dat1.B.sc.in, dat1.R.sc.in)
dat1.sc.t1<-full_join(dat1.B.sc.t1, dat1.R.sc.t1)
dat1.sc.tr<-full_join(dat1.B.sc.tr, dat1.R.sc.tr)
dat1.sc.pc<-full_join(dat1.B.sc.pc, dat1.R.sc.pc)
dat1.sc.bd<-full_join(dat1.B.sc.bd, dat1.R.sc.bd)
dat1.sc.int1<-full_join(dat1.sc.in, dat1.sc.t1)
dat1.sc.int1tr<-full_join(dat1.sc.int1,dat1.sc.tr)
dat1.sc.int1trpc<-full_join(dat1.sc.int1tr,dat1.sc.pc)
dat1.sc<-full_join(dat1.sc.int1trpc,dat1.sc.bd)

#create variables that will give the total # of shoots in the pot 
dat2 <- transmute(dat1.sc, pair = pair,temp=temp,clip=clip, trt = trt,
                  shoots.in = (no.shoots.B.in + no.shoots.R.in),
                  shoots.t1 = (no.shoots.B.t1 + no.shoots.R.t1),
                  shoots.tr = (no.shoots.B.tr + no.shoots.R.tr),
                  shoots.pc = (no.shoots.B.pc + no.shoots.R.pc),
                  shoots.bd = (no.shoots.B.bd + no.shoots.R.bd))
names(dat2)[names(dat2)=="trt"]  <- "Treatment"
dat2$Treatment <- factor(dat2$Treatment, level = c("A", "B", "C", "D"), labels = c("Control (35)", "Warmed (37)","Clipped (38)","Warmed + Clipped (33)"))
dat3<-melt(dat2,id.vars=c(1:4))

#stats
dat3$temp<-ifelse(dat2$temp == "cold", 0, ifelse(dat2$temp == "hot", 1,NA))
dat3$clip<-ifelse(dat2$clip == "none", 0, ifelse(dat2$clip == "clip", 1,NA))
#no shoots when heat turned off 
shoots.t1<-dat3[dat3$variable=="shoots.t1",]
hist(shoots.t1$value)
fit.shoots<-glm(value~temp*clip,data=shoots.t1,family=poisson)
dispersiontest(fit.shoots, trafo = 1)
fit.shoots<-glm(value~temp,data=shoots.t1,family=poisson)
summary(fit.shoots)
hist(resid(fit.shoots))
par(mfrow=c(2,2))
plot(fit.shoots)
par(mfrow=c(1,1))
anova(fit.shoots,test="Chi")
#no shoots 10 weeks post heating
shoots.tr<-dat3[dat3$variable=="shoots.tr",]
hist(shoots.tr$value)
fit.shoots<-glm(value~temp,data=shoots.tr,family=poisson)
dispersiontest(fit.shoots, trafo = 1)
fit.shoots<-glm(value~temp,data=shoots.tr,family=quasipoisson)
summary(fit.shoots)
hist(resid(fit.shoots))
par(mfrow=c(2,2))
plot(fit.shoots)
par(mfrow=c(1,1))
anova(fit.shoots,test="F")
#no shoots post- first clipping
shoots.pc<-dat3[dat3$variable=="shoots.pc",]
hist(shoots.pc$value)
fit.shoots<-glm(value~temp*clip,data=shoots.pc,family=poisson)
dispersiontest(fit.shoots, trafo = 1)
fit.shoots<-glm(value~temp*clip,data=shoots.pc,family=quasipoisson)
summary(fit.shoots)
hist(resid(fit.shoots))
par(mfrow=c(2,2))
plot(fit.shoots)
par(mfrow=c(1,1))
anova(fit.shoots,test="F")
#no shoots at breakdown 
shoots.bd<-dat3[dat3$variable=="shoots.bd",]
hist(shoots.bd$value)
fit.shoots<-glm(value~temp*clip,data=shoots.bd,family=poisson)
dispersiontest(fit.shoots, trafo = 1)
fit.shoots<-glm(value~temp*clip,data=shoots.bd,family=quasipoisson)
summary(fit.shoots)
hist(resid(fit.shoots))
par(mfrow=c(2,2))
plot(fit.shoots)
par(mfrow=c(1,1))
anova(fit.shoots,test="F")

#BIOMASS 
dat4 <- transmute(alive, pair = pair, trt = trt, temp=temp, clip=clip,color=color,tank=tank, 
                  above = old + new,
                  below = rhiz + roots) 

#create B- and R- specific shoot count variables for each time step
#above
dat4.B.above<- dat4 %>%
  filter(color == "B") %>%
  rename(above.B = above) %>%
  select_at(vars(pair, trt, tank,temp,clip, above.B))
dat4.R.above <- dat4 %>%
  filter(color == "R") %>%
  rename(above.R = above)%>%
  select_at(vars(pair, trt, tank,temp,clip, above.R))
dat4.above<-full_join(dat4.B.above, dat4.R.above)
#below
dat4.B.below<- dat4 %>%
  filter(color == "B") %>%
  rename(below.B = below) %>%
  select_at(vars(pair, trt, tank,temp,clip, below.B))
dat4.R.below <- dat4 %>%
  filter(color == "R") %>%
  rename(below.R = below)%>%
  select_at(vars(pair, trt,tank,temp,clip, below.R))
dat4.below<-full_join(dat4.B.below, dat4.R.below)
#combine all the new variables into one datasheet 
dat4.abbe<-full_join(dat4.above, dat4.below)

#create variables that will give the total biomass in the pot 
dat5 <- transmute(dat4.abbe, pair = pair,tank=tank,temp=temp,clip=clip, trt = trt,
                  above.pot = above.B + above.R,
                  below.pot = below.B + below.R)
dat2$clip <- factor(dat2$clip, levels = c("none", "clip"))

#stats
dat5$temp<-ifelse(dat5$temp == "cold", 0, ifelse(dat5$temp == "hot", 1,NA))
dat5$clip<-ifelse(dat5$clip == "none", 0, ifelse(dat5$clip == "clip", 1,NA))
#above
par(mfrow=c(1,1))
hist(dat5$above.pot)
fit.above<-glm(above.pot~temp*clip,data=dat5,family=Gamma(link="inverse"))
summary(fit.above)
hist(resid(fit.above))
par(mfrow=c(2,2))
plot(fit.above)
par(mfrow=c(1,1))
print(anova(fit.above,test="F"))
#below
fit.below<-glm(below.pot~temp*clip,data=dat5,family=Gamma(link="inverse"))
summary(fit.below)
hist(resid(fit.below))
par(mfrow=c(2,2))
plot(fit.below)
par(mfrow=c(1,1))
print(anova(fit.below,test="F"))

#Figure 1 

#setting position
pd <- position_dodge(0.1) 

#code for making minor tick marks 
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

#Figure 1A 
sumshoots<- ddply(dat3, c("Treatment","variable"), summarise, 
                  N=length(value),
                  mean=mean(value,na.rm=TRUE),
                  sd=sd(value,na.rm=TRUE),
                  se=sd/sqrt(N))
sumshoots$variable<-revalue(sumshoots$variable, c("shoots.in"="10/2/2017", "shoots.t1"="11/11/2017",
                      "shoots.tr"="01/26/2018","shoots.pc"="03/9/2018","shoots.bd"="04/23/2018"))
sumshoots$variable<-as.Date(sumshoots$variable,"%m/%d/%Y")
test<-as.Date(c("10/3/2017","11/12/2017","01/26/2018","03/09/2018"),"%m/%d/%Y")
test2<-as.Date(c("9/25/2017","11/03/2017","01/17/2018","02/28/2018"),"%m/%d/%Y")
custom_breaks <- seq(0, 24, 2)
shoot.plot<-ggplot(sumshoots,aes(variable,y=mean,shape=Treatment))+
  geom_vline(xintercept=test,col="gray30",lty=2,lwd=.5)+
  geom_line(position=pd,lwd=0.8)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=20, position=pd)+
  geom_point(position=pd, size=3)+
  scale_shape_manual(values=c(21,19,24,17))+
  annotate("text",x=test2,y=c(12,12,12,14),
           label=c("Heat on","Heat off","1st clipping","2nd clipping"), 
           angle=90,size=3)+
  ylab("No. Shoots per pot") +xlab(NULL)+
  scale_x_date(labels = date_format("%b-%y"),date_breaks = "1 month",
               limits=c(as.Date("09/15/2017","%m/%d/%Y"),as.Date("05/01/2018","%m/%d/%Y")))+
  scale_y_continuous(limits=c(2,24),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,2, inverse = TRUE))+
  theme_classic()+ 
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10),angle=60, hjust=1),
        legend.position=c(0.3,0.8),
        legend.title=element_blank(),
        legend.text=element_text(size=8))

#Figure 1B 
#above
dat6<-dat5 %>% drop_na(above.pot)
sumabove<- ddply(dat6, c("temp","clip"), summarise, 
                 N=length(above.pot),
                 mean=mean(above.pot),
                 sd=sd(above.pot),
                 se=sd/sqrt(N))
#below
dat7<-dat5 %>% drop_na(below.pot)
sumbelow<- ddply(dat7, c("temp","clip"), summarise, 
                 N=length(below.pot),
                 mean=mean(below.pot),
                 sd=sd(below.pot),
                 se=sd/sqrt(N))
sumbelow$mean<--sumbelow$mean
test<-rbind(sumabove,sumbelow)
test$type<-c("above","above","above","above","below","below",
             "below","below")
test$trt<-c("Control","Clipped","Warmed","Warmed + Clipped",
            "Control","Clipped","Warmed","Warmed + Clipped")
test$trt <- factor(test$trt, levels = c("Control","Warmed","Clipped","Warmed + Clipped"))
nabove<-test[1:4,1:3]
nbelow<-test[5:8,1:3]
custom_breaks <- seq(-4,5,1)
biomass.plot<-ggplot(test,aes(x=as.factor(clip),y=mean,fill=as.factor(temp)))+
  geom_bar(stat="identity", position = "dodge")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),  width=0.3,position=position_dodge(0.9))+
  scale_fill_manual(NULL,labels = c("Ambient", "Warmed"),values = c("gray","gray40")) +
  theme_classic()+ylab("Biomass (g dry mass) per pot")+ xlab(NULL)+
  scale_x_discrete(labels=c("Not clipped", "Clipped"))+
  geom_hline(yintercept=0)+
  scale_y_continuous(limits=c(-4,5),breaks = custom_breaks,
                     labels = every_nth(c(4,3,2,1,0,1,2,3,4,5),1, inverse = TRUE))+
  geom_text(data=nabove,aes(y=0.2,label=N),position=position_dodge(0.9),color="white",size=2.5)+
  geom_text(data=nbelow,aes(y=-0.2,label=N),position=position_dodge(0.9),color="white",size=2.5)+
  annotate("text", x = 2, y = 2.5, label = "Aboveground",size=3.5) +
  annotate("text", x = 2, y = -2.5,label = "Belowground",size=3.5)+
  theme(axis.ticks.length=unit(-0.25, "cm"),
      axis.text=element_text(size=8),
      axis.title.y = element_text(margin = margin(r = 6),size=10),
      axis.text.y=element_text(margin=margin(r=10)),
      axis.text.x=element_text(margin=margin(t=10)),
      legend.position=c(0.8,0.9), 
      legend.text=element_text(size=8))

#combining the two plots together
plot_grid(shoot.plot, biomass.plot,ncol=1, 
          align = "v", axis = "l",labels=c("A)","B)"),
          label_size=10)
ggsave("FE-2020-00384_Fig1.png",height=7,width=3.5,dpi=300)

ggsave("FE-2020-00384_Fig1.pdf",height=7,width=3.5,dpi=300)

##########
#Interaction outcome between genotypes 
##########

#isolate belowground biomass + potential covariates 
dat8 <- transmute(alive, pair = pair,tank=tank, trt = trt, temp=temp,clip=clip,color=color, 
                  start.len=start.len,
                  below=rhiz+roots,
                  diam=start.rhiz.diam/10, #divide by 10 to convert from cm to mm
                  rad=diam/2,
                  vol=pi*(rad^2)*start.rhiz.len)
#create B- and R- specific columns for each variable 
dat8.B.below <- dat8 %>%
  filter(color == "B") %>%
  rename(below.B = below) %>%
  rename(vol.B=vol)%>%
  rename(start.len.B=start.len)%>%
  select_at(vars(pair,trt,temp,clip,tank,below.B,vol.B,start.len.B))
dat8.R.below<- dat8 %>%
  filter(color == "R") %>%
  rename(below.R = below)%>%
  rename(vol.R = vol)%>%
  rename(start.len.R = start.len)%>%
  select_at(vars(pair,trt,temp,clip,tank,below.R,vol.R,start.len.R))
#combine all the new variables into one datasheet - there is probably
#an easier way to do this..but brute force won the day
dat8.below<-full_join(dat8.B.below, dat8.R.below)
dat9 <- transmute(dat8.below, pair = pair, tank=tank, trt = trt,temp=temp,clip=clip,
                  start.len=start.len.B + start.len.R,
                  below=below.B + below.R,
                  b.rel=round(below.B/(below.B + below.R),2),
                  r.rel=round(below.R/(below.B + below.R),2),
                  below.raw=below.B-below.R,
                  vol.tot=vol.R + vol.B,
                  start.diff=abs(start.len.R-start.len.B)/(start.len.R+start.len.B)*100,
                  vol.diff=abs(vol.R-vol.B)/(vol.tot)*100,
                  below.diff=abs(below.B-below.R)/(below.R+below.B)*100)
dat9$clip <- factor(dat9$clip, levels = c("none", "clip"))
dat9$trt<-factor(dat9$trt,labels=c("Control","Warmed","Clipped","Warmed + Clipped"))

#graphing potential covariates
#does initial diff in shoot length between genets predict final diff in biomass? NO
shoot.size<-ggplot(dat9,aes(x=start.diff,y=below.diff))+
  geom_point(cex=3)+
  theme_bw()+
  ylim(c(0,100))+
  xlim(c(0,80))+
  xlab("Initial % difference shoot length")+
  ylab("Final % difference belowground biomass")+
  facet_wrap(~trt)+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=10),
        axis.title.y = element_text(margin = margin(r = 6),size=12),
        axis.title.x = element_text(margin = margin(t = 6),size=12),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)))
#does initial diff in rhiz volume between genets predict final diff in biomass? NO
rhiz.vol<-ggplot(dat9,aes(x=vol.diff,y=below.diff))+
  geom_point(cex=3)+
  theme_bw()+
  ylim(c(0,100))+
  xlim(c(0,80))+
  xlab("Initial % difference rhizome volume")+
  ylab("Final % difference belowground biomass")+
  facet_wrap(~trt)+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=10),
        axis.title.y = element_text(margin = margin(r = 6),size=12),
        axis.title.x = element_text(margin = margin(t = 6),size=12),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)))

#combining the two panels together into a single figure S2
plot_grid(shoot.size, rhiz.vol,ncol=1, 
            align = "v", axis = "l",labels=c("A)","B)"),
            label_size=12)
ggsave("Fig_S2.png",height=8,width=6,dpi=300)


#isolate genets that have relative abundance of biomass > 0.5 in control
#to compare to what that genet's relative abundance of biomass is in the other 
#treatments 
dat10 <- select_at(dat9,vars(pair,trt,temp,clip,b.rel,r.rel))

#create a subset of your data that only corresponds to the controls
control<- dat10[dat10$trt=='Control',] 
#determine which genotype was the winner for each pair
control$winner <- ifelse(control$b.rel > control$r.rel, 'B', 'R') 
#append this back to the original dataframe
dat11 <- merge(dat10, control[,c(1,7)], by = 'pair') 
#calculate the difference between genotypes
dat11$max <- ifelse(dat11$winner == 'R', dat11$r.rel, dat11$b.rel)
dat11$clip <- factor(dat11$clip, levels = c("none","clip"))
#stats - note: permutation analyses result in different stats each time run
numperm=1000
print(anova(lm(max ~ trt, data = dat11)))
retval <- matrix(0, ncol = 1, nrow = numperm + 1)
retval[1, ] <- anova(lm(max ~ trt, data = dat11))[1, 4]
for (i in 1:numperm) {
  dat11$perm <- sample(dat11$trt, length(dat11$trt), replace = F)
  retval[i + 1, ] <- anova(lm(max ~ perm, data = dat11))[1, 4]
}
par(mfrow = c(1, 1))
probALL = NULL
for (i in 1) {
  prob <- sum(retval[-1, i] >= retval[1, i])/(dim(retval)[1] - 1)
  probALL <- rbind(probALL, prob)
  hist(retval[-1, i], main = paste("Observed F-statistic =", round(retval[1, i],2),";", "p =",prob),
       xlab = "Simulated F-statistic")
  points(c(retval[1, i], retval[1, i]), c(0, numperm), type = "l", lwd = 2, 
         col = "gray",lty=2)
}
PT<-pairwisePermutationTest(max~trt,data=dat11, method = "fdr")
cldList(comparison = PT$Comparison,
        p.value    = PT$p.adjust,
        threshold  = 0.05)
#figure 
anno<-data.frame(Treatment=c("Warmed","Warmed","Warmed","Warmed + Clipped", "Warmed + Clipped", "Warmed + Clipped"),
                 xstar = c(1,2,3,1,2,3),
                 ystar = c(0.45,0.45,0.45,0.45,0.45,0.45),
                 lab = c("a","a","b","ab","a","b")) 
custom_breaks<-seq(0,1,0.125)
ggplot(dat11,aes(x=clip,y=max,fill=temp))+
  geom_boxplot(outlier.colour=NA)+
  geom_dotplot(binaxis="y",binwidth=0.02,stackdir="center",position=position_dodge(0.8))+
  scale_fill_manual(NULL,labels = c("Ambient", "Warmed"),values = c("gray","gray40")) +
  theme_classic()+
  geom_hline(yintercept=0.5,lty=2,lwd=0.8)+
  ylab("Relative abundance (belowground biomass)")+
  xlab(NULL)+
  annotate(geom="text", x=c(0.8,1.2,1.8,2.2), y=c(1.1,1.1,1.1,1.1), label=c("A","AB","AB","B"),
           color="black",size=3)+
  scale_x_discrete(labels=c("Not Clipped","Clipped"))+
  scale_y_continuous(limits=c(0,1.1),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,2, inverse = TRUE))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
       axis.text=element_text(size=8),
       axis.title.y = element_text(margin = margin(r = 6),size=10),
       axis.text.y=element_text(margin=margin(r=10)),
       axis.text.x=element_text(margin=margin(t=10)),
       legend.text=element_text(size=8),
       legend.position=c(0.19,0.2))

ggsave("FE-2020-00384_Fig2.png", dpi=300, units="in", width=3.5, height=3.5)
ggsave("FE-2020-00384_Fig2.pdf", dpi=300, units="in", width=3.5, height=3.5)

#now do the same for each of the remaining treatments 
#Warming only
cooked<- dat10[dat10$trt=='Warmed',] 
cooked$winner <- ifelse(cooked$b.rel > cooked$r.rel, 'B', 'R') 
dat12 <- merge(dat10, cooked[,c(1,7)], by = 'pair')
dat12$max <- ifelse(dat12$winner == 'R', dat12$r.rel, dat12$b.rel) 
dat12$clip <- factor(dat12$clip, levels = c("none","clip"))
#stats
numperm=1000
print(anova(lm(max ~ trt, data = dat12)))
retval <- matrix(0, ncol = 1, nrow = numperm + 1)
retval[1, ] <- anova(lm(max ~ trt, data = dat12))[1, 4]
for (i in 1:numperm) {
  dat12$perm <- sample(dat12$trt, length(dat12$trt), replace = F)
  retval[i + 1, ] <- anova(lm(max ~ perm, data = dat12))[1, 4]
}
par(mfrow = c(1, 1))
probALL = NULL
for (i in 1) {
  prob <- sum(retval[-1, i] >= retval[1, i])/(dim(retval)[1] - 1)
  probALL <- rbind(probALL, prob)
  hist(retval[-1, i], main = paste("Observed F-statistic =", round(retval[1, i],2),";", "p =",prob),
       xlab = "Simulated F-statistic")
  points(c(retval[1, i], retval[1, i]), c(0, numperm), type = "l", lwd = 2, 
         col = "gray",lty=2)
}
PT<-pairwisePermutationTest(max~trt,data=dat12, method = "fdr")
cldList(comparison = PT$Comparison,
        p.value    = PT$p.adjust,
        threshold  = 0.05)
#figure
cook.plot<-ggplot(dat12,aes(x=clip,y=max,fill=temp))+
  geom_boxplot(outlier.colour=NA)+
  geom_dotplot(binaxis="y",binwidth=0.02,stackdir="center",position=position_dodge(0.8))+
  scale_fill_manual(NULL,labels = c("Ambient", "Warmed"),values = c("gray","gray40")) +
  theme_classic()+
  geom_hline(yintercept=0.5,lty=2,lwd=0.8)+
  ylab("Relative abundance (belowground biomass)")+
  xlab(NULL)+
  annotate(geom="text", x=c(0.8,1.2,1.8,2.2), y=c(1.1,1.1,1.1,1.1), label=c("A","B","A","A"),
           color="black",size=3)+
  scale_x_discrete(labels=c("Not Clipped","Clipped"))+
  scale_y_continuous(limits=c(0,1.1),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,2, inverse = TRUE))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position="none")

#clipping only
goose<- dat10[dat10$trt=='Clipped',] 
goose$winner <- ifelse(goose$b.rel > goose$r.rel, 'B', 'R') 
dat13 <- merge(dat10, goose[,c(1,7)], by = 'pair') 
dat13$max <- ifelse(dat13$winner == 'R', dat13$r.rel, dat13$b.rel) 
dat13$clip <- factor(dat13$clip, levels = c("none","clip"))
#stats
numperm=1000
print(anova(lm(max ~ trt, data = dat13)))
retval <- matrix(0, ncol = 1, nrow = numperm + 1)
retval[1, ] <- anova(lm(max ~ trt, data = dat13))[1, 4]
for (i in 1:numperm) {
  dat13$perm <- sample(dat13$trt, length(dat13$trt), replace = F)
  retval[i + 1, ] <- anova(lm(max ~ perm, data = dat13))[1, 4]
}
par(mfrow = c(1, 1))
probALL = NULL
for (i in 1) {
  prob <- sum(retval[-1, i] >= retval[1, i])/(dim(retval)[1] - 1)
  probALL <- rbind(probALL, prob)
  hist(retval[-1, i], main = paste("Observed F-statistic =", round(retval[1, i],2),";", "p =",prob),
       xlab = "Simulated F-statistic")
  points(c(retval[1, i], retval[1, i]), c(0, numperm), type = "l", lwd = 2, 
         col = "gray",lty=2)
}
PT<-pairwisePermutationTest(max~trt,data=dat13, method = "fdr")
cldList(comparison = PT$Comparison,
        p.value    = PT$p.adjust,
        threshold  = 0.05)
#figure
goose.plot<-ggplot(dat13,aes(x=clip,y=max,fill=temp))+
  geom_boxplot(outlier.colour=NA)+
  geom_dotplot(binaxis="y",binwidth=0.02,stackdir="center",position=position_dodge(0.8))+
  annotate(geom="text", x=c(0.8,1.2,1.8,2.2), y=c(1.1,1.1,1.1,1.1), label=c("A","A","B","A"),
           color="black",size=3)+
  scale_fill_manual(NULL,labels = c("Ambient", "Warmed"),values = c("gray","gray40")) +
  theme_classic()+
  geom_hline(yintercept=0.5,lty=2,lwd=0.8)+
  ylab("Relative abundance (belowground biomass)")+ 
  xlab(NULL)+
  scale_x_discrete(labels=c("Not Clipped","Clipped"))+
  scale_y_continuous(limits=c(0,1.1),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,2, inverse = TRUE))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position="none")

#Warming+Clipping
cookgoose<- dat10[dat10$trt=='Warmed + Clipped',]
cookgoose$winner <- ifelse(cookgoose$b.rel > cookgoose$r.rel, 'B', 'R') 
dat14 <- merge(dat10, cookgoose[,c(1,7)], by = 'pair') 
dat14$max <- ifelse(dat14$winner == 'R', dat14$r.rel, dat14$b.rel)
#stats
numperm=1000
print(anova(lm(max ~ trt, data = dat14)))
retval <- matrix(0, ncol = 1, nrow = numperm + 1)
retval[1, ] <- anova(lm(max ~ trt, data = dat14))[1, 4]
for (i in 1:numperm) {
  dat14$perm <- sample(dat14$trt, length(dat14$trt), replace = F)
  retval[i + 1, ] <- anova(lm(max ~ perm, data = dat14))[1, 4]
}
par(mfrow = c(1, 1))
probALL = NULL
for (i in 1) {
  prob <- sum(retval[-1, i] >= retval[1, i])/(dim(retval)[1] - 1)
  probALL <- rbind(probALL, prob)
  hist(retval[-1, i], main = paste("Observed F-statistic =", round(retval[1, i],2),";", "p =",prob),
       xlab = "Simulated F-statistic")
  points(c(retval[1, i], retval[1, i]), c(0, numperm), type = "l", lwd = 2, 
         col = "gray",lty=2)
}
PT<-pairwisePermutationTest(max~trt,data=dat14, method = "fdr")
cldList(comparison = PT$Comparison,
        p.value    = PT$p.adjust,
        threshold  = 0.05)
#figure
dat14$clip <- factor(dat14$clip, levels = c("none","clip"))
cookgoose.plot<-ggplot(dat14,aes(x=clip,y=max,fill=temp))+
  geom_boxplot(outlier.colour=NA)+
  geom_dotplot(binaxis="y",binwidth=0.02,stackdir="center",position=position_dodge(0.8))+
  annotate(geom="text", x=c(0.8,1.2,1.8,2.2), y=c(1.1,1.1,1.1,1.1), label=c("A","A","A","B"),
           color="black",size=3)+
  scale_fill_manual(NULL,labels = c("Ambient", "Warmed"),values = c("gray","gray40")) +
  theme_classic()+
  geom_hline(yintercept=0.5,lty=2,lwd=0.8)+
  ylab("Relative abundance (belowground biomass)")+ 
  xlab(NULL)+
  scale_x_discrete(labels=c("Not Clipped","Clipped"))+
  scale_y_continuous(limits=c(0,1.1),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,2, inverse = TRUE))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position="right")

#combining the panels into a single figure S2 
ggdraw() +
  draw_plot(cook.plot, x = 0.03, y = .5, width = .45, height = .5) +
  draw_plot(goose.plot, x = .48, y = .5, width = .45, height = .5) +
  draw_plot(cookgoose.plot, x = 0.03, y = -0.01, width = 0.625, height = 0.5) +
  draw_plot_label(label = c("A)", "B)", "C)"), size = 10,
                  x = c(0, 0.45, 0), y = c(1, 1, 0.5))    
ggsave("Fig_S3.png",dpi=300,units="in",width=6,height=6)

###################
#integrating trait values
###################

#dataset with trait data 
related<- read.csv("cookedgoose_traits_dat.csv")

# set of traits to be included in trait distance metric. This set is all the traits that differed significantly among genotypes, but you can add or remove them as you like. They are just the column headings
items<-c("Geno", "Pam_a","Pam_etrmax","ratio")
geno<-related$Geno

#Extracting that set of traits from the main file
related2=related[items]
names(related2) <- c("genotype", "a","etr","ratio")
pairs(related2[,2:4],lower.panel=NULL)

#PCA analysis
related2.pc<-prcomp(related2[,2:4],scale=T,center=T)
genos<-related2.pc$x
dat17<-data.frame("genotype"=as.character(related2$genotype),genos[,c(1:3)])
summary(related2.pc)


#isolate needed variables from main dataset
dat18 <- transmute(alive, pair = pair, trt = trt, temp=temp,clip=clip, color=color, 
                  genotype=genotype,below=rhiz+roots)
dat19<-full_join(dat18, dat17)
#create B- and R- specific columns for each variable 
dat19.B.below <- dat19 %>%
  filter(color == "B") %>%
  rename(below.B = below) %>%
  select_at(vars(pair,trt,temp,clip,below.B))
dat19.R.below<- dat19 %>%
  filter(color == "R") %>%
  rename(below.R = below)%>%
  select_at(vars(pair,trt,temp,clip,below.R))
#PC1
dat19.B.pc <- dat19 %>%
  filter(color == "B") %>%
  rename(pc.B = PC1) %>%
  select_at(vars(pair,trt,temp,clip,pc.B))
dat19.R.pc<- dat19 %>%
  filter(color == "R") %>%
  rename(pc.R = PC1)%>%
  select_at(vars(pair,trt,temp,clip,pc.R))
#PC2
dat19.B.pc2 <- dat19 %>%
  filter(color == "B") %>%
  rename(pc2.B = PC2) %>%
  select_at(vars(pair,trt,temp,clip,pc2.B))
dat19.R.pc2<- dat19 %>%
  filter(color == "R") %>%
  rename(pc2.R = PC2)%>%
  select_at(vars(pair,trt,temp,clip,pc2.R))
#PC3
dat19.B.pc3 <- dat19 %>%
  filter(color == "B") %>%
  rename(pc3.B = PC3) %>%
  select_at(vars(pair,trt,temp,clip,pc3.B))
dat19.R.pc3<- dat19 %>%
  filter(color == "R") %>%
  rename(pc3.R = PC3)%>%
  select_at(vars(pair,trt,temp,clip,pc3.R))
#combine all the new variables into one datasheet - there is probably
#an easier way to do this..but brute force won the day
dat19.below<-full_join(dat19.B.below, dat19.R.below)
dat19.pc<-full_join(dat19.B.pc,dat19.R.pc)
dat19.pc2<-full_join(dat19.B.pc2,dat19.R.pc2)
dat19.pc3<-full_join(dat19.B.pc3,dat19.R.pc3)
dat19.bpc<-full_join(dat19.below,dat19.pc)
dat19.bpc2<-full_join(dat19.bpc,dat19.pc2)
dat19.full<-full_join(dat19.bpc2,dat19.pc3)
dat20 <- transmute(dat19.full, pair = pair,trt = trt,temp=temp,clip=clip,
                  b.rel=round(below.B/(below.B + below.R),2),
                  pc.B=pc.B, 
                  pc2.B=pc2.B,
                  pc3.B=pc3.B,
                  bweight=b.rel*pc.B,
                  r.rel=round(below.R/(below.B + below.R),2),
                  pc.R=pc.R,
                  pc2.R=pc2.R,
                  pc3.R=pc3.R,
                  rweight=r.rel*pc.R,
                  trait=bweight+rweight)

#looking at trait value of genotype with relative abundance > 0.5

#genets in control
control<- dat20[dat20$trt=='A',] 
control$winner <- ifelse(control$b.rel > control$r.rel, 'B', 'R') 
control$max <- ifelse(control$winner == 'R', control$r.rel, control$b.rel) 
control$pc <- ifelse(control$winner == 'R', control$pc.R, control$pc.B)
control$pc2<- ifelse(control$winner == 'R', control$pc2.R, control$pc2.B)
control$pc3<- ifelse(control$winner == 'R', control$pc3.R, control$pc3.B)
control<-select_at(control,vars(pair,trt,temp,clip,winner,max,pc,pc2,pc3))
control<-na.omit(control)
#genets in warming only
cooked<- dat20[dat20$trt=='B',] 
cooked$winner <- ifelse(cooked$b.rel > cooked$r.rel, 'B', 'R') 
cooked$max <- ifelse(cooked$winner == 'R', cooked$r.rel, cooked$b.rel) 
cooked$pc <- ifelse(cooked$winner == 'R', cooked$pc.R, cooked$pc.B)
cooked$pc2 <- ifelse(cooked$winner == 'R', cooked$pc2.R, cooked$pc2.B)
cooked$pc3<-ifelse(cooked$winner == 'R', cooked$pc3.R, cooked$pc3.B)
cooked<-select_at(cooked,vars(pair,trt,temp,clip,winner,max,pc,pc2,pc3))
cooked<-na.omit(cooked)
#genets in clipping only
clipped<- dat20[dat20$trt=='C',] 
clipped$winner <- ifelse(clipped$b.rel > clipped$r.rel, 'B', 'R') 
clipped$max <- ifelse(clipped$winner == 'R', clipped$r.rel, clipped$b.rel) 
clipped$pc <- ifelse(clipped$winner == 'R', clipped$pc.R, clipped$pc.B) 
clipped$pc2 <- ifelse(clipped$winner == 'R', clipped$pc2.R, clipped$pc2.B) 
clipped$pc3 <- ifelse(clipped$winner == 'R', clipped$pc3.R, clipped$pc3.B) 
clipped<-select_at(clipped,vars(pair,trt,temp,clip,winner,max,pc,pc2,pc3))
clipped<-na.omit(clipped)
#genets in warming + clipping
cookgoose<- dat20[dat20$trt=='D',] 
cookgoose$winner <- ifelse(cookgoose$b.rel > cookgoose$r.rel, 'B', 'R')
cookgoose$max <- ifelse(cookgoose$winner == 'R', cookgoose$r.rel, cookgoose$b.rel) 
cookgoose$pc<-ifelse(cookgoose$winner == 'R', cookgoose$pc.R, cookgoose$pc.B) 
cookgoose$pc2<-ifelse(cookgoose$winner == 'R', cookgoose$pc2.R, cookgoose$pc2.B) 
cookgoose$pc3<-ifelse(cookgoose$winner == 'R', cookgoose$pc3.R, cookgoose$pc3.B) 
cookgoose<-select_at(cookgoose,vars(pair,trt,temp,clip,winner,max,pc,pc2,pc3))
cookgoose<-na.omit(cookgoose)
#combine all together
dominants<-rbind(control,cooked,clipped,cookgoose)

#looking at trait value of genotype with relative abundance < 0.5

#genets in control
controlL<- dat20[dat20$trt=='A',] 
controlL$loosers<- ifelse(controlL$b.rel < controlL$r.rel, 'B', 'R')
controlL$min <- ifelse(controlL$loosers == 'B', controlL$b.rel, controlL$r.rel) 
controlL$pc <- ifelse(controlL$loosers == 'B', controlL$pc.B, controlL$pc.R)
controlL$pc2<- ifelse(controlL$loosers == 'B', controlL$pc2.B, controlL$pc2.R)
controlL$pc3<- ifelse(controlL$loosers == 'B', controlL$pc3.B, controlL$pc3.R)
controlL<-select_at(controlL,vars(pair,trt,temp,clip,loosers,min,pc,pc2,pc3))
controlL<-na.omit(controlL)
#genets in warming only
cookedL<- dat20[dat20$trt=='B',] 
cookedL$loosers <- ifelse(cookedL$b.rel < cookedL$r.rel, 'B', 'R') 
cookedL$min <- ifelse(cookedL$loosers == 'B', cookedL$b.rel, cookedL$r.rel) 
cookedL$pc <- ifelse(cookedL$loosers == 'B', cookedL$pc.B, cookedL$pc.R)
cookedL$pc2 <- ifelse(cookedL$loosers == 'B', cookedL$pc2.B, cookedL$pc2.R)
cookedL$pc3<-ifelse(cookedL$loosers == 'B', cookedL$pc3.B, cookedL$pc3.R)
cookedL<-select_at(cookedL,vars(pair,trt,temp,clip,loosers,min,pc,pc2,pc3))
cookedL<-na.omit(cookedL)
#genets in clipping only
clippedL<- dat20[dat20$trt=='C',] 
clippedL$loosers <- ifelse(clippedL$b.rel < clippedL$r.rel, 'B', 'R') 
clippedL$min <- ifelse(clippedL$loosers == 'B', clippedL$b.rel, clippedL$r.rel) 
clippedL$pc <- ifelse(clippedL$loosers == 'B', clippedL$pc.B, clippedL$pc.R) 
clippedL$pc2 <- ifelse(clippedL$loosers == 'B', clippedL$pc2.B, clippedL$pc2.R) 
clippedL$pc3 <- ifelse(clippedL$loosers == 'B', clippedL$pc3.B, clippedL$pc3.R) 
clippedL<-select_at(clippedL,vars(pair,trt,temp,clip,loosers,min,pc,pc2,pc3))
clippedL<-na.omit(clippedL)
#genets in warming + clipping
cookgooseL<- dat20[dat20$trt=='D',] 
cookgooseL$loosers <- ifelse(cookgooseL$b.rel < cookgooseL$r.rel, 'B', 'R') 
cookgooseL$min <- ifelse(cookgooseL$loosers == 'B', cookgooseL$b.rel, cookgooseL$r.rel) 
cookgooseL$pc<-ifelse(cookgooseL$loosers == 'B', cookgooseL$pc.B, cookgooseL$pc.R) 
cookgooseL$pc2<-ifelse(cookgooseL$loosers == 'B', cookgooseL$pc2.B, cookgooseL$pc2.R) 
cookgooseL$pc3<-ifelse(cookgooseL$loosers == 'B', cookgooseL$pc3.B, cookgooseL$pc3.R) 
cookgooseL<-select_at(cookgooseL,vars(pair,trt,temp,clip,loosers,min,pc,pc2,pc3))
cookgooseL<-na.omit(cookgooseL)
#combine all together
loosers<-rbind(controlL,cookedL,clippedL,cookgooseL)

#winning genets have relative abundance of biomass > 0.7

dominants$trt<-factor(dominants$trt,labels=c("Control","Warmed","Clipped","Warmed + Clipped"))
bigwinners<-dominants[dominants$max>0.7,]
sumdat=ddply(bigwinners,c("trt"), summarize,
             Npc=length(pc),
             meanpc=mean(pc),
             sdpc=sd(pc),
             sepc=sdpc/sqrt(Npc),
             Npc2=length(pc2),
             meanpc2=mean(pc2),
             sdpc2=sd(pc2),
             sepc2=sdpc2/sqrt(Npc2))
sumdat$state<-c("Winner > 0.7","Winner > 0.7","Winner > 0.7","Winner > 0.7")
sumdat$trt <- factor(sumdat$trt, labels = c("Control (19)", "Warmed (22)","Clipped (23)","Warmed +\nClipped (15)"))
custom_breaks.x<-seq(-1,1,0.1)
custom_breaks.y<-round(seq(-0.8,1,0.1),2)
pc.plot<-ggplot(sumdat,aes(x=meanpc,y=meanpc2,shape=trt,fill=trt))+
  geom_point(cex=3)+
  theme_classic()+
  ylab("PC1: Photosynthetic efficiency (alpha)")+
  xlab(NULL)+
  scale_shape_manual(NULL,values=c(21,21,24,24))+
  scale_fill_manual(NULL,values=c("white","black","white","black"))+
  scale_color_manual(NULL,values=c("black","gray","white"))+
  labs(colour="Long title shortened\nwith wrapping") +
  geom_errorbar(aes(ymin=meanpc2-sepc2,ymax=meanpc2+sepc2))+
  geom_errorbarh(aes(xmax = meanpc + sepc, xmin = meanpc - sepc))+
  scale_x_continuous(limits=c(-0.8,0.6),breaks = custom_breaks.x,
                     labels = every_nth(custom_breaks.x,2, inverse = TRUE))+
  scale_y_continuous(limits=c(-0.8,0.8),breaks = custom_breaks.y,
                     labels = every_nth(custom_breaks.y,2, inverse = TRUE))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.title.x = element_text(margin = margin(t = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position=c(0.8,0.23))

sumdat$trt <- factor(sumdat$trt, labels = c("Control", "Warmed","Clipped","Warmed + Clipped"))


##big losers have a relative abundance of biomass < 0.3 
loosers$trt<-factor(loosers$trt,labels=c("Control","Warmed","Clipped","Warmed + Clipped"))
bigloosers<-loosers[loosers$min<0.3,]
sumlos=ddply(bigloosers,c("trt"), summarize,
             Npc=length(pc),
             meanpc=mean(pc),
             sdpc=sd(pc),
             sepc=sdpc/sqrt(Npc),
             Npc2=length(pc2),
             meanpc2=mean(pc2),
             sdpc2=sd(pc2),
             sepc2=sdpc2/sqrt(Npc2))
sumlos$state<-c("Loser < 0.3","Loser < 0.3","Loser < 0.3","Loser < 0.3")

##coexists 
medwinners<-dominants[dominants$max<0.7,]
medwinners<-medwinners %>% 
  rename(
    abund = max,
    even = winner)
medloosers<-loosers[loosers$min>0.3,]
medloosers<-medloosers %>% 
  rename(
    abund = min,
    even = loosers)
medall<-rbind(medwinners,medloosers)
medall$trt<-factor(medall$trt,labels=c("Control","Warmed","Clipped","Warmed + Clipped"))
sumeven=ddply(medall,c("trt"), summarize,
              Npc=length(pc),
              meanpc=mean(pc),
              sdpc=sd(pc),
              sepc=sdpc/sqrt(Npc),
              Npc2=length(pc2),
              meanpc2=mean(pc2),
              sdpc2=sd(pc2),
              sepc2=sdpc2/sqrt(Npc2))
sumeven$state<-c("Coexist = 0.3-0.7","Coexist = 0.3-0.7","Coexist = 0.3-0.7","Coexist = 0.3-0.7")

#combining them all together
all<-rbind(sumdat,sumlos,sumeven)
all$state <- ordered(all$state, 
  levels=c("Winner > 0.7","Coexist = 0.3-0.7", "Loser < 0.3"))
custom_breaks<-round(seq(-0.8,0.8,0.1),2)
all.plot<-ggplot(all,aes(x=meanpc,y=meanpc2,label=trt,fill=state))+
  ylab(NULL)+
  xlab(NULL)+
  scale_fill_manual(NULL,values=c("black","gray","white"))+
  scale_color_manual(NULL,values=c("black","gray","white"))+
  scale_x_continuous(limits=c(-0.81,0.81),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,4, inverse = TRUE))+
  scale_y_continuous(limits=c(-0.81,0.81),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,4, inverse = TRUE))+
  geom_errorbar(aes(ymin=meanpc2-sepc2,ymax=meanpc2+sepc2))+
  geom_errorbarh(aes(xmax = meanpc + sepc, xmin = meanpc - sepc))+
  geom_point(cex=3,shape=22)+
  facet_wrap(~trt)+
  geom_text(aes(label=Npc),hjust=2.6,vjust=2.57,size=1.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.length=unit(-0.1, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.title.x = element_text(margin = margin(t = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position="top",
        legend.spacing.x = unit(0.1, 'cm'))

#Figure 3
ggdraw() +
  draw_plot(pc.plot, x = 0.01, y = 0.1, width = .47, height = 0.9) +
  draw_plot(all.plot, x = .47, y = 0.1, width = .47, height = 0.9) +
  draw_text("PC2: Inverse (ETRmax + above:below)",size=10,x=0.5,y=0.07)+
  draw_plot_label(label = c("A)", "B)"), size = 10,
                  x = c(0, 0.47), y = c(1, 1))  
ggsave("FE-2020-00384_Fig3.png",height=4,width=7,dpi=300)
ggsave("FE-2020-00384_Fig3.pdf",height=4,width=7,dpi=300)

##########
#evenness
##########

#make a column "min" that has the relative abudnance for non-dominant genotype
for (i in 1:nrow(dat10)) {
  dat10$min[i] <- min(dat10[i,c("b.rel","r.rel")])
}
dat11<-dat10[abs(dat10$min)=="0",]
#stats (permutation analysis)
numperm=1000
print(anova(lm(min ~ temp*clip, data = dat10)))
retval <- matrix(0, ncol = 3, nrow = numperm + 1)
retval[1, ] <- anova(lm(min ~ temp*clip, data = dat10))[1:3, 4]
for (i in 1:numperm) {
  dat10$permtemp <- sample(dat10$temp, length(dat10$temp), replace = F)
  dat10$permclip <- sample(dat10$clip, length(dat10$clip), replace = F)
  retval[i + 1, ] <- anova(lm(min ~ permtemp*permclip, data = dat10))[1, 4]
}
par(mfrow = c(1, 3))
probALL = NULL
for (i in 1:3) {
  prob <- sum(retval[-1, i] >= retval[1, i])/(dim(retval)[1] - 1)
  probALL <- rbind(probALL, prob)
  hist(retval[-1, i],main = paste("Observed F-statistic =", round(retval[1, i],2),";", "p =",prob),
       xlab = "Simulated F-statistic")
  points(c(retval[1, i], retval[1, i]), c(0, numperm), type = "l", lwd = 2, 
         col = "gray",lty=2)
}
#figure 4 
custom_breaks<-seq(0,0.5,0.1)
ggplot(dat10,aes(x=clip,y=min,fill=as.factor(temp)))+
  geom_boxplot(outlier.colour=NA)+
  geom_dotplot(binaxis="y",binwidth=0.01,stackdir="center",position=position_dodge(0.8))+
  scale_fill_manual(NULL,labels = c("Ambient", "Warmed"),values = c("gray","gray40")) +
  theme_classic()+
  ylab("Evenness")+ 
  xlab(NULL)+
  scale_x_discrete(labels=c("Not Clipped","Clipped"))+
  scale_y_continuous(limits=c(0,0.5),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,1, inverse = TRUE))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.title.x = element_text(margin = margin(r = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position="top")
      
ggsave("FE-2020-00384_Fig4.png", dpi=300, units="in", width=3.5, height=3.5)
ggsave("FE-2020-00384_Fig4.pdf", dpi=300, units="in", width=3.5, height=3.5)

#########
#Conceptual figure 5
#########
con<-read.csv("cookedgoose_conceptual_dat.csv")
con<-melt(con,id.vars=c(1:2))

#both disturbances
con1<-con[con$plot=="1",]
pd <- position_dodge(0.1) 
p1<-ggplot(con1,aes(x=time,y=value,linetype=variable))+
  theme_classic()+
  geom_smooth(position=pd,lwd=1.2,color="black")+
  scale_linetype(labels=c("Genotype 1","Genotype 2"))+
  ylim(c(0,1))+
  ylab("Relative abundance")+
  xlab(NULL)+
  scale_x_continuous(limits=c(1.5,4.1),breaks=c(2,3), 
                     labels=c("Disturbance 1","Disturbance 2"))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position="none")

#disturbance 1 only
con2<-con[con$plot=="2",]
p2<-ggplot(con2,aes(x=time,y=value,linetype=variable))+
  theme_classic()+
  geom_smooth(position=pd,lwd=1.2,color="black")+
  scale_colour_manual(NULL,labels = c("Genotype 1", "Genotype 2"),values = c("gray","gray40")) +
  ylim(c(0,1))+
  ylab("Relative abundance")+
  xlab(NULL)+
  scale_x_continuous(limits=c(1.5,4.1),breaks=c(2), 
                     labels=c("Disturbance 1"))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position = "none")

#disturbance 2 only
con3<-con[con$plot=="3",]
p3<-ggplot(con3,aes(x=time,y=value,linetype=variable))+
  theme_classic()+
  geom_smooth(position=pd,lwd=1.2,color="black")+
  scale_colour_manual(NULL,labels = c("Genotype 1", "Genotype 2"),values = c("gray","gray40")) +
  ylim(c(0,1))+
  ylab("Relative abundance")+
  xlab(NULL)+
  scale_x_continuous(limits=c(1.5,4.1),breaks=c(3), 
                     labels=c("Disturbance 2"))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(r = 6),size=10),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.position = "none")

#combining all the panels into figure 5
plot_grid(p2,p3,p1,ncol=1, 
          align = "v", axis = "l",labels=c("A)","B)","C)"),label_size=10)
ggsave("FE-2020-00384_Fig5.png",dpi=300,units="in",width=3,height=6)
ggsave("FE-2020-00384_Fig5.pdf",dpi=300,units="in",width=3,height=6)

##########
#Temperature profile  
##########

#reading in the data
dat21<-read.csv("cookedgoose_temp_dat.csv")
dat21<-dat21[,(1:4)]
dat21$Date<-strptime(dat21$Date, format="%m/%d/%Y %H:%M")
dat21$Date<-as.POSIXct(dat21$Date)
dat21$Day<-format(dat21$Date, format="%m/%d/%Y")
dat21$Time<-format(dat21$Date,format="%H:%M")
dat21$Trt <- factor(dat21$Trt, levels = c("hot", "cold"))
avgtemp<-aggregate(dat21$Temp,
                   by=list(date=dat21$Date,trt=dat21$Trt),FUN=mean)
colnames(avgtemp)[colnames(avgtemp)=="x"] <- "mean"
avgtemp<-avgtemp[complete.cases(avgtemp),]
dat22<-avgtemp[avgtemp$date>="2017-10-03" & avgtemp$date<="2017-11-12",]
names(dat22)[3] <- "temp"
sumtemp<- ddply(dat22, c("trt"), summarise, 
                 N=length(temp),
                 mean=mean(temp),
                 min=min(temp),
                 max=max(temp))
custom_breaks<-seq(10,22,1)
ggplot(avgtemp,aes(x=date, y=mean,col=Trt))+
  geom_line(aes(col=trt))+
  scale_colour_manual(NULL,labels = c("Warmed", "Ambient"),values = c("black","gray70"))+
  theme_classic()+ 
  ylab("Mean temperature (degrees C)") +  xlab("Time")+
  scale_y_continuous(limits=c(10,22),breaks = custom_breaks,
                     labels = every_nth(custom_breaks,2, inverse = TRUE))+
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text=element_text(size=16),
        axis.title.y = element_text(margin = margin(r = 6),size=18),
        axis.title.x = element_text(margin = margin(r = 10),size=18),
        axis.text.y=element_text(margin=margin(r=10)),
        axis.text.x=element_text(margin=margin(t=10)),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.text=element_text(size=18))
ggsave("Fig_S1.png", dpi=300, units="in", width=8.5, height=11)

