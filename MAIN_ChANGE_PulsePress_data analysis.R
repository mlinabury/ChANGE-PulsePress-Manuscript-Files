####
####
####
####ChANGE Pulse-Press Analysis

library(plyr)
library(emmeans)
library(lmerTest)
library(lme4)
library(MuMIn)
library(tidyverse) #includes ggplot2, dplyr, reshape2
library(readxl)
library(codyn)
library(vegan)
library(nlme) #*n*lme to specify correlation structure
library(performance) #look at ANOVA assumptions
library(ggpubr)
library(MetBrewer)
library(multcomp)
library(car)
library(partR2)
library(ggeffects)
library(cowplot)



rm(list=ls())

#Set-up session:
options(contrasts=c('contr.sum','contr.poly'))

#Kim's home-built function for plotting, thank you Kim!
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  


setwd("C:/Users/mlinabu1/OneDrive - Colostate/Documents/Smith Lab/ChANGE/Data and Analysis/ChANGE pulse-press manuscript data and files")

#color palette: VanGogh1 from MetBrewer
  #0 -> 30g/m2: "gray60","#6f9954", "#89ab7c", "#87bcbd", "#969bc7", "#6b6ca3", "#434475", "#2c2d54"


####Most-abundant control species for background
######

#Identify the most abundant species during the study period with relative control abundance
sppcomp <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(site=="sgs",year>2020&year<2023,nitrogen==0,treatment=="nitrogen")%>%  #subset for 2020-2022 data
  dplyr::select(-treatment)

#Calculate bareground value so all total cover is above 100
#Sum total cover in plots
sppcomp_totalcover<-sppcomp %>% 
  dplyr::group_by(block, plot, site, year, nitrogen)%>%
  dplyr::summarise(totalcover=sum(cover))%>% 
  dplyr::as_tibble()

#Subtract totalcover value from 100
sppcomp_totalcover$bareground<-100-sppcomp_totalcover$totalcover

#Filter rows where bareground cover is more than or equal to 1
sppcomp_totalcover_pos<-filter(sppcomp_totalcover, bareground>=1)

#Gather the bareground col into a "species" column 
sppcomp_totalcover_pos<-gather(sppcomp_totalcover_pos,species,cover, bareground)

#delete total cover col for binding
sppcomp_totalcover_pos$totalcover<-NULL

#Bind the bareground dataset to the initial species composition dataset
sppcomp_bareground<-rbind(sppcomp,sppcomp_totalcover_pos)

#Check dimensions
dim(sppcomp);dim(sppcomp_totalcover_pos);dim(sppcomp_bareground) #Good, 14 bareground rows added

#Calculate totalcover again with new bareground variable, now there is no total cover value below 100
totalcover<-sppcomp_bareground %>% 
  dplyr::group_by(block, plot, site, year, nitrogen)%>%
  dplyr::summarise(plot_totalcover=sum(cover))%>% 
  dplyr::as_tibble()

#Merge focal species data to totalcover df
sppcomp_rel <- merge(sppcomp, totalcover, by=c("year", "site", "plot", "block", "nitrogen")) #merge documents based on block and plot

#Calculate relative cover by dividing cover by total cover (to get percent cover), then multiply by 100
sppcomp_rel$cover_rel<-(sppcomp_rel$cover/sppcomp_rel$plot_totalcover)*100

#Assess in table
abundance<-sppcomp_rel %>% 
  dplyr::group_by(site, species, nitrogen)%>%
  dplyr::summarise(cover_rel=sum(cover_rel)) %>%
  dplyr::as_tibble()

######

####Ammonium and Nitrate: analysis and visualization
######

nitrogen_x <- read.csv("MAIN_ChANGE_PulsePress_SGS_corrected ammonium nitrate_2021to2022.csv")

#Remove outlier
  #one nitrate value was off the charts without explanations: sgs, C23, 2021, 7.5, nitrogen plot
  #removed all 2021 C23 data, so calculations aren't effected (total N and percent change)
nitrogen<-nitrogen_x %>%
  filter(outlier!="x") %>%
  dplyr::select(-outlier)

#separate block_plot col into two
nitrogen<-separate(nitrogen, col=plot, into=c('block', 'plot'), sep=1)

#Rename treatments and channels for graphing ease
nitrogen$treatment[nitrogen$treatment == "n"]<-"Nitrogen"
nitrogen$treatment[nitrogen$treatment == "pp1"]<-"N+Deluge"
nitrogen$treatment[nitrogen$treatment == "pp2"]<-"N+Deluge"
nitrogen$channel[nitrogen$channel == "nitrate"]<-"Nitrate"
nitrogen$channel[nitrogen$channel == "ammonium"]<-"Ammonium"

nitrogen$treatment <- factor(nitrogen$treatment, levels=c("Nitrogen", "N+Deluge"))

#Spread data, so we can calculate total N
nitrogen_wide<-spread(nitrogen,channel,avgcalc_mg_l)
nitrogen_wide$"Total N"<-(nitrogen_wide$Ammonium+nitrogen_wide$Nitrate)

#Gather Nitrate and total N data back together into one col (also drop the ammonium col now)
nitrogen_plot<-nitrogen_wide%>%
  #dplyr::select(-Ammonium)%>%
  gather(Type, avgcalc_mg_l,c("Total N":"Nitrate"))

#fix structure
nitrogen_test<-nitrogen_plot %>% 
  dplyr::mutate_at(c("block","plot","site","nitrogen","sample_year","treatment","Type"),as.factor)
str(nitrogen_test) 

#Replicated mixed-effect model with a split-plot design and blocking, followed by type III ANOVA
  #Fixed effects: Nitrogen, Deluge
  #Random effects: Block, Plot
  #Note: We are only analyzing data within each year to account for the high MAP variability
  #Note: No unique "plotID" (i.e., block*plot) is needed since plots are already unique values (1-48), this unique term captures repeated measures

nitrogen_byyear <- lmer(avgcalc_mg_l~nitrogen*treatment #nitrogen and deluge treatment as fixed
                    + (1|block) + (1|plot), #block and plot as random, (1|plot) reflects the split-plot design
                    data=subset(nitrogen_test, sample_year==2022&Type=="Total N")) #subset CPER in a single year

check_model(nitrogen_byyear)
anova(nitrogen_byyear, ddf="Kenward-Roger") #ddf="Kenward-Roger" to calculate Type III test for mixed models

#Pairwise comparisons
nitrogen_by_treatment<-emmeans(nitrogen_byyear,pairwise~nitrogen|treatment)
cld(nitrogen_by_treatment$emmeans)

#Create a data frame for the cld
cld_n_trt_2021<-nitrogen_plot %>%
  filter(sample_year==2021)%>%
  ungroup() %>% 
  dplyr::select(sample_year,nitrogen,treatment,Type) %>%
  unique()%>%
  arrange(sample_year,treatment,Type,nitrogen)

cld_n_trt_2022<-nitrogen_plot %>%
  filter(sample_year==2022)%>%
  ungroup() %>% 
  dplyr::select(sample_year,nitrogen,treatment,Type) %>%
  unique()%>%
  arrange(sample_year,treatment,Type,nitrogen)


#manually add the cld letters to the data frame 
cld_n_trt_2021$cld_n_trt_2021<-c("a","a","a","a","a","a","a","a",   "a","ab","ab","ab","ab","ab","ab","b", #2021 N: Nitrate, Total N
                                   "a","ab","abc","bcd","abc","cd","d","d",   "a","ab","ab","bc","ab","cd","cd","d") #2021 PP: Nitrate, Total N

cld_n_trt_2022$cld_n_trt_2022<-c("a","a","a","a","a","a","a","a",   "a","a","a","a","ab","ab","b","c", #2022 N: Nitrate, Total N
                                   "a","a","ab","bcd","abc","cd","d","cd",   "a","a","a","abc","ab","bcd","cd","d") #2022 PP: Nitrate, Total N

#Separate data into years then merge cld data
nitrogen_plot_2021<-filter(nitrogen_plot, sample_year=="2021")
nitrogen_plot_2021_cld <- merge(nitrogen_plot_2021, cld_n_trt_2021, by=c("sample_year","treatment","nitrogen","Type"))

nitrogen_plot_2022<-filter(nitrogen_plot, sample_year=="2022")
nitrogen_plot_2022_cld <- merge(nitrogen_plot_2022, cld_n_trt_2022, by=c("sample_year","treatment","nitrogen","Type"))

#Plot absolute values
nplot_a_data<-barGraphStats(data=subset(nitrogen_plot_2021_cld),
                            variable="avgcalc_mg_l",
                            byFactorNames=c("nitrogen","sample_year","Type","treatment","cld_n_trt_2021")) %>%
  mutate(label_y = ifelse(Type == "Total N", mean+se+40, mean-se-40))  # adjust offset as needed

nplot_a<-ggplot(data=nplot_a_data, aes(x=as.factor(nitrogen), y=mean, color=Type)) +   geom_hline(yintercept = 0, size = 1, linetype="dashed") +
  geom_point(size=5) +
  facet_grid(sample_year~treatment) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5) +
  geom_text(aes(label = cld_n_trt_2021, y = label_y, color=Type), size=4) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        axis.line.x = element_line(linewidth = 1.5, colour="black"), 
        axis.line.y = element_line(linewidth = 1.5, colour="black"),
        axis.ticks=element_line(color="black",linewidth=1.5), 
        axis.text.x=element_text(color="black", size=12),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(color = "black", size = 25, hjust = 0),
        axis.title = element_blank(), 
        legend.position = c(.2,.7),
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1),
        strip.text = element_text(color="black",size = 20))+
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray50","darkorange3"))+ #alternate palette
  geom_vline(aes(xintercept= -Inf), size=2) + 
  geom_hline(yintercept = -Inf, size = 2) +
  labs(y=expression(Availibile~Nitrogen~(mg~l))) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ylim(-40,750) +
  ggtitle("a.")


nplot_b_data<-barGraphStats(data=subset(nitrogen_plot_2022_cld),
                            variable="avgcalc_mg_l",
                            byFactorNames=c("nitrogen","sample_year","Type","treatment","cld_n_trt_2022")) %>%
  mutate(label_y = ifelse(Type == "Total N", mean+se+50, mean-se-50))  # adjust offset as needed

nplot_b<-ggplot(data=nplot_b_data, aes(x=as.factor(nitrogen), y=mean, color=Type)) + 
  geom_hline(yintercept = 0, size = 1, linetype="dashed") +
  geom_point(size=5) +
  facet_grid(sample_year~treatment) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5, ) +
  geom_text(aes(label = cld_n_trt_2022, y = label_y, color=Type), size=4) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        axis.line.x = element_line(linewidth = 1.5, colour="black"), 
        axis.line.y = element_line(linewidth = 1.5, colour="black"),
        axis.ticks=element_line(color="black",linewidth=1.5), 
        axis.text.x=element_text(color="black", size=12),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(color = "black", size = 25, hjust = 0),
        axis.title.x = element_text(color="black",size = 20), 
        axis.title.y = element_blank(), 
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1),
        strip.text = element_text(color="black",size = 20))+
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray50","darkorange3"))+ #alternate palette
  geom_vline(aes(xintercept= -Inf), size=2) + 
  geom_hline(yintercept = -Inf, size = 2) +
  labs(y=expression(Available~Nitrogen~(mg~l))) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ylim(-40,865) +
  ggtitle("b.")



nplot<-ggarrange(nplot_a,nplot_b, ncol=1,nrow=2, heights = c(1,1.15))

annotate_figure(nplot,
                left = text_grob(expression(Available~Nitrogen~(mg~l^-1)),color="black",size=20, vjust = .4,rot = 90))

#Change magnitude
ddply(subset(nitrogen_plot, sample_year=="2022"), c("Type","treatment"), summarize, 
      n= length (avgcalc_mg_l),
      mean = mean(avgcalc_mg_l),
      sd = sd(avgcalc_mg_l),
      se = sd/(sqrt(n)))

######

####Green up: analysis and visualization
######

#Load in data
green <- read.csv("MAIN_ChANGE_PulsePress_SGS_calculated green up values_2021to2022.csv")%>%
  dplyr::select(-X,-format)

#Prep for plotting
green$treatment[green$treatment == "n"]<-"Nitrogen"
green$treatment[green$treatment == "pp1"]<-"N+Deluge"
green$treatment[green$treatment == "pp2"]<-"N+Deluge"

#Add first and last deluge application dates
green_covar<-green %>%
  filter(year=="2021"|year=="2022")%>%
  ungroup() %>% 
  dplyr::select(year) %>%
  unique()

green_covar$first_deluge_doy<-c(208,174) 
green_covar$last_deluge_doy<-c(217,179)

green<-merge(green, green_covar, by=c("year"))

#Define dates
green$date<-ymd(green$date)

green$deluge.date.1<-mdy(green$deluge.date.1)
green$deluge.date.1.doy <- strftime(green$deluge.date.1, format = "%j")
green$deluge.date.1.doy<-as.numeric(green$deluge.date.1.doy)

green$deluge.date.2<-mdy(green$deluge.date.2)
green$deluge.date.2.doy <- strftime(green$deluge.date.2, format = "%j")
green$deluge.date.2.doy<-as.numeric(green$deluge.date.2.doy)

#Fix structure
green<-green %>% 
  dplyr::mutate_at(c("site","block","plot","treatment","nitrogen","year"),as.factor)
str(green)

#Order factors
green$treatment <- factor(green$treatment, levels=c("Nitrogen", "N+Deluge"))

#Replicated mixed-effect model with a split-plot design and blocking, followed by type III ANOVA
  #Fixed effects: Nitrogen, Deluge
  #Random effects: Block, Plot
  #Note: We are only analyzing data within each year to account for the high MAP variability
  #Note: No unique "plotID" (i.e., block*plot) is needed since plots are already unique values (1-48), this unique term captures repeated measures
  #lme() is needed here to account for repeated measures and a correlation structure

func_anppSGS <-lme(gcc_mean~nitrogen*treatment*as.factor(doy), #nitrogen and year as fixed
                   data=subset(green, year==2022), #subset for year of interest
                   random=~1|block/plot, #block as random effect (with plot in block)
                   correlation = corCompSymm(form = ~doy | block/plot), #when assessing over time "year" could be here, but doy is the within year equivelent
                   control=lmeControl(returnObject=T)) #not doing anything probably

check_model(func_anppSGS)
summary(func_anppSGS)
anova.lme(func_anppSGS, type = "marginal")

#Pairwise comparisons
em_shortgrass_nitrogen_by_treatment<-emmeans(func_anppSGS,pairwise~treatment|doy)
pairs(em_shortgrass_nitrogen_by_treatment)

#Plotting
#Nitrogen is not significant in either model--average across nitrogen to block level to visualize deluge treatment effect
green_avg<-aggregate(green$gcc_mean,
                     by=list(site=green$site,
                             block=green$block,
                             treatment=green$treatment,
                             first_deluge_doy=green$first_deluge_doy,
                             last_deluge_doy=green$last_deluge_doy,
                             doy=green$doy,
                             year=green$year), FUN=mean)
names(green_avg)[names(green_avg)=="x"] <- "gcc_avg"

green_plot<-ggplot(data=barGraphStats(data=subset(green_avg), variable="gcc_avg", byFactorNames=c("treatment","year","doy","first_deluge_doy","last_deluge_doy")), aes(x=doy, y=mean, color=treatment)) + 
  geom_rect(aes(xmin=first_deluge_doy, xmax=last_deluge_doy, ymin=-Inf, ymax=Inf), color = NA, fill = "#87bcbd", alpha=.05) +
  geom_point(size=3.5) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=2) +
  #geom_smooth(method = "loess") + 
  facet_grid(year~.) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size=1.5, colour="black"), 
        axis.line.y = element_line(size=1.5, colour="black"),
        axis.ticks = element_line(color="black",size=1.5), 
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black",size=15),
        axis.title = element_text(color="black",size=20),
        plot.title = element_text(color="black", size=22, hjust=0),
        axis.title.y = element_text(vjust=1.2),
        strip.text.x = element_text(color="black",size=20),
        strip.text.y = element_text(color="black",size=20),
        legend.position = c(.25,.8),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1)) +
  scale_color_manual(values=c("gray60","#6f9954")) + #alternate palette
  #scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray60","#6f9954", "#89ab7c", "#87bcbd", "#969bc7", "#6b6ca3", "#434475", "#2c2d54"))+ #alternate palette
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  labs(y=expression(Greenness~(GCC))) + #y label
  labs(x=expression(Day~of~Year)) + #x label
  ggtitle("a.")

#plot average gcc by nitrogen treatment between deluge treatments to indicate interaction
green$treatment<-as.character(green$treatment)
green$treatment[green$treatment == "Nitrogen"]<-"N"
green$treatment[green$treatment == "N+Deluge"]<-"N+D"

green_in<-ggplot(data=barGraphStats(data=subset(green,year=="2022"), variable="gcc_mean", byFactorNames=c("nitrogen","treatment","year")), aes(x=nitrogen, y=mean, color=nitrogen)) +
  geom_point(size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5) +
  facet_grid(.~treatment,scales = "free") +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size=1.5, colour="black"), 
        axis.line.y = element_line(size=1.5, colour="black"),
        axis.ticks = element_line(color="black",size=1.5), 
        axis.text.x = element_text(color="black", size=8, angle = 90, vjust=.5),
        axis.text.y = element_text(color="black",size=8),
        plot.title = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(color="black",size=10),
        legend.position = "none", #c(.1,.78)
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1)) +
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray60","#6f9954", "#89ab7c", "#87bcbd", "#969bc7", "#6b6ca3", "#434475", "#2c2d54"))+ #alternate palette
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  labs(y=expression(GCC)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ggtitle("B.") +
  scale_y_continuous(breaks = seq(0.334,0.34, by = 0.0015), labels = scales::label_number(accuracy = 0.001))  # Adjust '0.2' to set the spacing of tick marks

green_final <-
  ggdraw() +
  draw_plot(green_plot) +
  draw_plot(green_in, x = 0.645, y = .2555, width = .3, height = .27)

green_final

######

####Soil Moisture: Analysis of nitrogen and N+Deluge plot differences
######

#Determine persistence of deluge treatment effects through t-tests

#1. Load in and filter for site, block, and date
#2. Test assumption: normally distributed data
  #If normally distributed, continue to step 3
  #If not normally distributed, do a non-parametric paired two-sample Wilcoxon test
#3. Test assumption: equal/unequal varience
  #If equal varience, do a pooled two-sample ttest
  #If unequal varience, do a welch-satterthwaithe ttest


#1. Load in and filter for site, year, and date (and block if interested)

vwc <- read.csv("MAIN_ChANGE_PulsePress_SGS_soil moisture_2021to2022.csv")%>%
  dplyr::select(-X)%>%
  filter(site=="sgs", year=="2022", days_since_deluge==1)

boxplot(vwc~treatment, data= vwc, ylab = "vwc")

#2. Test assumption: normally distributed data
shapiro.test(vwc$vwc) #Shapiro-Wilk normality test, null hypothesis: data are normally distributed.

#If normally distributed, continue to step 3
#If not normally distributed, do a non-parametric paired two-sample Wilcoxon test (with no assump of normality)
  #Null: no difference
wilcox.test(vwc~treatment, paired = TRUE, data = vwc)

#3. Test assumption: equal/unequal varience
var.test(vwc~treatment, data= vwc) #null: equal varience
leveneTest(vwc~treatment, data= vwc)  #null: equal varience
#If equal varience, do a pooled two-sample ttest
t.test(vwc~treatment, var.equal=TRUE, data= vwc)

#If unequal varience, do a welch-satterthwaithe ttestt
t.test(vwc~treatment, var.equal=FALSE, data= vwc)

ddply(subset(vwc), c("treatment"), summarize, 
      n= length (vwc),
      mean = mean(vwc),
      sd = sd(vwc),
      se = sd/(sqrt(n)))

######

####Soil Moisture: Visualization of water treatments
######

vwc <- read.csv("MAIN_ChANGE_PulsePress_SGS_soil moisture_2021to2022.csv")%>%
  dplyr::select(-X)%>%
  filter(site=="sgs")

#Simplify deluge treatment name--will separate PP1 and PP2 treatments later by year
vwc$treatment[vwc$treatment == "nitrogen"]<-"Nitrogen"
vwc$treatment[vwc$treatment == "pulse-press1"]<-"N + Deluge"
vwc$treatment[vwc$treatment == "pulse-press2"]<-"N + Deluge"

#create a column for duration of significant treatment differences
vwc$duration_of_diff<-ifelse(vwc$year==2021,15,25)

#fix structure
vwc<-vwc %>% 
  dplyr::mutate_at(c("site","block","plot","treatment","year"),as.factor)
str(vwc)

#Order factored treatments
vwc$treatment <- factor(vwc$treatment, levels=c("Nitrogen", "N + Deluge"))

#average vwc down to a single measure by day, by block
vwc_byday<-vwc %>% 
  dplyr::group_by(site,treatment,days_since_deluge,year,duration_of_diff)%>% #remove "Remove" column by leaving it out of the function
  dplyr::summarise(across(vwc,~mean(.x,na.rm=TRUE))) %>% 
  dplyr::as_tibble()

fig2b<-ggplot(data=subset(vwc_byday),  aes(x=days_since_deluge, y=vwc, color=treatment)) +
  geom_rect(aes(xmin=0, xmax=duration_of_diff, ymin=-Inf, ymax=Inf), color = NA, fill = "#FFCC99", alpha=.02) +
  geom_point(size=1.5) +
  geom_smooth(method = "gam", size=1.5) +
  facet_grid(year~.) +
  ylim(0,16) +
  xlim(-3,30) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=12, angle = 90),
        axis.text.y=element_text(color="black",size=15),
        axis.title=element_text(color="black",size=20),
        plot.title = element_text(color = "black", size = 22, hjust = 0),
        panel.grid.minor = element_line(colour = NA), 
        axis.title.y = element_text(vjust=1.2), 
        legend.position=c(.77,.88),
        legend.title = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1),
        strip.text.y = element_text(size = 20),
        legend.text = element_text(size=15)) +
  geom_vline(aes(xintercept=-Inf), size=3) + 
  geom_vline(aes(xintercept=0), size=2, linetype = "dashed") + 
  geom_hline(yintercept = -Inf, size = 2) +
  geom_hline(aes(yintercept=0), linewidth=1, linetype="dashed") +
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray50","#87bcbd","#87bcbd")) + #alternate palette
  labs(y="VWC (%)")+ #y label
  labs(x=expression(Days~since~deluge)) + #x label
  ggtitle("b.")

######

####Growing Season Precipitation: Visualization with deluge period
######

ppt <- read.csv("MAIN_ChANGE_PulsePress_SGS_precip_Jan1toAug31_2021to2022.csv")%>%
  dplyr::select(-X)

#Calculate seasonal totals 
ddply(subset(ppt), c("year"), summarize, 
      n= length (precipitation_mm),
      growing_season_total = sum(precipitation_mm),
      sd = sd(precipitation_mm),
      se = sd/(sqrt(n)))

#Create a dataset of the deluge periods then merge back to the main dataset
  #SGS- 2021: 7/27 to 8/7
  #SGS- 2022: 6/23 to 6/28
ppt_covar<-ppt %>%
  ungroup() %>% 
  dplyr::select(year) %>%
  unique()

ppt_covar$first_deluge_doy<-c(208,174) 
ppt_covar$last_deluge_doy<-c(219,179)

ppt<-merge(ppt, ppt_covar, by=c("year"))

#Create a second dataset for Deluge precipitation then add the 42 mm at point of deluge
ppt_trt<-ppt
ppt$treatment<-"Nitrogen"
ppt_trt$treatment<-"N + Deluge"

#Use an ifelse statement to add 42mm at the midpoint of deluge
  #2021 midpoint: 213.5 -> 214
  #2022 midpoint: 176.5 ->177
ppt_trt$precipitation_summed_mm<-ifelse(ppt_trt$year==2021&ppt_trt$doy>213,(ppt_trt$precipitation_summed_mm+42),(ppt_trt$precipitation_summed_mm=ppt_trt$precipitation_summed_mm))
ppt_trt$precipitation_summed_mm<-ifelse(ppt_trt$year==2022&ppt_trt$doy>176,(ppt_trt$precipitation_summed_mm+42),(ppt_trt$precipitation_summed_mm=ppt_trt$precipitation_summed_mm))

ppt_full<-full_join(ppt,ppt_trt)

#Use an ifelse statement to add a maximum precip value to the last day of the year
ppt_full$label<-ifelse(ppt_full$doy==243,(ppt_full$label=ppt_full$precipitation_summed_mm),(ppt_full$label=""))
ppt_full$label<-as.numeric(ppt_full$label)

#Plot
fig2a<-ggplot(data=subset(ppt_full,year!="2020"&year!="2023"&doy>122), aes(x=doy, y=precipitation_summed_mm, color=treatment)) + 
  geom_rect(aes(xmin=first_deluge_doy, xmax=last_deluge_doy, ymin=-Inf, ymax=Inf), color = NA, fill = "#87bcbd", alpha=.01) +
  geom_line(aes(y=precipitation_summed_mm),linewidth=1.6) +
  geom_text(aes(label = round(label, digits = 1)), hjust=0.7, vjust=-0.46, size=4, color="black") +
  guides(color = guide_legend(reverse=TRUE)) +
  facet_grid(year~.) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size=1.5, colour="black"), 
        axis.line.y = element_line(size=1.5, colour="black"),
        axis.ticks = element_line(color="black",size=1.5), 
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black",size=15),
        axis.title = element_text(color="black",size=20),
        plot.title = element_text(color="black", size=22, hjust=0),
        axis.title.y = element_text(vjust=1.2),
        strip.text.x = element_text(color="black",size=20),
        strip.text.y = element_text(color="black",size=20),
        legend.position = c(.5,.63),
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1)) +
        #scale_linetype_manual(values=c("dotted","solid"))+ #alternate palette
  scale_color_manual(values=c("#87bcbd","gray50")) + #alternate palette
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=0), linewidth=1, linetype="dashed") +
  labs(y=expression(Precipitation~(mm))) + #y label
  labs(x=expression(Day~of~Year)) + #x label
  ylim(0,310) +
  ggtitle("a.")

ggarrange(fig2a,fig2b,
                ncol = 2, nrow = 1,
                widths = c(1,1))

######

####Conceptual figure (base)
######

#Create a dataset for the conceptual figure

xlabel<-c("Nitrogen","Nitrogen","Deluge","Deluge","N+Deluge","N+Deluge","N+Deluge","N+Deluge")
type<-c("amb","n","amb","del","amb","n","del","syn")
value<-c(0,1,0,2,0,1,2,2)

concept <- data.frame(xlabel = xlabel, type = type, value = value)

#Fix structure to order factors
concept<-concept %>% 
  dplyr::mutate_at(c("xlabel","type"),as.factor)

concept$xlabel <- factor(concept$xlabel, levels=c("Ambient","Nitrogen","Deluge","N+Deluge"))
concept$type <- factor(concept$type, levels=c("amb","syn","del","n"))

str(concept) 

ggplot(data=subset(concept), aes(x=xlabel, y=value, fill=type)) + 
  geom_bar(stat="identity", position = "stack", width = 0.4) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        #panel.grid.major = element_line(colour = "gray", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks = element_line(color="black",size=1.5), 
        axis.text.x = element_text(color="black", size=20), # angle=40, hjust=1
        axis.text.y = element_text(color="black", size=15),
        plot.title = element_text(color="black",size=20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color="black",size=20),
        strip.text.x = element_text(size = 18),
        legend.position="none") +
  scale_fill_manual(values=c("black","gray50","#87bcbd","#434475"),name="Contribution by")+ #alternate palette
  geom_vline(aes(xintercept=-Inf), size=1.5) + 
  geom_hline(aes(yintercept=0), size=1, linetype="dashed") + 
  labs(y=expression(ANPP~Response)) +
  scale_y_continuous(breaks = seq(0,5,by=10), limits = c(0,7)) +
  ggtitle("a. Examples of interaction")

#Export plot to label in powerpoint

######




####Production: Quadratic Analysis and Visualization
######

#Load in dataset
anpp_all <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X)

#Correct facet wrap labels for consistency
anpp_all$treatment[anpp_all$treatment == "N"]<-"Nitrogen"
anpp_all$treatment[anpp_all$treatment == "PP1"]<-"N+Deluge"
anpp_all$treatment[anpp_all$treatment == "PP2"]<-"N+Deluge"

#fix structure
anpp_all<-anpp_all %>% 
  dplyr::mutate_at(c("block","plot","site","year","treatment"),as.factor)
str(anpp_all) 

#Replicated mixed-effect model with a split-plot design and blocking, followed by type III ANOVA
  #Fixed effects: Nitrogen, Deluge
  #Random effects: Block, Plot
  #Note: We are only analyzing data within each year to account for the high MAP variability
  #Note: No unique "plotID" (i.e., block*plot) is needed since plots are already unique values (1-48), this unique term captures repeated measures

mod_lin <- lmer(anpp~nitrogen*treatment #nitrogen and deluge treatment as fixed
                               + (1|block) + (1|plot), #block and plot as random, (1|plot) reflects the split-plot design
                               data=subset(anpp_all, year==2022), REML=FALSE) #subset CPER in a single year

#Center to reduce colinnearity between nitrogen and the quadratic nitrogen term
anpp_all$nitrogen_centered<-scale(as.numeric(as.character(anpp_all$nitrogen)), center=TRUE, scale=FALSE)

mod_quad <- lmer(anpp~(nitrogen_centered + I(nitrogen_centered^2))*treatment + #nitrogen and deluge treatment as fixed
                + (1|block) + (1|plot), #block and plot as random, (1|plot  ) reflects the split-plot design
                data=subset(anpp_all, year==2022), REML=TRUE) #subset CPER in a single year

#Quadratic model is a better fit, proceed with the quadratic model
anova(mod_lin, mod_quad)

#Analysis
check_model(mod_quad)
anova(mod_quad, ddf="Kenward-Roger") #ddf="Kenward-Roger" to calculate Type III test for mixed models

#Pairwise comparisons
emmeans(mod_quad,~treatment)
pairs(emmeans(mod_quad,~treatment))

#Calculate peaks in nitrogen
b1<-fixef(mod_quad)["nitrogen_centered"]
b2<-fixef(mod_quad)["I(nitrogen_centered^2)"]

anpp_2022<-subset(anpp_all, year==2022)
peak_centered<-(-b1/(2*b2))
peak_actual_N<-(peak_centered+mean(anpp_2022$nitrogen))
peak_actual_N

#Create a data frame for the main effects cld
cld_anpp_trt_years_main<-anpp_all %>%
  filter(site=="SGS")%>%
  ungroup() %>% 
  dplyr::select(site,treatment,year) %>%
  unique()%>%
  arrange(year,treatment)%>%
  filter(year==2021|treatment!="PP1")

cld_anpp_trt_years_main$cld_anpp_trt_years_main<-c("","", #2021 PP1 and N
                                         "b","a") #2022 PP2 and N

anpp_all <- merge(anpp_all, cld_anpp_trt_years_main, by=c("site","treatment","year"))

#Create a second dataframe for the main effects plot (prevents issues with factor naming)
anpp_all_main<-anpp_all
anpp_all_main$treatment<-as.character(anpp_all_main$treatment)
anpp_all_main$treatment[anpp_all_main$treatment == "Nitrogen"]<-"N"
anpp_all_main$treatment[anpp_all_main$treatment == "N+Deluge"]<-"N+D"

#Plotting
fig3_a<-ggplot(data=barGraphStats(data=subset(anpp_all_main, site=='SGS'), variable="anpp", byFactorNames=c("treatment","year","cld_anpp_trt_years_main")), aes(x=treatment, y=mean, fill=treatment)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  geom_text(aes(label = cld_anpp_trt_years_main, y = mean+se), vjust=-.5, size=4, color="black") +
  facet_grid(year~.,scales = "free") +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        axis.line.x = element_line(size=1.5, colour="black"), 
        axis.line.y = element_line(size=1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=20, vjust=0.5),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(size=25), #set as 25
        panel.grid.minor = element_line(colour=NA), 
        axis.title.y = element_text(color="black",size=20),
        axis.title.x = element_text(size=30), 
        legend.position="none",
        legend.text = element_text(size=15),
        legend.title = element_blank(),
        strip.text.x = element_text(color="black",size=20),
        strip.text.y = element_blank())+
  geom_vline(aes(xintercept=-Inf), size=1.5) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=0), linewidth=1, linetype="dashed") +
  ylim(0,185) +
  scale_fill_manual(values=c("gray60","#87bcbd"))+ #alternate palette
  labs(y=expression(ANPP~(g~m^-2)))+ #y label
  labs(x=" ") + #x label
  ggtitle("a.")

#Production across N treatments faceted by year and water treatment
anpp_all$treatment <- factor(anpp_all$treatment, levels=c("Nitrogen","N+Deluge"))

fig3_b<-ggplot(data=anpp_all, aes(x=nitrogen, y=anpp, fill=treatment, color=treatment)) +
  geom_smooth(method="lm", formula = y~x+I(x^2), se=TRUE) +
  facet_grid(year~.) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15, vjust = .4),
        axis.text.y=element_text(color="black",size=15),
        axis.title=element_text(color="black",size=20),
        plot.title = element_text(size=25),
        strip.text.x = element_text(color="black",size = 20),
        strip.text.y = element_text(color="black",size = 20),
        legend.position=c(.5,.62),
        legend.text = element_text(size=15),
        legend.title = element_blank()) +
  geom_vline(aes(xintercept=-Inf), size=1.5) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=0), linewidth=1, linetype="dashed") +
  #ylim(0,185) +
  scale_fill_manual(values=c("gray60","#87bcbd"))+ #alternate palette
  scale_color_manual(values=c("gray60","#87bcbd"))+ #alternate palette
  labs(y=expression(Predicted~ANPP~(g~m^-2)))+ #y label
  labs(x=expression(Nitrogen~(g~m^-2))) +  #x label
  ggtitle("b.")

#Arrange and annotate plots
anpp_main_fig<-ggarrange(fig3_a, fig3_b,
                ncol = 2, nrow = 1,
                widths = c(.5,1))
anpp_main_fig

#Calculate the magnitude increase of ANPP between control and deluge plots
anpp_all$fertilized<-ifelse(anpp_all$nitrogen=="0","control","nitrogen")
ddply(subset(anpp_all, year==2021), c("fertilized"), summarize, 
      n= length (anpp),
      mean = mean(anpp),
      sd = sd(anpp),
      se = sd/(sqrt(n)))

######

####Production: Analysis and visualization of absolute grass and forb production
######

#Load in dataset
anpp_all <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X)

#Correct labels for consistency
anpp_all$treatment[anpp_all$treatment == "N"]<-"N"
anpp_all$treatment[anpp_all$treatment == "PP1"]<-"N+D"
anpp_all$treatment[anpp_all$treatment == "PP2"]<-"N+D"

#fix structure
anpp_all<-anpp_all %>% 
  dplyr::mutate_at(c("block","plot","site","year","treatment"),as.factor)
str(anpp_all) 

#Replicated mixed-effect model with a split-plot design and blocking, followed by type III ANOVA
  #Fixed effects: Nitrogen, Deluge
  #Random effects: Block, Plot
  #Note: We are only analyzing data within each year to account for the high MAP variability
  #Note: No unique "plotID" (i.e., block*plot) is needed since plots are already unique values (1-48), this unique term captures repeated measures
  #ANPP is analyzed with a quadratic term, so grass and forb (the component parts) will be too

anpp_all$nitrogen_centered<-scale(as.numeric(as.character(anpp_all$nitrogen)), center=TRUE, scale=FALSE)

mod_quad <- lmer(grass~(nitrogen_centered + I(nitrogen_centered^2))*treatment + #nitrogen and deluge treatment as fixed
                   + (1|block) + (1|plot), #block and plot as random, (1|plot ) reflects the split-plot design
                 data=subset(anpp_all, year==2022), REML=TRUE) #subset CPER in a single year

#Analysis
check_model(mod_quad)
anova(mod_quad, ddf="Kenward-Roger") #ddf="Kenward-Roger" to calculate Type III test for mixed models

#Pairwise comparisons
emmeans(mod_quad,~treatment)
pairs(emmeans(mod_quad,~treatment))

#Reformat data for plot a
anpp_all <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X)

#Correct labels for consistency
anpp_all$treatment[anpp_all$treatment == "N"]<-"N"
anpp_all$treatment[anpp_all$treatment == "PP1"]<-"N+D"
anpp_all$treatment[anpp_all$treatment == "PP2"]<-"N+D"

#fix structure
anpp_all<-anpp_all %>% 
  dplyr::mutate_at(c("block","plot","site","year","treatment"),as.factor)

#Gather to long format
anpp_long<-anpp_all%>%
  gather(functional,anpp_func,"forb":"grass")

#Correct labels for consistency
anpp_long$functional[anpp_long$functional == "grass"]<-"Grass"
anpp_long$functional[anpp_long$functional == "forb"]<-"Forb"

#Need to fix error bars for stacked bar plots: grass is stacked on top, so we need to add the forb average to the grass
#The following code with generated by the help of the stat lab
data_func<-barGraphStats(data=subset(anpp_long, site=='SGS'), variable="anpp_func", byFactorNames=c("treatment","year","functional"))
data_func2 <- data_func %>%
  dplyr::select(-sd, -se) %>%
  tidyr::pivot_wider(names_from = functional, 
                     values_from = mean) %>%
  dplyr::mutate(Grass = Grass + Forb) %>%
  tidyr::pivot_longer(c(Forb, Grass), 
                      names_to = "functional", 
                      values_to = "mean_alt")

data_func_final <- dplyr::left_join(data_func, data_func2) %>%
  dplyr::mutate(functional = factor(functional, 
                                    levels = c("Grass", "Forb"))) %>%
  dplyr::mutate(treatment = factor(treatment, 
                                   levels = c("N", "N+D")))

#Only need the error bars from the one group (the grass since its stacked on top), so set the cld to nothing for the forb group
#in this type of plot, the errorbars of grass/forb mirror eachother, so we need to remove a set of errorbars
data_func_final$mean_alt<-ifelse(data_func_final$functional=="Forb",NA,data_func_final$mean_alt)
data_func_final$se<-ifelse(data_func_final$functional=="Forb",NA,data_func_final$se)

forbfig_a<-ggplot(data=subset(data_func_final), aes(x=treatment, y=mean, fill=functional)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean_alt-se, ymax=mean_alt+se), width=0.2) +
  facet_grid(year~"") +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        #panel.grid.major = element_line(colour = "gray"),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=20),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(size=25), #set as 25
        panel.grid.minor = element_line(colour = NA), 
        axis.title.y = element_text(color="black",size=20,vjust=1.2),
        axis.title.x = element_text(size=15,color="white"), 
        legend.position= c(.5,.4),
        legend.text = element_text(size = 18),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1),
        legend.title = element_blank(),
        strip.text.x = element_text(color="black",size = 20),
        strip.text.y = element_blank())+
  geom_vline(aes(xintercept=-Inf), size=1.5) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  ylim(0,185) +
  scale_fill_manual(values=c("#6f9954","#434475")) +
  labs(y=expression(ANPP~(g~m^-2)))+ #y label
  labs(x=expression(Nitrogen~(g~m^-2))) +  #x label
  ggtitle("a.")


#Plot across nitrogen levels
#Recreate dataset to prevent issues
anpp_all <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X)

#Correct facet wrap labels for consistency
anpp_all$treatment[anpp_all$treatment == "N"]<-"Nitrogen"
anpp_all$treatment[anpp_all$treatment == "PP1"]<-"N+Deluge"
anpp_all$treatment[anpp_all$treatment == "PP2"]<-"N+Deluge"

#fix structure
anpp_all<-anpp_all %>% 
  dplyr::mutate_at(c("block","plot","site","year","treatment"),as.factor)
str(anpp_all) 

#Plot fig b
anpp_plot <- anpp_all %>%
  pivot_longer(cols = c(grass, forb), names_to = "group", values_to = "anpp_func")

anpp_plot$treatment <- factor(anpp_plot$treatment, levels=c("Nitrogen", "N+Deluge"))

forbfig_b<-ggplot(data=subset(anpp_plot), aes(x=nitrogen, y=anpp_func, fill=group, color=group, group=group)) +
  geom_smooth(method="lm", formula = y~x+I(x^2), se=TRUE) +
  facet_grid(year~treatment,) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        #panel.grid.major = element_line(colour = "gray"),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=12, angle=90, vjust=.5),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(size=25), #set as 25
        panel.grid.minor = element_line(colour = NA), 
        axis.title = element_text(color="black",size=20,vjust=1.2),
        legend.position= c(.25,75),
        legend.text = element_text(size = 18),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1),
        legend.title = element_blank(),
        strip.text.x = element_text(color="black",size = 20),
        strip.text.y = element_text(color="black",size = 20))+
  geom_vline(aes(xintercept=-Inf), size=1.5) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  #ylim(0,185) +
  geom_hline(yintercept = 0, size = 1, linetype="dashed") +
  scale_fill_manual(values=c("#434475","#6f9954")) +
  scale_color_manual(values=c("#434475","#6f9954")) +
  labs(y=expression(Predicted~ANPP~(g~m^-2)))+ #y label
  labs(x=expression(Nitrogen~(g~m^-2))) +  #x label
  ggtitle("b.")

forbfig<-ggarrange(forbfig_a, forbfig_b,
                   ncol = 2, nrow = 1,
                   widths = c(.5,1))

forbfig

######

####Production: Visualization of observed v. expected 2022 response -- 2-panel plot
######

#The following is focusing on 2022 nitrogen/deluge effects to understand the nature of the interaction

#1. Calculate the nitrogen treatment effect
#must remove plot, or the spread function won't work (since plot is a unique value)
anpp_N <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X,-plot,-grass,-forb,-woody,-dead,-anpp_with_woody,-cactus)%>%
  filter(site=="SGS",treatment=="N",year==2022)

#reorganize data: spread the nitrogen column, then gather the 2.5-30 columns, creating a 0 g/m2 reference column
n_effect<-anpp_N%>%
  spread(nitrogen,anpp)%>%
  rename_at("0",~"control_anpp")%>%
  gather(nitrogen,anpp,"2.5":"30")
head(n_effect)

#create a treatment effect column of nitrogen-treated anpp - control anpp)
n_effect$n_effect_anpp<-(n_effect$anpp-n_effect$control_anpp)
head(n_effect)


#2. Calculate the deluge treatment effect
anpp_N_PP2 <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X,-plot,-grass,-forb,-woody,-dead,-anpp_with_woody,-cactus)%>%
  filter(site=="SGS",treatment!="PP1",year==2022)

#reorganize data: filter for only 0 nitrogen plots to isolate the deluge effect, spread the treatment columns
del_effect<-anpp_N_PP2%>%
  filter(nitrogen==0)%>%
  spread(treatment,anpp)
str(del_effect)
head(del_effect)

#create a treatment effect column of deluge - control (N only)
del_effect$deluge_effect_anpp<-(del_effect$PP2-del_effect$N)
head(del_effect)


#3. Calculate the predicted nitrogen + deluge treatment effect
#Drop unneeded cols so the merge function isn't confused
head(del_effect)
del_effect <- subset(del_effect, select = -c(nitrogen,N,PP2))
n_del_effect <- merge(n_effect, del_effect, by=c("block","site", "year"))
head(n_del_effect)

#Calculate the predictive additive effect of nitrogen + deluge
n_del_effect$predicted_n_del_anpp<-(n_del_effect$n_effect_anpp+n_del_effect$deluge_effect_anpp)


#4. Calculate production attributed to synergy
#filter for only PP2 treatments
  #exclude the 0g/m2 control since that was used to calculate the deluge effect and will give no new info
anpp_PP2 <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X,-plot,-grass,-forb,-woody,-dead,-anpp_with_woody,-cactus)%>%
  filter(site=="SGS",treatment=="PP2",year==2022, nitrogen!="0")%>%
  dplyr::rename("observed_n_del_anpp"="anpp")

#merge predicted ANPP to actual PP2 ANPP
n_del_effect_actual <- merge(n_del_effect, anpp_PP2, by=c("block","site", "year","nitrogen"))
head(n_del_effect_actual)

#create "synergy" anpp column of (observed_n_del_anpp - n_effect_anpp - deluge_effect_anpp)
n_del_effect_actual$Synergy<-(n_del_effect_actual$observed_n_del_anpp-n_del_effect_actual$n_effect_anpp-n_del_effect_actual$deluge_effect_anpp)
head(n_del_effect_actual)


#5. Determine significant differences between predicted and observed ANPP values
#create a new pared down dataset
pred_obs<-dplyr::select(n_del_effect_actual,c(block,nitrogen,n_effect_anpp,deluge_effect_anpp,Synergy))

#Calculate predicted and observed anpp totals
pred_obs$predicted=(pred_obs$n_effect_anpp+pred_obs$deluge_effect_anpp)
pred_obs$observed=(pred_obs$n_effect_anpp+pred_obs$deluge_effect_anpp+pred_obs$Synergy)

#Switch to long form
pred_obs_test <- gather(pred_obs, type, anpp, c(predicted,observed))

#Perform a t-test for main effects
  #data is normally distributed and have equal variance, so a pooled two-sample ttest is sufficient for all nitrogen levels
boxplot(anpp~type, data= pred_obs_test, ylab = "anpp increase")
t.test(anpp~type, var.equal=TRUE, data= pred_obs_test)

#Calculate the magnitude increase of ANPP between control and deluge plots
ddply(subset(pred_obs_test), c("type"), summarize, 
      n= length (anpp),
      mean = mean(anpp),
      sd = sd(anpp),
      se = sd/(sqrt(n)))


#6. Set-up data for plotting main effects
#Switch to long format
plot_data <- gather(n_del_effect_actual, contribution_by, anpp_increase, c(n_effect_anpp,deluge_effect_anpp,Synergy))
plot_data$contribution_by[plot_data$contribution_by == "n_effect_anpp"]<-"Nitrogen"
plot_data$contribution_by[plot_data$contribution_by == "deluge_effect_anpp"]<-"Deluge"

#Create two datasets for predicted and observed response
plot_data_pred<-subset(plot_data, plot_data$contribution_by!="Synergy")
plot_data_obs<-plot_data

#Create new columns for main effect group information
plot_data_pred$main<-"Pred."
plot_data_obs$main<-"Obs."

#Bind datasets back together
plot_data_main<-rbind(plot_data_pred,plot_data_obs)

#Create a lil dataset for later
sig_main<-plot_data_main %>%
  ungroup() %>% 
  dplyr::select(contribution_by,main) %>%
  dplyr::filter(contribution_by=="Synergy") %>%
  unique()


#7. Plot main effect plot
data_func<-barGraphStats(data=subset(plot_data_main), variable="anpp_increase", byFactorNames=c("contribution_by","main"))
data_func2 <- data_func %>%
  dplyr::select(-sd, -se) %>%
  tidyr::pivot_wider(names_from = contribution_by, 
                     values_from = mean) %>%
  dplyr::mutate(Deluge = Deluge + Nitrogen) %>%
  dplyr::mutate(Synergy = Synergy + Deluge) %>%
  tidyr::pivot_longer(c(Nitrogen, Deluge, Synergy), 
                      names_to = "contribution_by", 
                      values_to = "mean_alt")

data_func2 <- data_func2[complete.cases(data_func2), ] #remove the synergy row from the predicted group

data_func_final <- dplyr::left_join(data_func, data_func2) %>%
  dplyr::mutate(contribution_by = factor(contribution_by, 
                                         levels = c("Synergy","Deluge", "Nitrogen")))

data_func_final$main <- factor(data_func_final$main, levels = c("Pred.","Obs."))

#Remove se bars from the nitrogen and deluge sections since those aren't calculated for the observed data
data_func_final$se<-ifelse(data_func_final$contribution_by!="Synergy"&data_func_final$main=="Obs.",0,data_func_final$se)
data_func_final$mean_alt<-ifelse(data_func_final$contribution_by!="Synergy"&data_func_final$main=="Obs.",0,data_func_final$mean_alt)

pred_obs_effect_plot_main<-ggplot(data=subset(data_func_final), aes(x=main, y=mean, fill=contribution_by)) + 
  geom_bar(stat="identity", position = "stack", width = 0.6) +
  geom_errorbar(aes(ymin=mean_alt-se, ymax=mean_alt+se), width=0.4,linewidth=.5) +
  scale_y_continuous(limits = c(NA, 123), breaks = seq(0, max(pred_obs_test$anpp, na.rm = TRUE), by = 40)) +
  annotate("text", x=2, y=100, label ="*", size= 12, color="darkorange3") +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        #panel.grid.major = element_line(colour = "gray", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks = element_line(color="black",size=1.5), 
        axis.text.x = element_text(color="black", size=18),
        axis.text.y = element_text(color="black",size=15),
        plot.title = element_text(size=25), #set as 25
        axis.title.x = element_text(color="white",size=20, vjust = .4),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        legend.position=c(.45,.45),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1)) +
  scale_fill_manual(values=c("gray50","#87bcbd","#434475"),name="Contribution by")+ #alternate palette
  geom_vline(aes(xintercept=-Inf), size=1.5) + 
  geom_hline(aes(yintercept=0), size=1, linetype="dashed") + 
  labs(y=expression(ANPP~Response~(g~m^-2))) + #y label
  labs(x=expression(Nitrogen~(g~m^-2)))+
  ggtitle("a.")


#8. Plot predicted additive vs observed synergistic plot
  #Switch to long form
pred_obs_test <- gather(pred_obs, type, anpp, c(predicted,observed))
pred_obs_test$type[pred_obs_test$type == "observed"]<-"Observed Synergistic"
pred_obs_test$type[pred_obs_test$type == "predicted"]<-"Predicted Additive"

str(pred_obs_test)
pred_obs_test$nitrogen <- as.numeric(pred_obs_test$nitrogen)

pred_obs_effect_plot<-ggplot(data=pred_obs_test, aes(x=nitrogen, y=anpp, fill=type, color=type)) + 
  geom_smooth(method="lm", formula = y~x+I(x^2), se=TRUE) +
  #geom_bar(stat="identity", position = "stack", width = 0.6) +
  #geom_errorbar(aes(ymin=mean_alt-se, ymax=mean_alt+se), width=0.4,linewidth=.5) +
  scale_y_continuous(limits = c(NA,NA), breaks = seq(0, max(pred_obs_test$anpp, na.rm = TRUE), by = 40)) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        #panel.grid.major = element_line(colour = "gray", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks = element_line(color="black", size=1.5), 
        axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black",size=15),
        plot.title = element_text(size=25), #set as 25
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        strip.text.x = element_text(size = 18),
        legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(size = 16, hjust = 0.5),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1)) +
  scale_fill_manual(values=c("gray50","#87bcbd","#434475"), labels=function(x)str_wrap(x, width = 15))+ #alternate palette
  scale_color_manual(values=c("gray50","#87bcbd","#434475"), labels=function(x)str_wrap(x, width = 15)) +
  geom_vline(aes(xintercept=-Inf), size=1.5) + 
  geom_hline(aes(yintercept=0), size=1, linetype="dashed") + 
  labs(y=expression(Predicted~ANPP~Response~(g~m^-2))) + #y label
  labs(x=expression(Nitrogen~(g~m^-2)))+
  ggtitle("b.")

ggarrange(pred_obs_effect_plot_main,pred_obs_effect_plot,
                                ncol = 2, nrow = 1,
                                widths = c(1,2))

######

####Richness: Analysis, Visualization, and ANPP correlation
######

#Load in
sppcompdata_full <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")

#Reserves nitrogen treatment information--this will need to be merged back in after calculations
plotlist_trts<-sppcompdata_full %>%
  dplyr::select(block, plot, nitrogen) %>%
  unique()

#Calculate evar and richness
#This uses absolute cover, even though having cover <100 is not ideal (since we didn't collect a "bareground" metric)
sppcompdata_full <- sppcompdata_full%>% #create replicate var to distinguish two sites
  mutate(rep=paste(site, plot, treatment, sep='::'))
sppchange_noN<-community_structure(sppcompdata_full, time.var = "year", abundance.var = "cover", replicate.var = "rep",
                                   metric = c("Evar"))%>%
  separate(rep, c('site', 'plot',"treatment"), sep='::') #splits back apart the rep into site and plot
head(sppchange_noN) #Some evenness values contain NAs because there are plots with only one species

#Merge sppchange_noN dataset to nitrogen treatments
sppchange <- merge(sppchange_noN, plotlist_trts, by=c("plot"))

#correct facet wrap labels for consistency
sppchange$treatment[sppchange$treatment == "nitrogen"]<-"Nitrogen"
sppchange$treatment[sppchange$treatment == "pulse press 1"]<-"N+Deluge"
sppchange$treatment[sppchange$treatment == "pulse press 2"]<-"N+Deluge"

#Fix structure
sppchange_n<-sppchange %>% 
  dplyr::mutate_at(c("block","plot","site","nitrogen","year","treatment"),as.factor)
str(sppchange_n) 

#Order factors
sppchange_n$treatment <- factor(sppchange_n$treatment, levels=c("Nitrogen", "N+Deluge"))

#Replicated mixed-effect model with a split-plot design and blocking, followed by type III ANOVA
  #Fixed effects: Nitrogen, Deluge
  #Random effects: Block, Plot
  #Note: We are only analyzing data within each year to account for the high MAP variability
  #Note: No unique "plotID" (i.e., block*plot) is needed since plots are already unique values (1-48), this unique term captures repeated measures

richness_shortgrass_byyear <- lmer(richness~nitrogen*treatment #nitrogen and deluge treatment as fixed
                                   + (1|block) + (1|plot), #block and plot as random, (1|plot) reflects the split-plot design
                                   data=subset(sppchange_n, site=='sgs'&year==2022)) #subset CPER in a single year

check_model(richness_shortgrass_byyear)
anova(richness_shortgrass_byyear, ddf="Kenward-Roger") #ddf="Kenward-Roger" to calculate Type III test for mixed models

#Pairwise comparisons
em_richness_shortgrass_byyear<-emmeans(richness_shortgrass_byyear,pairwise~nitrogen)
cld(em_richness_shortgrass_byyear$emmeans)

#Create a data frame for the cld
cld_richness<-sppchange_n %>%
  filter(site=="sgs")%>%
  ungroup() %>% 
  dplyr::select(site,nitrogen,year) %>%
  unique()%>%
  arrange(year, nitrogen)

#manually add the cld letters to the data frame 
cld_richness$cld_richness<-c("ab","a","ab","ab","ab","ab","b","ab",  #2021
                             "a","ab","b","ab","ab","ab","b","b") #2022

sppchange_n <- merge(sppchange_n, cld_richness, by=c("site","nitrogen","year"))

#plot
rich_a<-ggplot(data=barGraphStats(data=subset(sppchange_n, site=='sgs'), variable="richness", byFactorNames=c("nitrogen","year","cld_richness")), aes(x=nitrogen, y=mean, color=nitrogen)) + 
  geom_point(size=5) +
  facet_grid(year~.) +
  geom_text(aes(label = cld_richness, y = mean+se), vjust=-.5, size=4, color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15, vjust = .4),
        axis.text.y=element_text(color="black",size=15),
        axis.title=element_text(color="black",size=20),
        plot.title = element_text(size=25),
        strip.text.x = element_text(color="black",size = 20),
        strip.text.y = element_blank(),
        legend.position="none") +
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray60","#6f9954", "#89ab7c", "#87bcbd", "#969bc7", "#6b6ca3", "#434475", "#2c2d54"))+ #alternate palette
  #scale_color_manual(values=met.brewer("VanGogh1", 8, direction=-1)) +
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  labs(y=expression(Richness)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ylim(0,11) +
  ggtitle("a.")


#Averaged across N treatments
#rename for plotting
sppchange$treatment[sppchange$treatment == "Nitrogen"]<-"N"
sppchange$treatment[sppchange$treatment == "N+Deluge"]<-"N+D"

#create a new cld
cld_richness_trt<-sppchange %>%
  filter(site=="sgs")%>%
  ungroup() %>% 
  dplyr::select(site,treatment,year) %>%
  unique()%>%
  arrange(year,treatment)

cld_richness_trt$cld_richness_trt<-c("","","a","b")
sppchange_trt <- merge(sppchange, cld_richness_trt, by=c("site","treatment","year"))

rich_b<-ggplot(data=barGraphStats(data=subset(sppchange_trt, site=='sgs'), variable="richness", byFactorNames=c("year","cld_richness_trt","treatment")), aes(x=treatment, y=mean, fill=treatment)) + 
  geom_bar(stat = "identity") +
  facet_grid(year~.) +
  geom_text(aes(label = cld_richness_trt, y = mean+se), vjust=-.5, size=4, color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15, vjust = .4),
        axis.text.y=element_text(color="black",size=15),
        axis.title.x=element_text(color="white",size=20),
        axis.title.y=element_blank(),
        plot.title = element_text(size=25),
        strip.text.x = element_text(color="black",size = 20),
        strip.text.y = element_text(color="black",size = 20),
        legend.position="none") +
  scale_fill_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray60","#87bcbd")) +
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  labs(y=expression(Richness)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ylim(0,8) +
  ggtitle("b.")

rich_ab<-ggarrange(rich_a, rich_b, ncol=2,nrow=1, widths = c(1,.6))
rich_ab


###Correlation of richness to anpp
#1. Load in production data
anpp_all <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X)

#Rename treatments for merging ease
anpp_all$treatment[anpp_all$treatment == "N"]<-"N"
anpp_all$treatment[anpp_all$treatment == "PP1"]<-"N+D"
anpp_all$treatment[anpp_all$treatment == "PP2"]<-"N+D"

#Reload the species change data
sppchange <- merge(sppchange_noN, trts, by=c("plot"))

#Rename treatments for merging ease
sppchange$treatment[sppchange$treatment == "nitrogen"]<-"N"
sppchange$treatment[sppchange$treatment == "pulse press 1"]<-"N+D"
sppchange$treatment[sppchange$treatment == "pulse press 2"]<-"N+D"

#merge species change and production datasets
head(anpp_all)
str(anpp_all)
head(sppchange)
str(sppchange)

sppchange_anpp <- merge(sppchange, anpp_all, by=c("plot","year","block","treatment","nitrogen"))

#Fix structure
sppchange_anpp<-sppchange_anpp %>% 
  dplyr::mutate_at(c("block","plot","nitrogen","year","treatment"),as.factor)
str(sppchange_anpp) 

#Order factors
sppchange_anpp$treatment <- factor(sppchange_anpp$treatment, levels=c("N", "N+D"))

#plot
ggplot(data=subset(sppchange_anpp),  aes(x=richness, y=anpp, color=nitrogen, group=year, linetype = treatment, shape = treatment)) +
  geom_point(size=4) +
  stat_cor(label.y = 215) +
  geom_smooth(aes(group=treatment),method = "lm", size=1.5, color = "darkorange3") +
  facet_grid(year~., scales = "free") +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15),
        axis.text.y=element_text(color="black",size=15),
        axis.title=element_text(color="black",size=20),
        plot.title = element_text(color = "black", size = 25, hjust = 0),
        panel.grid.minor = element_line(colour = NA), 
        axis.title.y = element_text(vjust=1.2),
        strip.text.y = element_text(size = 20),
        legend.position = c(.3,.9),
        #legend.direction= "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1)) +
  geom_hline(yintercept = 0, size = 1, linetype="dashed") +
  geom_hline(yintercept = -Inf, size = 1.5) +
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray60","#6f9954","#89ab7c","#87bcbd","#969bc7","#6b6ca3","#434475","#2c2d54"), guide="none") +
  labs(y=expression(ANPP~(g~m^-2)))+ #y label
  labs(x=expression(Richness)) + #x label
  ggtitle("C.")


#simple analysis
test_2021<-filter(sppchange_anpp,year==2021)
test_2022<-filter(sppchange_anpp,year==2022)
test_2021N<-filter(sppchange_anpp,year==2021&treatment=="N")
test_2021ND<-filter(sppchange_anpp,year==2021&treatment=="N+D")
test_2022N<-filter(sppchange_anpp,year==2022&treatment=="N")
test_2022ND<-filter(sppchange_anpp,year==2022&treatment=="N+D")

#Run analysis on dataset of interest
cor.test(test_2022ND$richness,test_2022ND$anpp)
model1<-lm(formula = anpp ~ richness, data = test_2022ND)
summary(model1)

######

####Evar: Analysis ANOVA and Visualization
######

#Load in dataset
sppcompdata_full <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")

#Reserves nitrogen treatment information--this will need to be merged back in after calculations
plotlist_trts<-sppcompdata_full %>%
  dplyr::select(block, plot, nitrogen) %>%
  unique()

#Calculate evar and richness
#This uses absolute cover, even though having cover <100 is not ideal (since we didn't collect a "bareground" metric)
sppcompdata_full <- sppcompdata_full%>% #create replicate var to distinguish two sites
  mutate(rep=paste(site, plot, treatment, sep='::'))
sppchange_noN<-community_structure(sppcompdata_full, time.var = "year", abundance.var = "cover", replicate.var = "rep",
                                   metric = c("Evar"))%>%
  separate(rep, c('site', 'plot',"treatment"), sep='::') #splits back apart the rep into site and plot
head(sppchange_noN) #Some evenness values contain NAs because there are plots with only one species

#Merge sppchange_noN dataset to nitrogen treatments
sppchange <- merge(sppchange_noN, plotlist_trts, by=c("plot"))

#correct facet wrap labels for consistency
sppchange$treatment[sppchange$treatment == "nitrogen"]<-"Nitrogen"
sppchange$treatment[sppchange$treatment == "pulse press 1"]<-"N + Deluge"
sppchange$treatment[sppchange$treatment == "pulse press 2"]<-"N + Deluge"

#Fix structure
sppchange<-sppchange %>% 
  dplyr::mutate_at(c("block","plot","site","nitrogen","year","treatment"),as.factor)
str(sppchange) 

#Order factors
sppchange$treatment <- factor(sppchange$treatment, levels=c("Nitrogen", "N + Deluge"))

#Replicated mixed-effect model with a split-plot design and blocking, followed by type III ANOVA
  #Fixed effects: Nitrogen, Deluge
  #Random effects: Block, Plot
  #Note: We are only analyzing data within each year to account for the high MAP variability
  #Note: No unique "plotID" (i.e., block*plot) is needed since plots are already unique values (1-48), this unique term captures repeated measures

evar_shortgrass_byyear <- lmer(Evar~nitrogen*treatment #nitrogen and deluge treatment as fixed
                               + (1|block) + (1|plot), #block and plot as random, (1|plot) reflects the split-plot design
                               data=subset(sppchange, site=='sgs'&year==2022)) #subset CPER in a single year

check_model(evar_shortgrass_byyear)
anova(evar_shortgrass_byyear, ddf="Kenward-Roger") #ddf="Kenward-Roger" to calculate Type III test for mixed models

#Pairwise comparisons
em_evar_shortgrass_byyear<-emmeans(evar_shortgrass_byyear,~nitrogen)
cld(em_evar_shortgrass_byyear)

#Create a data frame for the cld
cld_evar<-sppchange %>%
  filter(site=="sgs")%>%
  ungroup() %>% 
  dplyr::select(site,treatment,year) %>%
  unique()%>%
  arrange(year,treatment)

#manually add the cld letters to the data frame 
cld_evar$cld_evar<-c("","", #2021
                     "a","b") #2022

sppchange_evar <- merge(sppchange, cld_evar, by=c("site","treatment","year"))

#correct facet wrap labels for consistency
sppchange_evar$treatment<-as.character(sppchange_evar$treatment)
sppchange_evar$treatment[sppchange_evar$treatment == "Nitrogen"]<-"N"
sppchange_evar$treatment[sppchange_evar$treatment == "N + Deluge"]<-"N+D"
sppchange_evar$treatment<-factor(sppchange_evar$treatment)

#plot
evar_a<-ggplot(data=barGraphStats(data=subset(sppchange_n, site=='sgs'), variable="Evar", byFactorNames=c("nitrogen","year")), aes(x=nitrogen, y=mean, color=nitrogen)) + 
  geom_point(size=5) +
  facet_grid(year~.) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15, vjust = .4),
        axis.text.y=element_text(color="black",size=15),
        axis.title=element_text(color="black",size=20),
        plot.title = element_text(size=25),
        strip.text.x = element_text(color="black",size = 20),
        strip.text.y = element_blank(),
        legend.position="none") +
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray60","#6f9954", "#89ab7c", "#87bcbd", "#969bc7", "#6b6ca3", "#434475", "#2c2d54"))+ #alternate palette
  #scale_color_manual(values=met.brewer("VanGogh1", 8, direction=-1)) +
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  labs(y=expression(Evar)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ylim(0,.6) +
  ggtitle("a.")

evar_b<-ggplot(data=barGraphStats(data=subset(sppchange_evar, site=='sgs'), variable="Evar", byFactorNames=c("year","cld_evar","treatment")), aes(x=treatment, y=mean, fill=treatment)) + 
  geom_bar(stat="identity") +
  facet_grid(year~.) +
  geom_text(aes(label = cld_evar, y = mean+se), vjust=-.5, size=4, color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15),
        axis.text.y=element_text(color="black",size=15),
        axis.title.y=element_blank(),
        axis.title.x = element_text(color="white",size=20, vjust = .4),
        plot.title = element_text(size=25),
        strip.text.x = element_text(color="black",size = 20),
        strip.text.y = element_text(color="black",size = 20),
        legend.position="none") +
  scale_fill_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray60","#87bcbd"))+ #alternate palette
  #scale_color_manual(values=met.brewer("VanGogh1", 8, direction=-1)) +
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  labs(y=expression(Evar)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ylim(0,0.6) +
  ggtitle("b.")

evar<-ggarrange(evar_a, evar_b, ncol=2,nrow=1, widths = c(1,.6))
evar

#Calculate magnitude of change
ddply(subset(sppchange, site=='sgs'), c("year","treatment"), summarize, 
      n= length (Evar),
      mean = mean(Evar),
      sd = sd(Evar),
      se = sd/(sqrt(n)))

######




####Community composition: Analysis and Visualization
######

#Comparison across between plot and deluge subplot: does nitrogen affect how community composition responds to deluge in 2021 or 2022?
  #Adonis() can't handle the split plot design.
  #Instead, we can calculate distance between the centroids (using codyn) of the nitrogen plots and the corresponding pulse-press plots (so 1 value per plot)
  #We can then use our standard mixed model with the distance values as the response metric

#Load in an process each year of data separately
sppcompsgs <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(year==2021, site=="sgs")

sppcompsgs <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(year==2022, site=="sgs")

#Reserves nitrogen treatment information--this will need to be merged back in after calculations
plotlist_trts<-sppcompsgs %>%
  dplyr::select(block, plot, nitrogen) %>%
  unique()

#multivariate_difference: "Calculates the changes in composition/dispersion between treatments based off a Bray-Curtis dissimilarity matrix"
  #"Composition change is the pairwise distance between centroids of compared treatments and ranges from 0-1" 0 = identical, 1 = completely different communities
  #"Dispersion change is the difference between treatments in the dispersion of replicates, i.e. the average distance between a replicate and its centroid"
    #Dispersion change can be calculated without multiple replicates
#merge "treatment" and "plot" col to create a unique value
sppcompsgs$trt_plot<-paste(sppcompsgs$treatment, sppcompsgs$plot, sep="_")

#calculate composition_change
  #use "plot" as the "treatment variable"--this forces the function to calculate values for each plot
  #can also use the function multivariate_difference(), it calculates the same values but requires different inputs
comp_diff<-multivariate_change(
  sppcompsgs,
  species.var = "species",
  abundance.var = "cover",
  replicate.var = "trt_plot",
  treatment.var = "plot",
  time.var = "treatment"
)

#merge back in nitrogen treatment data:
comp_diff_n <- merge(comp_diff, plotlist_trts, "plot")

#Convert factors
comp_diff_n<-comp_diff_n %>% 
  dplyr::mutate_at(c("block","plot","nitrogen"),as.factor)
str(comp_diff_n)


#Replicated mixed-effect model with a split-plot design and blocking, followed by type III ANOVA
  #Fixed effects: Nitrogen, Deluge
  #Random effects: Block, Plot
  #Note: We are only analyzing data within each year to account for the high MAP variability
  #Note: No unique "plotID" (i.e., block*plot) is needed since plots are already unique values (1-48), this unique term captures repeated measures

#CPER
comp_shortgrass_byyear <- lmer(composition_change~nitrogen #nitrogen and deluge treatment as fixed
                               + (1|block), #block and plot as random, (1|plot) reflects the split-plot design
                               data=subset(comp_diff_n)) #subset CPER in a single year

check_model(comp_shortgrass_byyear)
anova(comp_shortgrass_byyear, ddf="Kenward-Roger") #ddf="Kenward-Roger" to calculate Type III test for mixed models

#Pairwise comparisons
em_shortgrass_nitrogen<-emmeans(comp_shortgrass_byyear,~nitrogen)
pairs(em_shortgrass_nitrogen)

#Create a dataframe for the cld
cld_comp<-comp_diff_n %>%
  ungroup() %>% 
  dplyr::select(nitrogen,treatment) %>%
  unique()%>%
  arrange(treatment, nitrogen)

#manually add the cld letters to the data frame 
cld_comp$cld_comp<-c("a","a","a","a","a","a","a","a")

comp_diff_n <- merge(comp_diff_n, cld_comp, by=c("nitrogen","treatment"))

#Production across N treatments faceted by year
n_pp1_2021<-ggplot(data=barGraphStats(data=subset(comp_diff_n), variable="composition_change", byFactorNames=c("nitrogen","cld_comp")), aes(x=nitrogen, y=mean, color=nitrogen)) + 
  geom_point(size=5) +
  facet_grid("2021"~.) +
  #geom_text(aes(label = cld_comp, y = mean+se), vjust=-.5, size=4, color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        #axis.text.x=element_text(color="black", size=12, angle = 90),
        axis.text.x=element_blank(),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(color="black",size=20),
        axis.title = element_blank(),
        #axis.title.y = element_text(color="black",size=20,vjust=1.2),
        strip.text.y = element_text(color="black",size = 20),
        legend.position="none") +
  scale_color_manual(values=c("gray60","#6f9954", "#89ab7c", "#87bcbd", "#969bc7", "#6b6ca3", "#434475", "#2c2d54"))+ #alternate palette
  #scale_color_manual(values=met.brewer("VanGogh1", 8, direction=-1)) +
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  labs(y=expression(Community~Difference)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ylim(0,0.7) +
  ggtitle("a.")

n_pp2_2022<-ggplot(data=barGraphStats(data=subset(comp_diff_n), variable="composition_change", byFactorNames=c("nitrogen","cld_comp")), aes(x=nitrogen, y=mean, color=nitrogen)) + 
  geom_point(size=5) +
  facet_grid("2022"~.) +
  #geom_text(aes(label = cld_comp, y = mean+se), vjust=-.5, size=4, color="black") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.5) +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white", size = 0.5),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=12, angle = 90, vjust = .5),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_blank(),
        axis.title.x  = element_text(color="black",size=20,vjust=1.2),
        axis.title.y = element_blank(),
        strip.text.y = element_text(color="black",size = 20),
        legend.position="none") +
  scale_color_manual(values=c("gray60","#6f9954", "#89ab7c", "#87bcbd", "#969bc7", "#6b6ca3", "#434475", "#2c2d54"))+ #alternate palette
  #scale_color_manual(values=met.brewer("VanGogh1", 8, direction=-1)) +
  geom_vline(aes(xintercept=-Inf), linewidth=2) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  labs(y=expression(Community~Difference)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) + #x label
  ylim(0,0.7) 

temp<-ggarrange(n_pp1_2021,n_pp2_2022,
          ncol = 1, nrow = 2,
          heights = c(1,1.15))

annotate_figure(temp,
                left = text_grob("Centroid Difference",color="black",size=20, rot = 90))

######

####Community composition: Simper to identify influential species
######

#Load in and assess each year of data separately
sppcomp_sgs <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(year==2021,site=="sgs")

sppcomp_sgs <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(year==2022,site=="sgs")

#Change the nitrogen treatment name, so the wide format isn't confused
sppcomp_sgs$treatment[sppcomp_sgs$treatment == "nitrogen"]<-"press"

#Adjust structure
sppcomp_sgs<-sppcomp_sgs %>% 
  dplyr::mutate_at(c("site","year","plot","nitrogen","block","treatment"),as.factor)
str(sppcomp_sgs)

#Convert to wide format and create matrix
compwide_sgs <- spread(sppcomp_sgs, species, cover, fill = 0)
maxcol<-ncol(compwide_sgs)
matrix<-(compwide_sgs[,7:maxcol]) #Columns that contain the species names

#Understand community change through Simper analysis
sim_SGS<-with(compwide_sgs,simper(matrix,nitrogen)) #Currently assessing deluge treatment, could switch to nitrogen for slightly diff. results, but species and order of importance is not impacted

#Identify the top three/four species that contributed most to overall dissimilarity
#Summarize overall dissimilarity
species<-index$species
average_column=c()
species_column=c()
contrast_column=c()
contrast_names=names(summary(sim_SGS))

for(contrast in contrast_names){
  index<-as.data.frame(do.call(cbind, sim_SGS[[contrast]]))
  index$average <- as.numeric(index$average)
  average_column=append(average_column,index$average)
  species_column=append(species_column,index$species)
  contrast_column=append(contrast_column, rep(contrast,33)) #number of species
}

species_avg_SGS=data.frame(species=species_column,average=average_column,contrast=contrast_column)

dissim_avg_SGS<-ddply(species_avg_SGS, c("species"), summarize,
                      n= length (average),
                      mean = mean(average),
                      sd = sd(average),
                      se = sd/(sqrt(n)))
#Results based on overall dissimilarity

#Create a dataframe with the names of the top contributing species for easier graphing later
top_4_sgs_2021_del<-c("Elymus elymoides", "Vulpia octoflora", "Chenopodium sp.", "Bromus tectorum")
top_4_sgs_2022_del<-c("Elymus elymoides", "Chenopodium sp.", "Salsola tragus", "Sphaeralcea coccinea")
top_3_sgs_2021_del<-c("Elymus elymoides", "Vulpia octoflora", "Bromus tectorum")
top_3_sgs_2022_del<-c("Elymus elymoides", "Chenopodium sp.", "Salsola tragus")

######

####Community composition: Analysis + Visualization of influential species
######

#See code above to see how these were determined
top_4_sgs_2021_del<-c("Elymus elymoides", "Vulpia octoflora", "Chenopodium sp.", "Bromus tectorum")
top_4_sgs_2022_del<-c("Elymus elymoides", "Chenopodium sp.", "Salsola tragus", "Sphaeralcea coccinea")
top_3_sgs_2021_del<-c("Elymus elymoides", "Vulpia octoflora", "Bromus tectorum")
top_3_sgs_2022_del<-c("Elymus elymoides", "Chenopodium sp.", "Salsola tragus")

#Load in data
sppcompdata_PP1 <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(year==2021,site=="sgs", species%in%top_3_sgs_2021_del)

sppcompdata_PP2 <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(year==2022,site=="sgs", species%in%top_3_sgs_2022_del)

sppcomp_sgs<-full_join(sppcompdata_PP1,sppcompdata_PP2)

#Change the nitrogen treatment name, so the wide format isn't confused
sppcomp_sgs$treatment[sppcomp_sgs$treatment == "nitrogen"]<-"Nitrogen"
sppcomp_sgs$treatment[sppcomp_sgs$treatment == "pulse press 1"]<-"N+Deluge"
sppcomp_sgs$treatment[sppcomp_sgs$treatment == "pulse press 2"]<-"N+Deluge"


#Analyze most-contributing species
#Replicated mixed-effect model with a split-plot design and blocking, followed by type III ANOVA
  #Fixed effects: Nitrogen, Deluge
  #Random effects: Block, Plot
  #Note: We are only analyzing data within each year to account for the high MAP variability
  #Note: No unique "plotID" (i.e., block*plot) is needed since plots are already unique values (1-48), this unique term captures repeated measures
anpp_byyear <- lmer(cover~nitrogen*treatment #nitrogen and deluge treatment as fixed
                    + (1|block) + (1|plot), #block and plot as random, (1|plot) reflects the split-plot design
                    data=subset(sppcomp_sgs, year==2022&species=="Salsola tragus")) #subset CPER in a single year

check_model(anpp_byyear)
anova(anpp_byyear, ddf="Kenward-Roger") #ddf="Kenward-Roger" to calculate Type III test for mixed models

#Pairwise comparisons
em_shortgrass_nitrogen_by_treatment<-emmeans(anpp_byyear,pairwise~nitrogen|treatment)
pairs(em_shortgrass_nitrogen_by_treatment)

#Calculate the magnitude increase of ANPP between control and deluge plots
sppcomp_sgs$fertilized<-ifelse(sppcomp_sgs$nitrogen=="0","control","nitrogen")
ddply(subset(sppcomp_sgs, year==2021&species=="Bromus tectorum"), c("treatment"), summarize, 
      n= length (cover),
      mean = mean(cover),
      sd = sd(cover),
      se = sd/(sqrt(n)))

#prep for plotting
sppcomp_sgs_temp<-sppcomp_sgs
sppcomp_sgs_temp$treatment[sppcomp_sgs_temp$treatment == "Nitrogen"]<-"N"
sppcomp_sgs_temp$treatment[sppcomp_sgs_temp$treatment == "N+Deluge"]<-"N+D"

#Order factored treatments
sppcomp_sgs$nitrogen<-as.factor(sppcomp_sgs$nitrogen)
sppcomp_sgs$treatment <- factor(sppcomp_sgs$treatment, levels=c("Nitrogen", "N+Deluge"))

#Nitrogen effect
abun_a<-ggplot(data=barGraphStats(data=subset(sppcomp_sgs), variable="cover", byFactorNames=c("species","year","nitrogen")), aes(x=nitrogen, y=mean, color=species, group=species)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  geom_point(size=5) +
  geom_line(size=1) +
  facet_grid(year~., scales="free")+
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=12, angle = 90, vjust = 0.5),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(color = "black", size = 25, hjust = 0),
        axis.title.x = element_text(color="black",size=15),
        axis.title.y = element_text(color="black",size=20),
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=15),
        strip.text = element_blank()) +
  geom_hline(aes(yintercept=-Inf), size=2) +
  geom_hline(aes(yintercept=0), size=1, linetype="dashed") + 
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray50","darkorange3","#6f9954","#2c2d54","#969bc7"))+ #alternate palette
  labs(y=expression(Abundance)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) +
  ggtitle("C.")

#plot deluge effect
abun_b<-ggplot(data=barGraphStats(data=subset(sppcomp_sgs_temp), variable="cover", byFactorNames=c("species","year","treatment")), aes(x=treatment, y=mean, color=species, group=species)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) +
  geom_point(size=5) +
  geom_line(size=1) +
  facet_grid(year~.)+
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(color = "black", size = 25, hjust = 0),
        axis.title.x = element_text(color="white",size=15),
        axis.title.y = element_blank(),
        legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=12, face = "italic"),
        strip.text = element_text(color="black",size = 20)) +
  geom_hline(aes(yintercept=0), size=1, linetype="dashed") + 
  geom_hline(aes(yintercept=-Inf), size=2) +
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray50","darkorange3","#6f9954","#2c2d54","#969bc7"))+ #alternate palette
  labs(y=expression(Abundance)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) +
  ggtitle("") +
  guides(color = guide_legend(label.theme = element_text(size = 12, lineheight = 1.2), label.position = "top",  title.position = "top"))

ggarrange(abun_a,abun_b, ncol=2, nrow=1, widths = c(.65,1))

######

####Multiple regression model: predictors of ANPP for 2022 deluge data only
######

#Use model selection and multiple regression to understand the best predictors of anpp

#1. Create one large dataset with all potential predictors
  #predictors must be univariate in a single column
  #not currently including nitrogen and deluge explicitly in the model because we don't want to overlap treatments with associated responses
  #anpp predictors: deluge treatment (represented as soil moisture), nitrogen (represented as nitrate/total N/ammonium), chenopodium cover/proportion forbs/evenness, green-up

#a. Nitrogen availability dataframe
nitrogen_x <- read.csv("MAIN_ChANGE_PulsePress_SGS_corrected ammonium nitrate_2021to2022.csv")

#remove outlier
#one nitrate value was off the charts without explanations: sgs, C23, 2021, 7.5, nitrogen plot
#removed all 2021 C23 data, so calculations aren't effected (total N and percent change)
nitrogen<-nitrogen_x %>%
  dplyr::filter(outlier!="x") %>%
  dplyr::select(-outlier) %>%
  dplyr::rename(year=sample_year)

#separate block_plot col into two
nitrogen<-tidyr::separate(nitrogen, col=plot, into=c('block', 'plot'), sep=1)

#rename treatments and channels for merging/plotting ease
nitrogen$site[nitrogen$site == "sgs"]<-"SGS"
nitrogen$treatment[nitrogen$treatment == "n"]<-"Nitrogen"
nitrogen$treatment[nitrogen$treatment == "pp1"]<-"N+Deluge"
nitrogen$treatment[nitrogen$treatment == "pp2"]<-"N+Deluge"
nitrogen$channel[nitrogen$channel == "nitrate"]<-"Nitrate"
nitrogen$channel[nitrogen$channel == "ammonium"]<-"Ammonium"

#spread data, so we can calculate total N
nitrogen_aic<-tidyr::spread(nitrogen,channel,avgcalc_mg_l)
nitrogen_aic$"Total_N"<-(nitrogen_aic$Ammonium+nitrogen_aic$Nitrate)


#b. Soil moisture dataframe
vwc <- read.csv("MAIN_ChANGE_PulsePress_SGS_soil moisture_2021to2022.csv")%>%
  dplyr::select(-X)%>%
  dplyr::filter(site=="sgs",days_since_deluge>0&days_since_deluge<26)

#simplify deluge treatment name--will separate PP1 and PP2 treatments later by year
vwc$treatment[vwc$treatment == "nitrogen"]<-"Nitrogen"
vwc$treatment[vwc$treatment == "pulse-press1"]<-"N+Deluge"
vwc$treatment[vwc$treatment == "pulse-press2"]<-"N+Deluge"
vwc$site[vwc$site == "sgs"]<-"SGS"

#average vwc down to a single measure 30-day average for each plot/treatment
vwc_aic<-vwc %>% 
  dplyr::group_by(plot,block,site,treatment,year)%>% #remove "Remove" column by leaving it out of the function
  dplyr::summarise(across(vwc,~mean(.x,na.rm=TRUE))) %>% 
  dplyr::as_tibble()


#c. Green-up dataframe
green_2021 <- read.csv("MAIN_ChANGE_PulsePress_SGS_calculated green up values_2021to2022.csv")%>%
  dplyr::select(-X,-format)%>%
  dplyr::filter(year==2021)

green_2022 <- read.csv("MAIN_ChANGE_PulsePress_SGS_calculated green up values_2021to2022.csv")%>%
  dplyr::select(-X,-format)%>%
  dplyr::filter(treatment!="pp1",year==2022,doy<211&doy>179) #NOTE: filter to data from after the deluge to day 210

str(green_2022)
green<-rbind(green_2021,green_2022)

#Clarify factors
green$block<-toupper(green$block)
green$site[green$site == "sgs"]<-"SGS"
green$treatment[green$treatment == "n"]<-"Nitrogen"
green$treatment[green$treatment == "n"]<-"Nitrogen"
green$treatment[green$treatment == "pp1"]<-"N+Deluge"
green$treatment[green$treatment == "pp2"]<-"N+Deluge"

#average vwc down to a single measure 30-day average for each plot/treatment
green_aic<-green %>% 
  dplyr::group_by(plot,block,site,treatment,year,nitrogen)%>% #remove "Remove" column by leaving it out of the function
  dplyr::summarise(across(gcc_mean,~mean(.x,na.rm=TRUE))) %>% 
  dplyr::as_tibble()


#d. Chenopodium cover dataframe
sppchange <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  dplyr::filter(species=="Chenopodium sp.")%>%
  dplyr::rename(cheno_cover=cover)

#correct facet wrap labels for consistency
sppchange$site[sppchange$site == "sgs"]<-"SGS"
sppchange$treatment[sppchange$treatment == "nitrogen"]<-"Nitrogen"
sppchange$treatment[sppchange$treatment == "pulse press 1"]<-"N+Deluge"
sppchange$treatment[sppchange$treatment == "pulse press 2"]<-"N+Deluge"

sppchange_aic<-sppchange %>%
  dplyr::select(-species)


#e. Proportion of forbs dataset (includes anpp)
anpp_PP <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X,-woody,-dead,-anpp_with_woody,-cactus)

#rename treatments for clarity
anpp_PP$treatment[anpp_PP$treatment == "N"]<-"Nitrogen"
anpp_PP$treatment[anpp_PP$treatment == "PP1"]<-"N+Deluge"
anpp_PP$treatment[anpp_PP$treatment == "PP2"]<-"N+Deluge"

#create a "forb percentage" column by diving forbs by total anpp
anpp_PP$forb_prct<-(anpp_PP$forb/anpp_PP$anpp)

#remove grass and fob col for easier merge
anpp_forbs_aic<-anpp_PP %>%
  dplyr::select(-grass,-forb)


#f. Total cover dataframe
sppchange <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  dplyr::rename(total_cover=cover)

#correct factors for consistency
sppchange$site[sppchange$site == "sgs"]<-"SGS"
sppchange$treatment[sppchange$treatment == "nitrogen"]<-"Nitrogen"
sppchange$treatment[sppchange$treatment == "pulse press 1"]<-"N+Deluge"
sppchange$treatment[sppchange$treatment == "pulse press 2"]<-"N+Deluge"

#average vwc down to a single measure 30-day average for each plot/treatment
total_cover_aic<-sppchange %>% 
  dplyr::group_by(plot,block,site,treatment,year,nitrogen)%>% #remove "Remove" column by leaving it out of the function
  dplyr::summarise(across(total_cover,~mean(.x,na.rm=TRUE))) %>% 
  dplyr::as_tibble()


#2. Merge all datasets so each predictor is in its own col

#start with "anpp+proportion of forbs" dataset since it includes the response value (ANPP)
#add available N data
head(anpp_forbs_aic) #8cols
head(nitrogen_aic) #9cols
anpp_forb_nitrogen_aic<-merge(anpp_forbs_aic,nitrogen_aic,all.x=TRUE) #11cols, good

#add soil moisture
head(anpp_forb_nitrogen_aic) #11cols
head(vwc_aic) #should add 1col
anpp_forb_nitrogen_vwc_aic<-merge(anpp_forb_nitrogen_aic,vwc_aic,all.x=TRUE) #12cols, good

#add green-up
head(anpp_forb_nitrogen_vwc_aic) #12cols
head(green_aic) #should add 1col
anpp_forb_nitrogen_vwc_green_aic<-merge(anpp_forb_nitrogen_vwc_aic,green_aic,all.x=TRUE) #13cols, good

#add cheno cover
head(anpp_forb_nitrogen_vwc_green_aic) #13cols
head(sppchange_aic) #should add 1col
anpp_forb_nitrogen_vwc_green_cheno_aic<-merge(anpp_forb_nitrogen_vwc_green_aic,sppchange_aic,all.x=TRUE) #14cols, good

#add total cover
head(anpp_forb_nitrogen_vwc_green_cheno_aic) #13cols
head(total_cover_aic) #should add 1col
total_aic<-merge(anpp_forb_nitrogen_vwc_green_cheno_aic,total_cover_aic,all.x=TRUE) #15cols, good


#3. Subset dataframe for data of interest
#subset data to the 2022 deluge treatment only
  #we are trying to figure out what predictors drove the synergistic interaction
#years should always be assessed separately due to high variability
aic2022_all<-subset(total_aic,year==2022&treatment=="N+Deluge")


#4. Select the best predictor among the variables that are correlated
#a. generate aic values for nitrate, totalN
nitrogen_pred1<-lm(anpp~Nitrate, data=aic2022_all)
extractAIC(nitrogen_pred1)
nitrogen_pred2<-lm(anpp~Total_N, data=aic2022_all)
extractAIC(nitrogen_pred2)

AIC(nitrogen_pred1,nitrogen_pred2)

#Total_N is the best predictor, remove others
aic2022_all<-aic2022_all %>%
  dplyr::select(-Nitrate,-Ammonium,-nitrogen)

#b. generate aic values for forb proportion, cheno cover, total cover
#note: cheno cover has na values, convert to zeros
aic2022_all$cheno_cover[is.na(aic2022_all$cheno_cover)] <- 0
comm_pred1<-lm(anpp~cheno_cover, data=aic2022_all)
extractAIC(comm_pred1)
comm_pred2<-lm(anpp~forb_prct, data=aic2022_all)
extractAIC(comm_pred2)
comm_pred3<-lm(anpp~total_cover, data=aic2022_all)
extractAIC(comm_pred3)

AIC(comm_pred1,comm_pred2,comm_pred3)

#forb percent is the best predictor, remove others
aic2022_all<-aic2022_all %>%
  dplyr::select(-cheno_cover,-total_cover)


#5. Create the model and analyze with AIC
#First assess for correlation
library(corrplot)
corrdata<-aic2022_all[,6:10]
pairs(corrdata)
cor.test(corrdata$vwc,corrdata$gcc_mean)
corr_matrix<-cor(corrdata)
corrplot(corr_matrix)

#There is a high correlation of vwc and gcc_mean
comm_pred1<-lm(anpp~gcc_mean, data=aic2022_all)
extractAIC(comm_pred1)
comm_pred2<-lm(anpp~vwc, data=aic2022_all)
extractAIC(comm_pred2)

AIC(comm_pred1,comm_pred2)
  #Drop ggc_mean from further analysis

#Run the full model
model_aic<-lmer(anpp~Total_N+vwc+forb_prct
              + (1|block),
              data=aic2022_all,REML=FALSE)
Anova(model_aic,type=3)

#Conduct model selection
  #na's are not great for model selection--stop model if na's present
options(na.action="na.fail")
selected_model<-dredge(model_aic,rank="AIC")
selected_model
sw(selected_model) #sliding window is the updated "importance" function
  #Add R2 to table
mod<-get.models(selected_model, subset = TRUE)
R2_list<-lapply(mod, r.squaredGLMM)
selected_model$R2_marginal <- sapply(R2_list, function(x) x[1])
selected_model$R2_conditional <- sapply(R2_list, function(x) x[2])
selected_model

#Extract estimates, uncertainty, and p values from each model
  #Only the top model
top_model<-get.models(selected_model, 1)[[1]]
summary(top_model)
confint(top_model, method="Wald")

  #Model-average of those with deltaAIC below 2
#Get model-averaged estimates for all models with ΔAIC < 2 models
top_models<-subset(selected_model, delta < 2)
avg_est<-model.avg(top_models)
summary(avg_est)
confint(avg_est)


#Individual model assessment:
r2<-partR2(model_temp, data=aic2022_all, partvars=c("forb_prct"), R2_type= "marginal",nboot=20)
r2

model_temp<-lmer(anpp~vwc
                + (1|block),
                data=aic2022_all,REML=FALSE)
Anova(model_temp,type=2)
anova(model_temp)

######

####Mechanisms of production synergy: compile Figure 5
######

#FigA: Visualize proportion of forb prop to total anpp
#Generate ANPP dataframe
anpp_PP <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X,-woody,-dead,-anpp_with_woody,-cactus)

#Rename treatments for faceting ease
anpp_PP$treatment[anpp_PP$treatment == "N"]<-"Nitrogen"
anpp_PP$treatment[anpp_PP$treatment == "PP1"]<-"N+Deluge"
anpp_PP$treatment[anpp_PP$treatment == "PP2"]<-"N+Deluge"

#Create a "forb percentage" column by diving forbs by total anpp
anpp_PP$forb_prct<-(anpp_PP$forb/anpp_PP$anpp)

#fix str
anpp_PP<-anpp_PP %>% 
  dplyr::mutate_at(c("block","plot","nitrogen","treatment"),as.factor)
str(anpp_PP)

#Run correlation tests for manuscript
test_2021_N<-filter(anpp_PP,year==2021&treatment=="Nitrogen")
test_2021_PP<-filter(anpp_PP,year==2021&treatment=="N+Deluge")
test_2022_N<-filter(anpp_PP,year==2022&treatment=="Nitrogen")
test_2022_PP<-filter(anpp_PP,year==2022&treatment=="N+Deluge")
test_2021<-filter(anpp_PP,year==2021)
test_2022<-filter(anpp_PP,year==2022)

cor.test(test_2022_PP$forb_prct,test_2022_PP$anpp)

#Plotting
anpp_PP$treatment <- factor(anpp_PP$treatment, levels=c("Nitrogen", "N+Deluge"))
mechs_a<-ggplot(data=subset(anpp_PP),  aes(x=forb_prct, y=anpp, color=nitrogen, group=year, linetype = treatment, shape = treatment)) +
  geom_point(size=2) +
  geom_smooth(data = subset(anpp_PP, !(treatment == "N+Deluge"&year=="2021")), aes(group = treatment), method = "lm", size = 1, color = "darkorange3") + 
  #stat_cor(label.y = 215) +
  facet_grid(year~., scales = "free") +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=12),
        axis.text.y=element_text(color="black",size=15),
        axis.title=element_text(color="black",size=15),
        plot.title = element_text(color = "black", size = 20, hjust = 0),
        panel.grid.minor = element_line(colour = NA), 
        strip.text.y = element_text(size = 20),
        legend.position = c(.4,.35),
        legend.title = element_text(size=12),
        #legend.direction= "horizontal",
        legend.text = element_text(size=12),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1)) +
  geom_hline(yintercept = 0, size = 1, linetype="dashed") +
  geom_hline(yintercept = -Inf, size = 1.5) +
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("gray60","#6f9954","#89ab7c","#87bcbd","#969bc7","#6b6ca3","#434475","#2c2d54")) +
  labs(y=expression(ANPP~(g~m^-2)))+ #y label
  labs(x=expression(Proportion~Forbs)) + #x label
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 1)) +  # Set y-axis limits and breaks
  guides(color =  guide_legend(position = "right",title = expression(N~(g~m^-2))),
         label.theme = element_text(size = 12, lineheight = 1.2), 
         title.position = "top",
         title.theme = element_text(size = 15, lineheight = 1.2),
         linetype = guide_legend(title = NULL),
         shape = guide_legend(title = NULL)) +
  ggtitle("a.")


#FigB: proportion plot
#Generate ANPP dataframe
anpp_all <- read.csv("MAIN_ChANGE_PulsePress_SGS_ANPP_2021to2022.csv")%>%
  dplyr::select(-X,-woody,-dead,-anpp_with_woody,-cactus)

#Create a forb:grass value
anpp_all$forb_prop<-(anpp_all$forb/anpp_all$anpp)

#Correct facet wrap labels for consistency
anpp_all$treatment[anpp_all$treatment == "PP1"]<-"N+D"
anpp_all$treatment[anpp_all$treatment == "PP2"]<-"N+D"

#Gather to long format
anpp_long<-anpp_all%>%
  gather(functional,anpp_func,"forb":"grass")

#Create a "Percent Forbs" col
anpp_all$Forb<-(anpp_all$forb/anpp_all$anpp)
anpp_all$Grass<-(anpp_all$grass/anpp_all$anpp)

#Gather to long format
anpp_long<-anpp_all%>%
  gather(functional,proportion,"Forb":"Grass")

#Fix str
anpp_long<-anpp_long %>% 
  dplyr::mutate_at(c("site","year","plot","nitrogen","block"),as.factor)
str(anpp_long)
anpp_long$treatment <- factor(anpp_long$treatment, levels=c("N", "N+D"))

#Calculate magnitudes
ddply(subset(anpp_all), c("year"), summarize, 
      n= length (forb_prop),
      mean = mean(forb_prop),
      sd = sd(forb_prop),
      se = sd/(sqrt(n)))

#Need to fix error bars for stacked bar plots: grass is stacked on top, so we need to add the forb average to the grass
#The following code with generated by the help of the stat lab
data_func<-barGraphStats(data=subset(anpp_long, site=='SGS'), variable="proportion", byFactorNames=c("treatment","year","functional"))
data_func2 <- data_func %>%
  dplyr::select(-sd, -se) %>%
  tidyr::pivot_wider(names_from = functional, 
                     values_from = mean) %>%
  dplyr::mutate(Grass = Grass + Forb) %>%
  tidyr::pivot_longer(c(Forb, Grass), 
                      names_to = "functional", 
                      values_to = "mean_alt")

data_func_final <- dplyr::left_join(data_func, data_func2) %>%
  dplyr::mutate(functional = factor(functional, 
                                    levels = c("Grass", "Forb")))

#Only need the error bars from the Forb group, so set the mean's and se's to zero
#in this type of plot, the errorbars of grass/forb mirror eachother, so we need to remove a set of errorbars
data_func_final$se_alt<-ifelse(data_func_final$functional=="Grass",data_func_final$se=="",data_func_final$se)
data_func_final$mean_alt<-ifelse(data_func_final$functional=="Grass",data_func_final$mean_alt=="",data_func_final$mean_alt)

mechs_b<-ggplot(data=subset(data_func_final), aes(x=treatment, y=mean, fill=functional)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean_alt-se_alt, ymax=mean_alt+se_alt), width=0.2) +
  facet_grid(year~., scales = "free") +
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(size=20), #set as 25
        panel.grid.minor = element_line(colour = NA), 
        axis.title.y = element_text(size=15,color = "black"),
        axis.title.x = element_blank(), 
        legend.position= c(.5,.2),
        legend.text = element_text(size = 15),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = 1),
        legend.title = element_blank(),
        strip.text.x = element_text(color="black",size = 20),
        strip.text.y = element_text(color="black",size = 20))+
  geom_vline(aes(xintercept=-Inf), size=1.5) +
  geom_hline(aes(yintercept=-Inf), linewidth=2) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 1)) +  # Set y-axis limits and breaks
  scale_fill_manual(values=c("#6f9954","#434475")) +
  labs(y=expression(Functional~Group~Proportion)) + #y label
  labs(x=expression(N(g~m^-2))) +  #x label
  ggtitle("b.")


#FigC: abundance plots
#See code above to see how these were determined
top_3_sgs_2021_del<-c("Elymus elymoides", "Vulpia octoflora", "Bromus tectorum")
top_3_sgs_2022_del<-c("Elymus elymoides", "Chenopodium sp.", "Salsola tragus")

#Filter to year of interest
sppcompdata_PP1 <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(year==2021,site=="sgs", species%in%top_3_sgs_2021_del)

sppcompdata_PP2 <- read.csv("MAIN_ChANGE_PulsePress_SGS_sppcomp_2021to2022.csv")%>%
  filter(year==2022&treatment!="pulse press 1",site=="sgs", species%in%top_3_sgs_2022_del)

sppcomp_sgs<-full_join(sppcompdata_PP1,sppcompdata_PP2)

#Change the nitrogen treatment name, so the wide format isn't confused
sppcomp_sgs$treatment[sppcomp_sgs$treatment == "nitrogen"]<-"Nitrogen"
sppcomp_sgs$treatment[sppcomp_sgs$treatment == "pulse press 1"]<-"N+Deluge"
sppcomp_sgs$treatment[sppcomp_sgs$treatment == "pulse press 2"]<-"N+Deluge"

#prep for plotting
sppcomp_sgs_temp<-sppcomp_sgs
sppcomp_sgs_temp$treatment[sppcomp_sgs_temp$treatment == "Nitrogen"]<-"N"
sppcomp_sgs_temp$treatment[sppcomp_sgs_temp$treatment == "N+Deluge"]<-"N+D"

#Order factored treatments
sppcomp_sgs$nitrogen<-as.factor(sppcomp_sgs$nitrogen)
sppcomp_sgs$treatment <- factor(sppcomp_sgs$treatment, levels=c("Nitrogen", "N+Deluge"))

#abundance plots
mechs_c_pt1<-ggplot(data=barGraphStats(data=subset(sppcomp_sgs), variable="cover", byFactorNames=c("species","year","nitrogen")), aes(x=nitrogen, y=mean, color=species, group=species)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5) +
  geom_point(size=3) +
  geom_line(size=1) +
  facet_grid(year~., scales="free")+
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=12, angle = 90, vjust = 0.5),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(color = "black", size = 20, hjust = 0),
        axis.title = element_text(color="black",size=15),
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=15),
        strip.text = element_blank()) +
  geom_hline(aes(yintercept=-Inf), size=2) +
  geom_hline(aes(yintercept=0), size=1, linetype="dashed") + 
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("#2c2d54","darkorange3","#6f9954","#969bc7","gray50"))+ #alternate palette
  labs(y=expression(Abundance)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) +
  ggtitle("c.")

#plot deluge effect
sppcomp_sgs_temp$species[sppcomp_sgs_temp$species == "Bromus tectorum"]<-"B.tectorum"
sppcomp_sgs_temp$species[sppcomp_sgs_temp$species == "Chenopodium sp."]<-"Chenopodium sp."
sppcomp_sgs_temp$species[sppcomp_sgs_temp$species == "Elymus elymoides"]<-"E.elymoides"
sppcomp_sgs_temp$species[sppcomp_sgs_temp$species == "Salsola tragus"]<-"S.tragus"
sppcomp_sgs_temp$species[sppcomp_sgs_temp$species == "Vulpia octoflora"]<-"V.octoflora"

mechs_c_pt2<-ggplot(data=barGraphStats(data=subset(sppcomp_sgs_temp), variable="cover", byFactorNames=c("species","year","treatment")), aes(x=treatment, y=mean, color=species, group=species)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3) +
  geom_point(size=3) +
  geom_line(size=1) +
  facet_grid(year~.)+
  theme(panel.background = element_rect(fill = "white",colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        axis.line.x = element_line(size = 1.5, colour="black"), 
        axis.line.y = element_line(size = 1.5, colour="black"),
        axis.ticks=element_line(color="black",size=1.5), 
        axis.text.x=element_text(color="black", size=15),
        axis.text.y=element_text(color="black",size=15),
        plot.title = element_text(color = "black", size = 20, hjust = 0),
        axis.title.x = element_text(color="white",size=15),
        axis.title.y = element_blank(),
        legend.position="right",
        legend.title = element_blank(),
        legend.text = element_text(color="black",size=12, face = "italic"),
        strip.text = element_text(color="black",size = 20)) +
  geom_hline(aes(yintercept=0), size=1, linetype="dashed") + 
  geom_hline(aes(yintercept=-Inf), size=2) +
  scale_color_manual(name=expression(Nitrogen~(g~m^-2)),values=c("#2c2d54","darkorange3","#6f9954","#969bc7","gray50"))+ #alternate palette
  labs(y=expression(Abundance)) + #y label
  labs(x=expression(Nitrogen~(g~m^-2))) +
  ggtitle("")
#guides(color = guide_legend(label.theme = element_text(size = 15, lineheight = 1.2), label.position = "top",  title.position = "top"))

#arrange
mechs_ab<-ggarrange(mechs_a, mechs_b, ncol=2, nrow=1, widths = c(1.7,1))
mechs_c<-ggarrange(mechs_c_pt1,mechs_c_pt2, ncol=2, nrow=1, widths = c(1,1.4))
ggarrange(mechs_ab,mechs_c,ncol=1, nrow=2)

######


