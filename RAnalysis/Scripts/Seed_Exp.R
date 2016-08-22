#Effects of Ocean Acidification on Phenotype and DNA Methylation in Juvenile Geoduck
#Data published in
#Title:
#Contact: Hollie Putnam hollieputnam@gmail.com ; Steven Roberts sr320@u.washington.edu
#Supported by: NOAA OA
#last modified 20160708 H Putnam
#See Readme file for project details 
#See metadata file for details on data files and equipment 

rm(list=ls()) # removes all prior objects

#Read in required libraries
library(seacarb) #seawater carbonate chemistry
library(reshape) #reshape data
library(plotrix) #functions in tapply
library(ggplot2)
library(plyr)
library(gridExtra)

#Required Data files
#Flow_Juvenile_Geoduck.csv
#Avtech_data_Juvenile_Geoduck_Exp1.csv
#Avtech_Data_Juvenile_Geoduck_ICG.csv
#Hobo_Temperature_Juvenile_Geoduck_OCG.csv
#Avtech_data_Juvenile_Geoduck_Exp2.csv
#pH_Calibration_Files/
#SW_Chem_Juvenile_Geoduck.csv
#Cell_Counts_Juvenile_Geoduck_Exp1.csv
#Cell_Counts_Juvenile_Geoduck_Exp2.csv
#Size_Juvenile_Geoduck.csv

#############################################################
setwd("/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/Data/") #set working directory
mainDir<-'/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/' #set main directory

##### DESCRIPTIVE STATISTICS FOR EXPERIMENTAL DESIGN #####
Flow1 <- read.csv("Flow1_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
#flow values measured in ml/10 sec
Flow1$Rate <- (((Flow1$Rate1 + Flow1$Rate2)/2)*6)/1000*60
mean(Flow1$Rate)
sd(Flow1$Rate)/sqrt(length(Flow1$Rate))
Flow1.Rate <- cat("Flow (mean±sem, n) =", mean(Flow1$Rate),  sd(Flow1$Rate)/sqrt(length(Flow1$Rate)), length(Flow1$Rate))

FL1 <- aov(Rate ~ Treatment, data=Flow1)
anova(FL1)
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
hist(FL1$residuals)
boxplot(FL1$residuals)
plot(FL1)

Flow2 <- read.csv("Flow2_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
#flow values measured in ml/10 sec
Flow2$Rate <- ((Flow2$Rate1)*4)/1000*60
mean(Flow2$Rate)
sd(Flow2$Rate)/sqrt(length(Flow2$Rate))
Flow2.Rate <- cat("Flow (mean±sem, n) =", mean(Flow2$Rate),  sd(Flow2$Rate)/sqrt(length(Flow2$Rate)), length(Flow2$Rate))

FL2 <- aov(Rate ~ Treatment, data=Flow2)
anova(FL2)
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
hist(FL2$residuals)
boxplot(FL2$residuals)
plot(FL2)

##### CONTINUOUS WATERBATH TEMPERATURE AND PH HEADER TANK EXPOSURE 1 #####
#Avtech_data_Juvenile_Geoduck_Exp1.csv
Exp1.pH <- read.csv("Avtech_data_Juvenile_Geoduck_Exp1.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
Exp1.pH$Date.Time <-as.POSIXct(Exp1.pH$Date.Time, format="%m/%d/%y %H:%M")
Exp1.pH$Date <- as.Date(Exp1.pH$Date.Time) #already got this one from the answers above
Exp1.pH$Time <- format(as.POSIXct(Exp1.pH$Date.Time) ,format = "%H:%M:%S") 

pH.low.m <- aggregate(Header2.pH ~ Date, data = Exp1.pH, mean)
pH.low.se <- aggregate(Header2.pH ~ Date, data = Exp1.pH, std.error)
pH.med.m <- aggregate(Header1.pH ~ Date, data = Exp1.pH, mean)
pH.med.se <- aggregate(Header1.pH ~ Date, data = Exp1.pH, std.error)
pH.amb.m <- aggregate(Reservoir.pH ~ Date, data = Exp1.pH, mean)
pH.amb.se <- aggregate(Reservoir.pH ~ Date, data = Exp1.pH, std.error)
pH.low <- cbind(pH.low.m, pH.low.se$Header2.pH)
pH.low$Treatment <- "Super.Low"
colnames(pH.low) <- c("Date", "mean", "se", "Treatment")
pH.med <- cbind(pH.med.m, pH.med.se$Header1.pH)
pH.med$Treatment <- "Low"
colnames(pH.med) <- c("Date", "mean", "se", "Treatment")
pH.amb <- cbind(pH.amb.m, pH.amb.se$Reservoir.pH)
pH.amb$Treatment <- "Ambient"
colnames(pH.amb) <- c("Date", "mean", "se", "Treatment")
daily.pH <- rbind(pH.low, pH.med, pH.amb)
daily.pH

Fig.E1.daily.pH <- ggplot(daily.pH, aes(x=Date, y=mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=daily.pH$mean-daily.pH$se, ymax=daily.pH$mean+daily.pH$se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 2, position = position_dodge(width = 0.05)) +
  xlab("Time") +
  ylab("pH Total Scale") +
  ylim(6.8,8.1) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position=c(.3, .65), #set legend location
        legend.text = element_text(size = 5),
        legend.key = element_blank(),
        legend.title = element_text(size=6, face="bold")) +
  ggtitle("A) Exposure 1") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.E1.daily.pH

#Exposure 1 Waterbath Temperatures
Exp1.Temp <- Exp1.pH
Temp.low.m <- aggregate(Header2.Temp ~ Date, data = Exp1.Temp, mean)
Temp.low.se <- aggregate(Header2.Temp ~ Date, data = Exp1.Temp, std.error)
Temp.med.m <- aggregate(Header1.Temp ~ Date, data = Exp1.Temp, mean)
Temp.med.se <- aggregate(Header1.Temp ~ Date, data = Exp1.Temp, std.error)
Temp.amb.m <- aggregate(Reservoir.Temp ~ Date, data = Exp1.Temp, mean)
Temp.amb.se <- aggregate(Reservoir.Temp ~ Date, data = Exp1.Temp, std.error)
Temp.low <- cbind(Temp.low.m, Temp.low.se$Header2.Temp)
Temp.low$Treatment <- "Super.Low"
colnames(Temp.low) <- c("Date", "mean", "se", "Treatment")
Temp.med <- cbind(Temp.med.m, Temp.med.se$Header1.Temp)
Temp.med$Treatment <- "Low"
colnames(Temp.med) <- c("Date", "mean", "se", "Treatment")
Temp.amb <- cbind(Temp.amb.m, Temp.amb.se$Reservoir.Temp)
Temp.amb$Treatment <- "Ambient"
colnames(Temp.amb) <- c("Date", "mean", "se", "Treatment")
daily.Temp <- rbind(Temp.low, Temp.med, Temp.amb)
daily.Temp

Fig.E1.daily.Temp <- ggplot(daily.Temp, aes(x=Date, y=mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=daily.Temp$mean-daily.Temp$se, ymax=daily.Temp$mean+daily.Temp$se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 2, position = position_dodge(width = 0.05)) +
  xlab("Time") +
  ylab("Temperature °C") +
  ylim(10, 20) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),#Set the plot background
        legend.position=c(.4, .7), #set legend location
        legend.text = element_text(size = 5),
        legend.key = element_blank(),
        legend.title = element_text(size=6, face="bold")) +
  ggtitle("A) Exposure 1") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.E1.daily.Temp

##### CONTINUOUS TEMPERATURE AND PH OF INDOOR COMMON GARDEN #####
#Avtech_Data_Juvenile_Geoduck_ICG.csv
ICG <- read.csv("Avtech_Data_Juvenile_Geoduck_ICG.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
ICG$Date.Time <-as.POSIXct(ICG$Date.Time, format="%m/%d/%y %H:%M")
ICG$Date <- as.Date(ICG$Date.Time) 
ICG$Time <- format(as.POSIXct(ICG$Date.Time) ,format = "%H:%M:%S") 

ICG.Temps.m <- aggregate(Temp.C ~ Date, data=ICG, mean)
ICG.Temps.se <- aggregate(Temp.C ~ Date, data=ICG, std.error)
ICG.Temps <- cbind(ICG.Temps.m, ICG.Temps.se$Temp.C)
colnames(ICG.Temps) <- c("Date", "mean", "se")

Fig.ICG.Temp <- ggplot(ICG.Temps, aes(x=Date, y=mean)) + 
  geom_errorbar(aes(ymin=ICG.Temps$mean-ICG.Temps$se, ymax=ICG.Temps$mean+ICG.Temps$se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(x=Date, y=mean), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(x=Date, y=mean), size = 2, position = position_dodge(width = 0.05)) +
  xlab("Time") +
  ylab("Temperature °C") +
  ylim(10,20) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.key = element_blank()) + #remove legend background
  ggtitle("B) Indoor Common Garden") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.ICG.Temp

ICG.pH.m <- aggregate(pH.Total ~ Date, data=ICG, mean)
ICG.pH.se <- aggregate(pH.Total ~ Date, data=ICG, std.error)
ICG.pH <- cbind(ICG.pH.m, ICG.pH.se$pH.Total)
colnames(ICG.pH) <- c("Date", "mean", "se")

Fig.ICG.pH <- ggplot(ICG.pH, aes(x=Date, y=mean)) + 
  geom_errorbar(aes(ymin=ICG.pH$mean-ICG.pH$se, ymax=ICG.pH$mean+ICG.pH$se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(x=Date, y=mean), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(x=Date, y=mean), size = 2, position = position_dodge(width = 0.05)) +
  xlab("Time") +
  ylab("pH Total Scale") +
  ylim(6.8, 8.1) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.key = element_blank()) + #remove legend background
  ggtitle("B) Indoor Common Garden") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.ICG.pH

##### CONTINUOUS TANK TEMPERATURE OUTDOOR COMMON GARDEN #####
#Hobo_Temperature_Juvenile_Geoduck_OCG.csv
OCG.Temp <- read.csv("Hobo_Temperature_Juvenile_Geoduck_OCG.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
OCG.Temp$Date.Time <-as.POSIXct(OCG.Temp$Date.Time, format="%m/%d/%y %H:%M")
OCG.Temp$Date <- as.Date(OCG.Temp$Date.Time) 
OCG.Temp$Time <- format(as.POSIXct(OCG.Temp$Date.Time) ,format = "%H:%M:%S") 

OCG.Temps.m <- aggregate(Temp.C ~ Date*Treatment, data=OCG.Temp, mean)
OCG.Temps.se <- aggregate(Temp.C ~  Date*Treatment, data=OCG.Temp, std.error)
OCG.Temps <- cbind(OCG.Temps.m, OCG.Temps.se$Temp.C)
colnames(OCG.Temps) <- c("Date", "Treatment", "mean", "se")

Fig.OCG.Temp <- ggplot(OCG.Temps, aes(x=Date, y=mean), group=Treatment) + 
  geom_errorbar(aes(ymin=OCG.Temps$mean-OCG.Temps$se, ymax=OCG.Temps$mean+OCG.Temps$se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 2, position = position_dodge(width = 0.05)) +
  xlab("Time") +
  ylab("Temperature °C") +
  ylim(10, 20) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position = "none") + #remove legend background
  ggtitle("B) Indoor Common Garden") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.OCG.Temp

#Blank plot for plotting grid
Fig.OCG.pH <- ggplot(OCG.Temps, aes(x=Date, y=mean), group=Treatment) + #blank plot for plotting grid
  geom_blank() +
  xlab("Time") +
  ylab("pH Total Scale") +
  ylim(6.8, 8.0) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position = "none") + #Removes the legend
  ggtitle("C) Outdoor Common Garden") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.OCG.pH

##### CONTINUOUS WATERBATH TEMPERATURE AND PH HEADER TANK EXPOSURE 2 #####
#Avtech_data_Juvenile_Geoduck_Exp2.csv
#Exposure 2 Waterbath pH
Exp2.data <- read.csv("Avtech_data_Juvenile_Geoduck_Exp2.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
Exp2.data$Date.Time <-as.POSIXct(Exp2.data$Date.Time, format="%m/%d/%y %H:%M")
Exp2.data$Date <- as.Date(Exp2.data $Date.Time) #already got this one from the answers above
Exp2.data$Time <- format(as.POSIXct(Exp2.data $Date.Time) ,format = "%H:%M:%S") 

E2.pH.low.m <- aggregate(pH.Low ~ Date, data = Exp2.data, mean)
E2.pH.low.se <- aggregate(pH.Low ~ Date, data = Exp2.data, std.error)
E2.pH.amb.m <- aggregate(pH.Ambient ~ Date, data = Exp2.data, mean)
E2.pH.amb.se <- aggregate(pH.Ambient ~ Date, data = Exp2.data, std.error)
E2.pH.low <- cbind(E2.pH.low.m, E2.pH.low.se$pH.Low)
E2.pH.low$Treatment <- "Low"
colnames(E2.pH.low) <- c("Date", "mean", "se", "Treatment")
E2.pH.amb <- cbind(E2.pH.amb.m, E2.pH.amb.se$pH.Ambient)
E2.pH.amb$Treatment <- "Ambient"
colnames(E2.pH.amb) <- c("Date", "mean", "se", "Treatment")
E2.daily.pH <- rbind(E2.pH.low, E2.pH.amb)
E2.daily.pH

Fig.E2.daily.pH <- ggplot(E2.daily.pH, aes(x=Date, y=mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=E2.daily.pH$mean-E2.daily.pH$se, ymax=E2.daily.pH$mean+E2.daily.pH$se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 2, position = position_dodge(width = 0.05)) +
  xlab("Time") +
  ylab("pH Total Scale") +
  ylim(6.8,8.1) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position="none") +
  ggtitle("D) Exposure 2") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.E2.daily.pH

#Exposure 2 Waterbath Temperatures
E2.Temp.low.m <- aggregate(Temp.Low ~ Date, data = Exp2.data, mean)
E2.Temp.low.se <- aggregate(Temp.Low ~ Date, data = Exp2.data, std.error)
E2.Temp.amb.m <- aggregate(Temp.Ambient ~ Date, data = Exp2.data, mean)
E2.Temp.amb.se <- aggregate(Temp.Ambient ~ Date, data = Exp2.data, std.error)
E2.Temp.low <- cbind(E2.Temp.low.m, E2.Temp.low.se$Temp.Low)
E2.Temp.low$Treatment <- "Low"
colnames(E2.Temp.low) <- c("Date", "mean", "se", "Treatment")
E2.Temp.amb <- cbind(E2.Temp.amb.m, E2.Temp.amb.se$Temp.Ambient)
E2.Temp.amb$Treatment <- "Ambient"
colnames(E2.Temp.amb) <- c("Date", "mean", "se", "Treatment")
E2.daily.Temp <- rbind(E2.Temp.low, E2.Temp.amb)
E2.daily.Temp

Fig.E2.daily.Temp <- ggplot(E2.daily.Temp, aes(x=Date, y=mean, group=Treatment)) + 
  geom_errorbar(aes(ymin=E2.daily.Temp$mean-E2.daily.Temp$se, ymax=E2.daily.Temp$mean+E2.daily.Temp$se), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 2, position = position_dodge(width = 0.05)) +
  xlab("Time") +
  ylab("Temperature °C") +
  ylim(10,20) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position="none") +
  ggtitle("D) Exposure 2") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.E2.daily.Temp

##### SEAWATER CHEMISTRY ANALYSIS FOR DISCRETE MEASUREMENTS #####

##### pH Tris Calibration Curves
#Data to calculate conversion equations from mV to total scale using tris standard for pH probe
path <-("/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/Data/pH_Calibration_Files/")

#list all the file names in the folder to get only get the csv files
file.names<-list.files(path = path, pattern = "csv$")

pH.cals <- data.frame(matrix(NA, nrow=length(file.names), ncol=4, dimnames=list(file.names,c("Date", "Intercept", "Slope","R2")))) #generate a 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Calib.Data <-read.table(file.path(path,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
  model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
  coe <- coef(model) #extracts the coeffecients
  R <- summary(model)$r.squared #extracts the R2
  pH.cals[i,2:3] <- coe #inserts coef in the dataframe
  pH.cals[i,4] <- R #inserts R2 in the dataframe
  pH.cals[i,1] <- substr(file.names[i],1,8) #stores the file name in the Date column
}

colnames(pH.cals) <- c("Calib.Date",  "Intercept",  "Slope", "R2")
pH.cals

# read in total alkalinity, temperature, and salinity
SW.chem <- read.csv("SW_Chem_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

#merge with Seawater chemistry file
SW.chem <- merge(pH.cals, SW.chem, by="Calib.Date")

#constants for use in pH calculation 
R <- 8.31447215 #gas constant in J mol-1 K-1 
F <-96485.339924 #Faraday constant in coulombs mol-1

mvTris <- SW.chem$Temperature*SW.chem$Slope+SW.chem$Intercept #calculate the mV of the tris standard using the temperature mv relationships in the measured standard curves 
STris<-27.5 #salinity of the Tris from local batch made for Manchester WA salinity conditions in Spring 2016. Batch date = 
phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(SW.chem$Temperature+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(SW.chem$Temperature+273.15)-0.11149858*(SW.chem$Temperature+273.15) #calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem$pH.Total<-phTris+(mvTris/1000-SW.chem$pH.MV/1000)/(R*(SW.chem$Temperature+273.15)*log(10)/F) #calculate the pH on the total scale (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem <- na.omit(SW.chem)

##### SEACARB CALCULATIONS #####

#Calculate CO2 parameters using seacarb
carb.output <- carb(flag=8, var1=SW.chem$pH.Total, var2=SW.chem$TA/1000000, S= SW.chem$Salinity, T=SW.chem$Temperature, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb

carb.output$ALK <- carb.output$ALK*1000000 #convert to µmol kg-1
carb.output$CO2 <- carb.output$CO2*1000000 #convert to µmol kg-1
carb.output$HCO3 <- carb.output$HCO3*1000000 #convert to µmol kg-1
carb.output$CO3 <- carb.output$CO3*1000000 #convert to µmol kg-1
carb.output$DIC <- carb.output$DIC*1000000 #convert to µmol kg-1

carb.output <- cbind(SW.chem$Measure.Date,  SW.chem$Exposure, SW.chem$Tank,  SW.chem$Treatment, carb.output) #combine the sample information with the seacarb output
colnames(carb.output) <- c("Date", "Exposure", "Tank",  "Treatment",	"flag",	"Salinity",	"Temperature",	"Pressure",	"pH",	"CO2",	"pCO2",	"fCO2",	"HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat", 	"Calcite.Sat") #Rename columns to describe contents
carb.output <- subset(carb.output, select= c("Exposure","Date",  "Tank",  "Treatment",	"Salinity",	"Temperature",		"pH",	"CO2",	"pCO2",	"HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat"))

##### SEAWATER CHEMISTRY DESCRIPTIVE STATISTICS EXPOSURE 1 AND 2 #####

ccarb <- melt(carb.output[,c(1,3,4,5:14)], id.vars=c("Exposure", "Treatment", "Tank")) #reshape data into long format

#Calculate descriptive stats by Tank
Exp1 <-subset(ccarb, Exposure == "Exposure1") #separate out exposure 1 for all data
Exp2 <-subset(ccarb, Exposure == "Exposure2") #separate out exposure 2 for all data

SWC.Tanks <- ddply(ccarb, c("Exposure", "Tank", "variable"), summarise,
      N = length(na.omit(value)), #count the sample size removing NA
      mean = mean(value), #calculate average 
      sem = sd(value)/sqrt(N)) #calcualte the standard error of the mean

#Test for tank and treatment differneces in Temperature and Total Alkalinity in Exposure 1
Exp1.Temp <-subset(Exp1, variable=="Temperature") #separate out exposure 1 for all data
temp1.tank <- aov(value ~Tank, data=Exp1.Temp)
anova(temp1.tank)
par(mfrow=c(3,3))
hist(temp1.tank$residuals)
boxplot(temp1.tank$residuals)
plot(temp1.tank)

temp1.trt <- aov(value ~Treatment, data=Exp1.Temp)
anova(temp1.trt)
par(mfrow=c(3,3))
hist(temp1.trt$residuals)
boxplot(temp1.trt$residuals)
plot(temp1.trt)

Exp1.TA <-subset(Exp1, variable=="TA") #separate out exposure 1 for all data
TA1.tank <- aov(value ~Tank, data=Exp1.TA)
anova(TA1.tank)
par(mfrow=c(3,3))
hist(TA1.tank$residuals)
boxplot(TA1.tank$residuals)
plot(TA1.tank)

TA1.trt <- aov(value ~Treatment, data=Exp1.TA)
anova(TA1.trt)
par(mfrow=c(3,3))
hist(TA1.trt$residuals)
boxplot(TA1.trt$residuals)
plot(TA1.trt)

#Test for tank and treatment differneces in Temperature and Total Alkalinity in Exposure 2
Exp2.Temp <-subset(Exp2, variable=="Temperature") #separate out exposure 2 for all data
temp2.tank <- aov(value ~Tank, data=Exp2.Temp)
anova(temp2.tank)
par(mfrow=c(3,3))
hist(temp2.tank$residuals)
boxplot(temp2.tank$residuals)
plot(temp2.tank)

temp2.trt <- aov(value ~Treatment, data=Exp2.Temp)
anova(temp2.trt)
par(mfrow=c(3,3))
hist(temp2.trt$residuals)
boxplot(temp2.trt$residuals)
plot(temp2.trt)

Exp2.TA <-subset(Exp2, variable=="TA") #separate out exposure 2 for all data
TA2.tank <- aov(value ~Tank, data=Exp2.TA)
anova(TA2.tank)
par(mfrow=c(3,3))
hist(TA2.tank$residuals)
boxplot(TA2.tank$residuals)
plot(TA2.tank)

TA2.trt <- aov(value ~Treatment, data=Exp2.TA)
anova(TA2.trt)
par(mfrow=c(3,3))
hist(TA2.trt$residuals)
boxplot(TA2.trt$residuals)
plot(TA2.trt)

#Calculate descriptive stats by Treatment
SWC.Treatments <- ddply(ccarb, c("Exposure", "Treatment", "variable"), summarise,
      N = length(na.omit(value)), #count the sample size removing NA
      mean = mean(value), #calculate average 
      sem = sd(value)/sqrt(N)) #calcualte the standard error of the mean

Exposure1 <-subset(SWC.Treatments, Exposure == "Exposure1") #separate out exposure 1
Exposure2 <-subset(SWC.Treatments, Exposure == "Exposure2") #separate out exposure 2

Exposure1.long <- reshape(Exposure1, idvar="Treatment", direction="wide", timevar = "variable", drop = c("Exposure", "N")) #reshape data format for table layout
Exposure2.long <- reshape(Exposure2, idvar="Treatment", direction="wide", timevar = "variable", drop = c("Exposure", "N")) #reshape data format for table layout

write.table (Exposure1.long, file="/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_oa/RAnalysis/Output/Seawater_chemistry_table_Output_Seed_exposure1.csv", sep=",", row.names = FALSE) #save data to output file
write.table (Exposure2.long, file="/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_oa/RAnalysis/Output/Seawater_chemistry_table_Output_Seed_exposure2.csv", sep=",", row.names = FALSE) #save data to output file

##### ALGAL FEED DENSITY COUNTS EXPOSURE 1 #####
#Cell_Counts_Juvenile_Geoduck_Exp1.csv
cells.1 <- read.csv("Cell_Counts_Juvenile_Geoduck_Exp1.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
cells.1$Avg.Cells <- rowMeans(cells.1[,c("Count1",  "Count2")], na.rm = TRUE) #calculate average of counts
cells.1$cell.num <- cells.1$Avg.Cells/cells.1$Volume.Counted #calculate density

avg.cells.tank.1 <- aggregate(cell.num ~ Tank, data=cells.1, mean)
se.cells.tank.1 <- aggregate(cell.num ~ Tank, data=cells.1, std.error)
avg.cells.trt.1 <- aggregate(cell.num ~ Treatment, data=cells.1, mean)
se.cells.trt.1 <- aggregate(cell.num ~ Treatment, data=cells.1, std.error)
Cell.Counts.1 <- cbind(avg.cells.tank.1, se.cells.tank.1$cell.num, avg.cells.trt.1$cell.num, se.cells.trt.1$cell.num)
colnames(Cell.Counts.1) <- c("Tank", "tank.avg", "tank.se", "trt.avg", "trt.se")

Fig.Exp1.cells <- ggplot(Cell.Counts.1, aes(x=Tank, y=tank.avg)) + 
  geom_errorbar(aes(ymin=Cell.Counts.1$tank.avg-Cell.Counts.1$tank.se, ymax=Cell.Counts.1$tank.avg+Cell.Counts.1$tank.se), colour="black", width=.1, position = position_dodge(width = 0.6)) +
  #geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.6)) +   
  geom_point(aes(), size = 3, position = position_dodge(width = 0.6)) +
  xlab("Tanks") +
  ylab("cells per ml") +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank()) 
Fig.Exp1.cells

cells1.tank <- aov(cell.num ~Tank, data=cells.1)
anova(cells1.tank)
par(mfrow=c(3,2))
hist(cells1.tank$residuals)
boxplot(cells1.tank$residuals)
plot(cells1.tank)

cells1.trt <- aov(cell.num ~Treatment, data=cells.1)
anova(cells1.trt)
par(mfrow=c(3,2))
hist(cells1.trt$residuals)
boxplot(cells1.trt$residuals)
plot(cells1.trt)

##### ALGAL FEED DENSITY COUNTS EXPOSURE 2 #####
#Cell_Counts_Juvenile_Geoduck_Exp2.csv
cells.2 <- read.csv("Cell_Counts_Juvenile_Geoduck_Exp2.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
cells.2$cell.num <- rowMeans(cells.2[,c("Count1",  "Count2", "Count3")], na.rm = TRUE) #calculate average of counts
avg.cells.tank <- aggregate(cell.num ~ Tank, data=cells.2, mean)
se.cells.tank <- aggregate(cell.num ~ Tank, data=cells.2, std.error)
avg.cells.trt <- aggregate(cell.num ~ Treatment, data=cells.2, mean)
se.cells.trt <- aggregate(cell.num ~ Treatment, data=cells.2, std.error)
Cell.Counts <- cbind(avg.cells.tank, se.cells.tank$cell.num, avg.cells.trt$cell.num, se.cells.trt$cell.num)
colnames(Cell.Counts) <- c("Tank", "tank.avg", "tank.se", "trt.avg", "trt.se")

Fig.Exp2.cells <- ggplot(Cell.Counts, aes(x=Tank, y=tank.avg)) + 
  geom_errorbar(aes(ymin=Cell.Counts$tank.avg-Cell.Counts$tank.se, ymax=Cell.Counts$tank.avg+Cell.Counts$tank.se), colour="black", width=.1, position = position_dodge(width = 0.6)) +
  #geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.6)) +   
  geom_point(aes(), size = 3, position = position_dodge(width = 0.6)) +
  xlab("Tanks") +
  ylab("cells per ml") +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank()) 
Fig.Exp2.cells

cells2.tank <- aov(cell.num ~Tank, data=cells.2)
anova(cells2.tank)
par(mfrow=c(3,2))
hist(cells2.tank$residuals)
boxplot(cells2.tank$residuals)
plot(cells2.tank)

cells2.trt <- aov(cell.num ~Treatment, data=cells.2)
anova(cells2.trt)
par(mfrow=c(3,2))
hist(cells2.trt$residuals)
boxplot(cells2.trt$residuals)
plot(cells2.trt)

##### JUVENILE GEODUCK SHELL SIZE #####

seed.size <- read.csv("Size_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

Initial <- subset(seed.size, Timepoint=="Initial", select = Date:Area)
Initial <- subset(seed.size, Component=="Exp1", select = Date:Area)
Initial$Ratio <- Initial$Length/Initial$Width

#calculating mean for Exposure 1 normalizations
T0.norms <- aggregate(Area ~ Day*Treatment, data=Initial, mean)

T0.norm.area.amb <- subset(T0.norms, Day=="1" & Treatment=="Ambient", select = Area)
T0.norm.area.med <- subset(T0.norms, Day=="1" & Treatment=="Medium", select = Area)
T0.norm.area.high <- subset(T0.norms, Day=="1" & Treatment=="High", select = Area)

T0.norms <- function(x) { 
  if(x == "Ambient") y <- T0.norm.area.amb
  if(x == "Medium") y <- T0.norm.area.med
  if(x == "High") y <- T0.norm.area.high
  return(y)
}

Initial$A.norm <- as.numeric(sapply(Initial$Treatment,T0.norms))
Initial$A.rel <- Initial$Area/Initial$A.norm
Init.avg.area <- aggregate(Area ~ Day*Treatment, data=Initial, mean)
Init.se.area <- aggregate(Area ~ Day*Treatment, data=Initial, std.error)
Init.avg.Len <- aggregate(Length ~ Day*Treatment, data=Initial, mean)
Init.se.Len <- aggregate(Length ~ Day*Treatment, data=Initial, std.error)
Init.avg.Wid <- aggregate(Width ~ Day*Treatment, data=Initial, mean)
Init.se.Wid <- aggregate(Width ~ Day*Treatment, data=Initial, std.error)
Init.avg.Ratio <- aggregate(Ratio ~ Day*Treatment, data=Initial, mean)
Init.se.Ratio<- aggregate(Ratio ~ Day*Treatment, data=Initial, std.error)
Init.avg.ARel <- aggregate(A.rel ~ Day*Treatment, data=Initial, mean)
Init.se.ARel<- aggregate(A.rel ~ Day*Treatment, data=Initial, std.error)
Init.Shell.size <- cbind(Init.avg.area, Init.se.area$Area, Init.avg.Len$Length, Init.se.Len$Length, Init.avg.Wid$Width, Init.se.Wid$Width,Init.avg.Ratio$Ratio, Init.se.Ratio$Ratio, Init.avg.ARel$A.rel, Init.se.ARel$A.rel)
colnames(Init.Shell.size) <- c("Day", "Treatment", "avg.area", "se.area","avg.len", "se.len", "avg.wid", "se.wid","avg.ratio", "se.ratio","avg.Arel", "se.Arel")

#calculating mean for Common Garden normalizations
E1.norm.area.amb <- subset(Init.Shell.size, Day=="23" & Treatment=="Ambient", select = avg.area)
E1.norm.area.med <- subset(Init.Shell.size, Day=="23" & Treatment=="Medium", select = avg.area)
E1.norm.area.high <- subset(Init.Shell.size, Day=="23" & Treatment=="High", select = avg.area)

CG <- subset(seed.size, Component=="Common.Garden", select = Date:Area)
CG$Ratio <- CG$Length/CG$Width
#ifelse(CG$Treatment == "Medium", CG$A.norm <- E1.norm.area.amb, NA)

E1.norms <- function(x) { 
  if(x == "Ambient") y <- E1.norm.area.amb
  if(x == "Medium") y <- E1.norm.area.med
  if(x == "High") y <- E1.norm.area.high
  return(y)
}

CG$A.norm <- as.numeric(sapply(CG$Treatment,E1.norms))
CG$A.rel <- CG$Area/CG$A.norm
CG.avg.area <- aggregate(Area ~ Day*Treatment, data=CG, mean)
CG.se.area <- aggregate(Area ~ Day*Treatment, data=CG, std.error)
CG.avg.Len <- aggregate(Length ~ Day*Treatment, data=CG, mean)
CG.se.Len <- aggregate(Length ~ Day*Treatment, data=CG, std.error)
CG.avg.Wid <- aggregate(Width ~ Day*Treatment, data=CG, mean)
CG.se.Wid <- aggregate(Width ~ Day*Treatment, data=CG, std.error)
CG.avg.Ratio <- aggregate(Ratio ~ Day*Treatment, data=CG, mean)
CG.se.Ratio<- aggregate(Ratio ~ Day*Treatment, data=CG, std.error)
CG.avg.ARel <- aggregate(A.rel ~ Day*Treatment, data=CG, mean)
CG.se.ARel<- aggregate(A.rel ~ Day*Treatment, data=CG, std.error)
CG.Shell.size <- cbind(CG.avg.area, CG.se.area$Area, CG.avg.Len$Length, CG.se.Len$Length, CG.avg.Wid$Width, CG.se.Wid$Width,CG.avg.Ratio$Ratio, CG.se.Ratio$Ratio, CG.avg.ARel$A.rel, CG.se.ARel$A.rel)
colnames(CG.Shell.size) <- c("Day", "Treatment", "avg.area", "se.area","avg.len", "se.len", "avg.wid", "se.wid","avg.ratio", "se.ratio","avg.Arel", "se.Arel")


#calculating mean for Exposure 2 normalizations
norm.area.amb <- subset(CG.Shell.size, Day=="135" & Treatment=="Ambient", select = avg.area)
norm.area.med <- subset(CG.Shell.size, Day=="135" & Treatment=="Medium", select = avg.area)
norm.area.high <- subset(CG.Shell.size, Day=="135" & Treatment=="High", select = avg.area)

ReExp <- subset(seed.size, Timepoint=="Reexposure", select = Date:Area)
ReExp$Ratio <- ReExp$Length/ReExp$Width

norms <- function(x) { 
  if(x == "Ambient") y <- norm.area.amb
  if(x == "Medium") y <- norm.area.med
  if(x == "High") y <- norm.area.high
  return(y)
}

ReExp$A.norm <- as.numeric(sapply(ReExp$Treatment,norms))
ReExp$A.rel <- ReExp$Area/ReExp$A.norm
ReExp.avg.area <- aggregate(Area ~ Day*Treatment*Secondary, data=ReExp, mean)
ReExp.se.area <- aggregate(Area ~ Day*Treatment*Secondary, data=ReExp, std.error)
ReExp.avg.Len <- aggregate(Length ~ Day*Treatment*Secondary, data=ReExp, mean)
ReExp.se.Len <- aggregate(Length ~ Day*Treatment*Secondary, data=ReExp, std.error)
ReExp.avg.Wid <- aggregate(Width ~ Day*Treatment*Secondary, data=ReExp, mean)
ReExp.se.Wid <- aggregate(Width ~ Day*Treatment*Secondary, data=ReExp, std.error)
ReExp.avg.Ratio <- aggregate(Ratio ~ Day*Treatment*Secondary, data=ReExp, mean)
ReExp.se.Ratio<- aggregate(Ratio ~ Day*Treatment*Secondary, data=ReExp, std.error)
ReExp.avg.ARel <- aggregate(A.rel ~ Day*Treatment*Secondary, data=ReExp, mean)
ReExp.se.ARel<- aggregate(A.rel ~ Day*Treatment*Secondary, data=ReExp, std.error)
ReExp.Shell.size <- cbind(ReExp.avg.area, ReExp.se.area$Area, ReExp.avg.Len$Length, ReExp.se.Len$Length, ReExp.avg.Wid$Width, ReExp.se.Wid$Width,ReExp.avg.Ratio$Ratio, ReExp.se.Ratio$Ratio, ReExp.avg.ARel$A.rel, ReExp.se.ARel$A.rel)
colnames(ReExp.Shell.size) <- c("Day", "Treatment", "Secondary", "avg.area", "se.area","avg.len", "se.len", "avg.wid", "se.wid","avg.ratio", "se.ratio","avg.Arel", "se.Arel")
ReExp.Shell.size$combos <- c("Ambient-Ambient", "Ambient-Ambient", "High-Ambient","High-Ambient","Medium-Ambient","Medium-Ambient","Ambient-Medium","Ambient-Medium","High-Medium","High-Medium","Medium-Medium","Medium-Medium")

# Significance testing for initial exposure
Init <- aov(log10(A.rel) ~ Treatment*Day, data=Initial)
Init.res <- anova(Init)
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
hist(Init$residuals)
boxplot(Init$residuals)
plot(Init)

# Significance testing for common garden
CG <- aov(log10(A.rel) ~ Treatment*Day, data=CG)
CG.res <- anova(CG)
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
hist(CG$residuals)
boxplot(CG$residuals)
plot(CG)

#posthoc <- TukeyHSD(x=X,conf.level=0.95)
#posthoc

Day1 <- subset(Initial, Day==1, select = Date:Area)
Day10 <- subset(Initial, Day==10, select = Date:Area)
Day23 <- subset(Initial, Day==23, select = Date:Area)
Day51 <- subset(CG, Day==51, select = Date:Area)
Day135 <- subset(CG, Day==135, select = Date:Area)

D1 <- aov(Area ~ Treatment, data=Day1)
anova(D1)
hist(D1$residuals)

D10 <- aov(Area ~ Treatment, data=Day10)
anova(D10)
hist(D10$residuals)

D23 <- aov(Area ~ Treatment, data=Day23)
anova(D23)
hist(D23$residuals)
posthoc.23 <- TukeyHSD(x=D23,conf.level=0.95)
posthoc.23

D51 <- aov(Area ~ Treatment, data=Day51)
anova(D51)
hist(D51$residuals)

D135 <- aov(Area ~ Treatment, data=Day135)
anova(D135)
hist(D135$residuals)
posthoc.135 <- TukeyHSD(x=D135,conf.level=0.95)
posthoc.135

#rect <- data.frame(xmin=3.1, xmax=Inf, ymin=-Inf, ymax=Inf) #identify background shading placement

Fig.Exp1.size <- ggplot(Init.Shell.size, aes(x=Day, y=avg.Arel, group=Treatment)) + 
  geom_errorbar(aes(ymin=Init.Shell.size$avg.Arel-Init.Shell.size$se.Arel, ymax=Init.Shell.size$avg.Arel+Init.Shell.size$se.Arel), colour="black", width=.1, position = position_dodge(width = 0.6)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.6)) +   
  geom_point(aes(shape=Treatment), size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(breaks=seq(0,160,10)) +
  xlab("Days") +
  ylab(expression(paste("Relative Size"))) +
  ylim(0.5,3.5) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),  #remove legend background
        legend.position=c(.2, .7)) + #set legend location
  ggtitle("A) Initial Exposure") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.Exp1.size


Fig.CG.size <- ggplot(CG.Shell.size, aes(x=Day, y=avg.Arel, group=Treatment)) + 
  geom_errorbar(aes(ymin=CG.Shell.size$avg.Arel-CG.Shell.size$se.Arel, ymax=CG.Shell.size$avg.Arel+CG.Shell.size$se.Arel), colour="black", width=.1, position = position_dodge(width = 0.6)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.6)) +   
  geom_point(aes(shape=Treatment), size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(breaks=seq(0,160,10)) +
  xlab("Days") +
  ylab(expression(paste("Relative Size"))) +
  ylim(0.5,3.5) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.position='none') + #remove legend background
  ggtitle("B) Common Garden") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))
Fig.CG.size

# Significance testing for secondary exposure
Rexp <- aov(log10(Area) ~ Treatment*Secondary*Day, data=ReExp)
Rexp.res <- anova(Rexp)
#Check Assumptions
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
hist(Rexp$residuals)
boxplot(Rexp$residuals)
plot(Rexp)

# Rexp.posthoc <- TukeyHSD(x=Rexp,conf.level=0.95)
# Rexp.posthoc

Day145 <- subset(ReExp.Shell.size, Day==145, select = Day:se.Arel)
Day158 <- subset(ReExp.Shell.size, Day==158, select = Day:se.Arel)

Fig.Exp2.D10.size <- ggplot(Day145, aes(x=Secondary, y=avg.Arel, group=Treatment)) + 
  geom_errorbar(aes(ymin=Day145$avg.Arel-Day145$se.Arel, ymax=Day145$avg.Arel+Day145$se.Arel), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 3, position = position_dodge(width = 0.05)) +
  xlab("Secondary Treatment") +
  ylab("Relative Size") +
  ylim(0.5,1.5) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none') + #remove legend background
  ggtitle("C) Secondary Exposure Day 10") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))

Fig.Exp2.D10.size


Fig.Exp2.D23.size <- ggplot(Day158, aes(x=Secondary, y=avg.Arel, group=Treatment)) + 
  geom_errorbar(aes(ymin=Day158$avg.Arel-Day158$se.Arel, ymax=Day158$avg.Arel+Day158$se.Arel), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 3, position = position_dodge(width = 0.05)) +
  xlab("Secondary Treatment") +
  ylab("Relative Size") +
  ylim(0.5,1.5) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position='none') + #remove legend background
  ggtitle("D) Secondary Exposure Day 23") +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0))

Fig.Exp2.D23.size


#all size over time
All.avg.area <- aggregate(Area ~ Day*Treatment*Secondary, data=seed.size, mean)
All.se.area <- aggregate(Area ~ Day*Treatment*Secondary, data=seed.size, std.error)

All <- cbind(All.avg.area, All.se.area$Area)
colnames(All) <- c("Day",  "Treatment",	"Secondary",	"mean", "se")
#All$grps <- c("pH8.0", "pH8.0","pH8.0","pH8.0","pH8.0","pH8.0","pH7.4", "pH7.4", "pH7.4", "pH7.4", "pH7.4", "pH7.4", "None", "None", "None", "None", "None","None","None", "None","None", "None", "None","None","None","None","None")

# #grouped plotting over time
# Fig.Time <- ggplot(All, aes(x=Day, y=mean, group=Treatment)) + 
#   geom_errorbar(aes(ymin=All$mean-All$se, ymax=All$mean+All$se), colour="black", width=.1, position = position_dodge(width = 0.6)) +
#   geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.6)) +   
#   geom_point(aes(shape=Treatment), size = 3, position = position_dodge(width = 0.6)) +
#   #scale_shape_manual(values = c(16, 21, 15, 22, 17,24)) +
#   scale_x_continuous(breaks=seq(0,160,10)) +
#   xlab("Days") +
#   ylab(expression(paste("Relative Size"))) +
#   #ylim(0.5,3.5) +
#   theme_bw() + #Set the background color
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
#         axis.line = element_line(color = 'black'), #Set the axes color
#         panel.border = element_blank(), #Set the border
#         panel.grid.major = element_blank(), #Set the major gridlines
#         panel.grid.minor = element_blank(), #Set the minor gridlines
#         plot.background=element_blank(),  #Set the plot background
#         legend.key = element_blank(),  #remove legend background
#         legend.position=c(.2, .6)) + #set legend location
#   ggtitle("E) Exposure 2 over time") +
#   theme(plot.title = element_text(face = 'bold', 
#                                   size = 12, 
#                                   hjust = 0))
# Fig.Time


##### CAPTURE STATISTICAL OUTPUT AND FIGURES TO FILE #####
setwd(file.path(mainDir, 'Output'))
capture.output(Init.res, CG.res, Rexp.res, file="Geoduck_Juvenile_Statistical_Results.txt")

#Capture Figures to File
Figure1.Size <- arrangeGrob(Fig.Exp1.size,Fig.CG.size, Fig.Exp2.D10.size, Fig.Exp2.D23.size, ncol=2)
ggsave(file="Geoduck_Size.pdf", Figure1.Size, width = 6.5, height = 8, units = c("in"))

FigureS1 <- arrangeGrob(Fig.E1.daily.pH, Fig.ICG.pH, Fig.OCG.pH, Fig.E2.daily.pH, ncol=4)
ggsave(file="Fig.S1.pH.pdf", FigureS1, width = 11, height = 4, units = c("in"))

FigureS2 <- arrangeGrob( Fig.E1.daily.Temp, Fig.ICG.Temp, Fig.OCG.Temp, Fig.E2.daily.Temp, ncol=4)
ggsave(file="Fig.S1.Temp.pdf", FigureS2, width = 11, height = 4, units = c("in"))

setwd(file.path(mainDir, 'Data'))
