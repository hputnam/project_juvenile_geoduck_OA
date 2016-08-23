#Effects of Ocean Acidification on Phenotype and DNA Methylation in Juvenile Geoduck
#Data published in
#Title:
#Contact: Hollie Putnam hollieputnam@gmail.com ; Steven Roberts sr320@u.washington.edu
#Supported by: NOAA OA
#last modified 20160822 H Putnam
#See Readme file for project details 
#See metadata file for details on data files and equipment 

rm(list=ls()) # removes all prior objects

#Read in required libraries
library(seacarb) 
library(reshape) 
library(plotrix) 
library(ggplot2) 
library(plyr)
library(gridExtra)
library(multcompView)
library(lsmeans)

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
Flow1$Rate <- (((Flow1$Rate1 + Flow1$Rate2)/2)*6)/1000*60 #calculate flow rate in liters per hour
Flow1.Rate <- cat("Flow (mean±sem, n) =", mean(Flow1$Rate),  sd(Flow1$Rate)/sqrt(length(Flow1$Rate)), length(Flow1$Rate)) #calculate, concatenate and print average±sem flow rate

FL1 <- aov(Rate ~ Treatment, data=Flow1) #test for differences in flow rate between treatments
FL1.res <- anova(FL1) #display anova results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(FL1$residuals) #plot histogram of residuals
boxplot(FL1$residuals) #plot boxplot of residuals
plot(FL1) #display residuals versus fitter, normal QQ plot, leverage plot

Flow2 <- read.csv("Flow2_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
#flow values measured in ml/15 sec
Flow2$Rate <- ((Flow2$Rate1)*4)/1000*60 #calculate flow rate in liters per hour
Flow2.Rate <- cat("Flow (mean±sem, n) =", mean(Flow2$Rate),  sd(Flow2$Rate)/sqrt(length(Flow2$Rate)), length(Flow2$Rate)) #calculate, concatenate and print average±sem flow rate

FL2 <- aov(Rate ~ Treatment, data=Flow2) #test for differences in flow rate between treatments
FL2.res <-anova(FL2) #display anova results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(FL2$residuals) #plot histogram of residuals
boxplot(FL2$residuals) #plot boxplot of residuals
plot(FL2) #display residuals versus fitter, normal QQ plot, leverage plot

#Feeding
#Exposure1
#Cell_Counts_Juvenile_Geoduck_Exp1.csv
cells.1 <- read.csv("Cell_Counts_Juvenile_Geoduck_Exp1.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
cells.1$Avg.Cells <- rowMeans(cells.1[,c("Count1",  "Count2")], na.rm = TRUE) #calculate average of counts
cells.1$cell.num <- cells.1$Avg.Cells/cells.1$Volume.Counted #calculate density

avg.cells.tank.1 <- do.call(data.frame,aggregate(cell.num ~ Tank, data = cells.1, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each Tank
avg.cells.trt.1 <- do.call(data.frame,aggregate(cell.num ~ Treatment, data = cells.1, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each Treatment
colnames(avg.cells.tank.1) <- c("Tank", "mean", "se") #rename columns 
colnames(avg.cells.trt.1) <- c("Treatment", "mean", "se") #rename columns 

cells1.tank <- aov(cell.num ~Tank, data=cells.1) #test the hypothesis feeding does not differ between tanks
cells1.tank.res <-anova(cells1.tank) #display anova results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(cells1.tank$residuals) #plot histogram of residuals
boxplot(cells1.tank$residuals) #plot boxplot of residuals
plot(cells1.tank) #display residuals versus fitter, normal QQ plot, leverage plot

cells1.trt <- aov(cell.num ~Treatment, data=cells.1) #test the hypothesis feeding does not differ between treatments
cells1.trt.res <-anova(cells1.trt) #display anova results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(cells1.trt$residuals) #plot histogram of residuals
boxplot(cells1.trt$residuals) #plot boxplot of residuals
plot(cells1.trt) #display residuals versus fitter, normal QQ plot, leverage plot

#Exposure2
#Cell_Counts_Juvenile_Geoduck_Exp2.csv
cells.2 <- read.csv("Cell_Counts_Juvenile_Geoduck_Exp2.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
cells.2$cell.num <- rowMeans(cells.2[,c("Count1",  "Count2", "Count3")], na.rm = TRUE) #calculate average of counts
avg.cells.tank.2 <- do.call(data.frame,aggregate(cell.num ~ Tank, data = cells.2, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each Tank
avg.cells.trt.2 <- do.call(data.frame,aggregate(cell.num ~ Treatment, data = cells.2, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each Treatment
colnames(avg.cells.tank.2) <- c("Tank", "mean", "se") #rename columns 
colnames(avg.cells.trt.2) <- c("Treatment", "mean", "se") #rename columns 

cells2.tank <- aov(cell.num ~Tank, data=cells.2) #test the hypothesis feeding does not differ between tanks
cells2.tank.res <-anova(cells2.tank) #display anova results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(cells2.tank$residuals) #plot histogram of residuals
boxplot(cells2.tank$residuals) #plot boxplot of residuals
plot(cells2.tank) #display residuals versus fitter, normal QQ plot, leverage plot

cells2.trt <- aov(cell.num ~Treatment, data=cells.2) #test the hypothesis feeding does not differ between treatments
cells2.trt.res <-anova(cells2.trt) #display anova results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(cells2.trt$residuals) #plot histogram of residuals
boxplot(cells2.trt$residuals) #plot boxplot of residuals
plot(cells2.trt) #display residuals versus fitter, normal QQ plot, leverage plot

##### CONTINUOUS EXPERIMENTAL PH DATA #####
#Avtech_pH_data.csv
pH <- read.csv("Avtech_pH_data.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
pH$Date.Time <-as.POSIXct(pH$Date.Time, format="%m/%d/%y %H:%M") #convert date format
pH$Date <- as.Date(pH$Date.Time) #convert Date only
pH$Time <- format(as.POSIXct(pH$Date.Time) ,format = "%H:%M:%S") #convert Time only

pH.low <- do.call(data.frame,aggregate(Low ~ Date, data = pH, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Day
pH.med <- do.call(data.frame,aggregate(Medium ~ Date, data = pH, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Day
pH.amb <- do.call(data.frame,aggregate(Ambient ~ Date, data = pH, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Day
pH.low$Treatment <- "Low" #Add treatment Information
colnames(pH.low) <- c("Date", "mean", "se", "Treatment") #rename columns to generic format
pH.med$Treatment <- "Medium" #Add treatment Information
colnames(pH.med) <- c("Date", "mean", "se", "Treatment") #rename columns to generic format
pH.amb$Treatment <- "Ambient" #Add treatment Information
colnames(pH.amb) <- c("Date", "mean", "se", "Treatment") #rename columns to generic format
daily.pH <- rbind(pH.amb, pH.med, pH.low) #bind treatment data 
daily.pH #view data

# Plot daily averages of pH data for the complete experiment
All.pH <- ggplot(daily.pH, aes(x=Date, y=mean, group=Treatment)) + #set up plot information
  geom_errorbar(aes(ymin=daily.pH$mean-daily.pH$se, ymax=daily.pH$mean+daily.pH$se), colour="black", width=.1, position = position_dodge(width = 0.05)) + #add standard error bars about the mean
  geom_point(aes(shape=Treatment), size = 2, position = position_dodge(width = 0.05)) + #include points in the shape of the treatments
  annotate("text", x=as.Date("2016-03-25"), y=7.95, label = "No Data") + #add text to the graphic where data are missing
  annotate("text", x=as.Date("2016-06-15"), y=7.95, label = "No Data") + #add text to the graphic where data are missing
  xlab("Time") + #label x axis
  ylab("pH Total Scale") + # label y axis
  ylim(6.8,8.1) + # set y axis scale
  geom_vline(xintercept = 16899, linetype="dotted", color = "gray", size=1) + #add vertical line
  geom_vline(xintercept = 16928, linetype="dashed", color = "gray", size=1) + #add vertical line
  geom_vline(xintercept = 17013, linetype="dotted", color = "gray", size=1) + #add vertical line
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position=c(.65, .25), #set legend location
        legend.text = element_text(size = 8), #set the legend text size
        legend.key = element_blank(), #remove the legend background
        legend.title = element_text(size=8, face="bold")) + #set legend title attributes
  ggtitle("A) Experimental pH") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
All.pH #view plot

##### CONTINUOUS EXPERIMENTAL TEMPERATURE DATA #####
#Avtech_temp_data.csv
temp <- read.csv("Avtech_temp_data.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
temp$Date.Time <-as.POSIXct(temp$Date.Time, format="%m/%d/%y %H:%M") #convert date format
temp$Date <- as.Date(temp$Date.Time) #convert Date only
temp$Time <- format(as.POSIXct(temp$Date.Time) ,format = "%H:%M:%S") #convert time only

temp.low <- do.call(data.frame,aggregate(Low ~ Date, data = temp, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Day
temp.med <- do.call(data.frame,aggregate(Medium ~ Date, data = temp, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Day
temp.amb <- do.call(data.frame,aggregate(Ambient ~ Date, data = temp, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem of each treatment by Day
temp.low$Treatment <- "Low" #Add treatment Information
colnames(temp.low) <- c("Date", "mean", "se", "Treatment") #rename columns to generic format
temp.med$Treatment <- "Medium" #Add treatment Information
colnames(temp.med) <- c("Date", "mean", "se", "Treatment") #rename columns to generic format
temp.amb$Treatment <- "Ambient" #Add treatment Information
colnames(temp.amb) <- c("Date", "mean", "se", "Treatment") #rename columns to generic format
daily.temp <- rbind(temp.amb, temp.med, temp.low) #bind treatment data 
daily.temp #view data

All.temp <- ggplot(daily.temp, aes(x=Date, y=mean, group=Treatment)) + #set up plot information
  geom_errorbar(aes(ymin=daily.temp$mean-daily.temp$se, ymax=daily.temp$mean+daily.temp$se), colour="black", width=.1, position = position_dodge(width = 0.05)) + #add standard error bars about the mean
  geom_point(aes(shape=Treatment), size = 2, position = position_dodge(width = 0.05)) + #include points in the shape of the treatments
  annotate("text", x=as.Date("2016-05-25"), y=14, label = "No Data") + #add text to the graphic where data are missing
  xlab("Time") + #label x axis
  ylab("Temperature °C") + #label y axis
  ylim(0,20) + # set y axis scale
  geom_vline(xintercept = 16899, linetype="dotted", color = "gray", size=1) + #add vertical line
  geom_vline(xintercept = 16928, linetype="dashed", color = "gray", size=1) + #add vertical line
  geom_vline(xintercept = 17013, linetype="dotted", color = "gray", size=1) + #add vertical line
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(), #Set the plot background
        legend.position = "none") + #remove legend
  ggtitle("B) Experimental Temperature") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0)) #set title attributes
All.temp #view plot

##### SEAWATER CHEMISTRY ANALYSIS FOR DISCRETE MEASUREMENTS #####

#pH Tris Calibration Curves
#Data to calculate conversion equations from mV to total scale using tris standard for pH probe
path <-("/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/Data/pH_Calibration_Files/") #set path to calibration file folder
file.names<-list.files(path = path, pattern = "csv$") #list all the file names with csv 
pH.cals <- data.frame(matrix(NA, nrow=length(file.names), ncol=4, dimnames=list(file.names,c("Date", "Intercept", "Slope","R2")))) #generate an empty 3 column dataframe with specific column names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Calib.Data <-read.table(file.path(path,file.names[i]), header=TRUE, sep=",", na.string="NA", as.is=TRUE) #reads in the data files
  model <-lm(mVTris ~ TTris, data=Calib.Data) #runs a linear regression of mV as a function of temperature
  coe <- coef(model) #extracts the coeffecients
  R <- summary(model)$r.squared #extracts the R2
  pH.cals[i,2:3] <- coe #inserts coef in the dataframe
  pH.cals[i,4] <- R #inserts R2 in the dataframe
  pH.cals[i,1] <- substr(file.names[i],1,8) #stores the file name in the Date column
}

colnames(pH.cals) <- c("Calib.Date",  "Intercept",  "Slope", "R2") #names the columns of the dataframe
pH.cals #view data

#Prepare input file for seacarb
SW.chem <- read.csv("SW_Chem_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
SW.chem <- merge(pH.cals, SW.chem, by="Calib.Date") #merge pH calibrations with Seawater chemistry data
R <- 8.31447215 #gas constant in J mol-1 K-1 
F <-96485.339924 #Faraday constant in coulombs mol-1
mvTris <- SW.chem$Temperature*SW.chem$Slope+SW.chem$Intercept #calculate the mV of the tris standard using the temperature mv relationships in the measured standard curves 
STris<-27.5 #salinity of the Tris from local batch made for Manchester WA salinity conditions in Spring 2016. Batch date = 20160214
phTris<- (11911.08-18.2499*STris-0.039336*STris^2)*(1/(SW.chem$Temperature+273.15))-366.27059+ 0.53993607*STris+0.00016329*STris^2+(64.52243-0.084041*STris)*log(SW.chem$Temperature+273.15)-0.11149858*(SW.chem$Temperature+273.15) #calculate the pH of the tris (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem$pH.Total<-phTris+(mvTris/1000-SW.chem$pH.MV/1000)/(R*(SW.chem$Temperature+273.15)*log(10)/F) #calculate the pH on the total scale (Dickson A. G., Sabine C. L. and Christian J. R., SOP 6a)
SW.chem <- na.omit(SW.chem) #remove NAs

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

#Calculate descriptive stats for seawater chemistry by Exposure
ccarb <- melt(carb.output[,c(1,3,4,5:14)], id.vars=c("Exposure", "Treatment", "Tank")) #reshape data into long format
Exp1 <-subset(ccarb, Exposure == "Exposure1") #separate out exposure 1 for all data
Exp2 <-subset(ccarb, Exposure == "Exposure2") #separate out exposure 2 for all data

SWC.Tanks <- ddply(ccarb, c("Exposure", "Tank", "variable"), summarise, #apply functions to sewater chem data
      N = length(na.omit(value)), #count the sample size removing NA
      mean = mean(value), #calculate average 
      sem = sd(value)/sqrt(N)) #calcualte the standard error of the mean

#Test for tank and treatment differences in Temperature and Total Alkalinity in Exposure 1
Exp1.Temp <-subset(Exp1, variable=="Temperature") #separate out exposure 1 for all data
temp1.tank <- aov(value ~Tank, data=Exp1.Temp) #test the hypothesis there is no difference in temperature between tanks
temp1.tank.res <-anova(temp1.tank) #view results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(temp1.tank$residuals) #plot histogram of residuals
boxplot(temp1.tank$residuals) #plot boxplot of residuals
plot(temp1.tank) #display residuals versus fitter, normal QQ plot, leverage plot

temp1.trt <- aov(value ~Treatment, data=Exp1.Temp) #test the hypothesis there is no difference in temperature between treatments
temp1.trt.res <- anova(temp1.trt) #statistical results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(temp1.trt$residuals) #plot histogram of residuals
boxplot(temp1.trt$residuals) #plot boxplot of residuals
plot(temp1.trt) #display residuals versus fitter, normal QQ plot, leverage plot

Exp1.TA <-subset(Exp1, variable=="TA") #separate out exposure 1 for all data
TA1.tank <- aov(value ~Tank, data=Exp1.TA) #test the hypothesis there is no difference in total alkalinity between tanks
TA1.tank.res <- anova(temp1.trt) #statistical results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(TA1.tank$residuals) #plot histogram of residuals
boxplot(TA1.tank$residuals) #plot boxplot of residuals
plot(TA1.tank) #display residuals versus fitter, normal QQ plot, leverage plot

TA1.trt <- aov(value ~Treatment, data=Exp1.TA) #test the hypothesis there is no difference in total alkalinity between treatments
TA1.trt.res <- anova(temp1.trt) #statistical results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(TA1.trt$residuals) #plot histogram of residuals
boxplot(TA1.trt$residuals) #plot boxplot of residuals
plot(TA1.trt) #display residuals versus fitter, normal QQ plot, leverage plot

#Test for tank and treatment differences in Temperature and Total Alkalinity in Exposure 2
Exp2.Temp <-subset(Exp2, variable=="Temperature") #separate out exposure 2 for all data
temp2.tank <- aov(value ~Tank, data=Exp2.Temp) #test the hypothesis there is no difference in temperature between tanks
temp2.tank.res <-anova(temp2.tank) #view results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(temp2.tank$residuals) #plot histogram of residuals
boxplot(temp2.tank$residuals) #plot boxplot of residuals
plot(temp2.tank) #display residuals versus fitter, normal QQ plot, leverage plot

temp2.trt <- aov(value ~Treatment, data=Exp2.Temp) #test the hypothesis there is no difference in temperature between treatments
temp2.trt.res <- anova(temp2.trt) #statistical results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(temp2.trt$residuals) #plot histogram of residuals
boxplot(temp2.trt$residuals) #plot boxplot of residuals
plot(temp2.trt) #display residuals versus fitter, normal QQ plot, leverage plot

Exp2.TA <-subset(Exp2, variable=="TA") #separate out exposure 2 for all data
TA2.tank <- aov(value ~Tank, data=Exp2.TA) #test the hypothesis there is no difference in total alkalinity between tanks
TA2.tank.res <- anova(temp2.trt) #statistical results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(TA2.tank$residuals) #plot histogram of residuals
boxplot(TA2.tank$residuals) #plot boxplot of residuals
plot(TA2.tank) #display residuals versus fitter, normal QQ plot, leverage plot

TA2.trt <- aov(value ~Treatment, data=Exp2.TA) #test the hypothesis there is no difference in total alkalinity between treatments
TA2.trt.res <- anova(temp2.trt) #statistical results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(TA2.trt$residuals) #plot histogram of residuals
boxplot(TA2.trt$residuals) #plot boxplot of residuals
plot(TA2.trt) #display residuals versus fitter, normal QQ plot, leverage plot

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

##### JUVENILE GEODUCK SHELL SIZE #####

seed.size <- read.csv("Size_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA

#Exposure 2 Experimental Component
Initial <- subset(seed.size, Component=="Exp1", select = Date:Area) #subset dataframe to exposure 1 data
Initial$Ratio <- Initial$Length/Initial$Width #calculate the length:width of the shell

#calculating shell area means for Exposure 1 normalizations
E1.norms <- aggregate(Area ~ Day*Treatment, data=Initial, mean) #calculate mean area by Day and Treatment
E1.norm.area.amb <- subset(E1.norms, Day=="Day1" & Treatment=="Ambient", select = Area) #subset mean area for day 1 ambient condition
E1.norm.area.med <- subset(E1.norms, Day=="Day1" & Treatment=="Medium", select = Area) #subset mean area for day 1 medium condition
E1.norm.area.low <- subset(E1.norms, Day=="Day1" & Treatment=="Low", select = Area) #subset mean area for day 1 low condition

E1.norms <- function(x) {  #write function
  if(x == "Ambient") y <- E1.norm.area.amb #if Treatment equals Ambient assign day 1 Ambient mean as normalization factor
  if(x == "Medium") y <- E1.norm.area.med #if Treatment equals Medium assign day 1 Medium mean as normalization factor
  if(x == "Low") y <- E1.norm.area.low #if Treatment equals Low assign day 1 Low mean as normalization factor
  return(y) #return result
}

Initial$A.norm <- as.numeric(sapply(Initial$Treatment,E1.norms)) #add a column with the normalization values 
Initial$A.rel <- Initial$Area/Initial$A.norm #normalize the area to be relative size
Init.area <- do.call(data.frame,aggregate(Area ~ Day*Treatment, data = Initial, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
Init.Len <- do.call(data.frame,aggregate(Length ~ Day*Treatment, data = Initial, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
Init.Wid <- do.call(data.frame,aggregate(Width ~ Day*Treatment, data = Initial, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
Init.Ratio <- do.call(data.frame,aggregate(Ratio ~ Day*Treatment, data = Initial, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
Init.ARel <- do.call(data.frame,aggregate(A.rel ~ Day*Treatment, data = Initial, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
Init.Shell.size <- cbind(Init.area, Init.Len$Length.mean, Init.Len$Length.se, Init.Wid$Width.mean, Init.Wid$Width.se, Init.Ratio$Ratio.mean, Init.Ratio$Ratio.se, Init.ARel$A.rel.mean, Init.ARel$A.rel.se) #bind columns into dataframe
colnames(Init.Shell.size) <- c("Day", "Treatment", "avg.area", "se.area","avg.len", "se.len", "avg.wid", "se.wid","avg.ratio", "se.ratio","avg.Arel", "se.Arel") #rename columns

# Significance testing for initial exposure
Init <- aov(log10(A.rel) ~ Treatment*Day, data=Initial) #test the hypothesis relative size does not differ between treatments and time
Init.res <-anova(Init) #display anova results
par(mfrow=c(3,2)) #set plotting configuration
par(mar=c(1,1,1,1)) #set margins for plots
hist(Init$residuals) #plot histogram of residuals
boxplot(Init$residuals) #plot boxplot of residuals
plot(Init) #display residuals versus fitter, normal QQ plot, leverage plot

#PostHoc Tukey adjustment comparison of least-squares means
E1.lsm <- lsmeans(Init, ~ Treatment*Day, adjust="tukey") #compute least-squares means for Treatment*Day from ANOVA model
E1.pairs <- cld(E1.lsm, alpha=.05, Letters=letters) #list pairwise tests and letter display
E1.pairs #view results

#Exposure 1 Plotting
Init.Shell.size$Day.Num <- c(1,10,23,1,10,23,1,10,23) #set days as numeric for plotting

Fig.Exp1.size <- ggplot(Init.Shell.size, aes(x=Day.Num, y=avg.Arel, group=Treatment)) + 
  geom_errorbar(aes(ymin=Init.Shell.size$avg.Arel-Init.Shell.size$se.Arel, ymax=Init.Shell.size$avg.Arel+Init.Shell.size$se.Arel), colour="black", width=.1, position = position_dodge(width = 0.6)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.6)) +   
  geom_point(aes(shape=Treatment), size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(breaks=seq(0,160,10)) +
  annotate("text", x=9, y=0.73, label = "a", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=1, y=1.2, label = "abc", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=11, y=1, label = "ab", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=c(9,22.5), y=c(1.2,0.95), label = "b", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=22.5, y=1.25, label = "bc", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=22, y=1.6, label = "c", size = 3) + #add text to the graphic for posthoc letters
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

#Common Garden Experimental Component
CG.data <- subset(seed.size, Component=="Common.Garden", select = Date:Area) #subset dataframe to common garden data
CG.data$Ratio <- CG.data$Length/CG.data$Width #calculate the length:width of the shell

#calculating shell area means for Common Garden normalizations
CG.norm.area.amb <- subset(Init.Shell.size, Day=="Day23" & Treatment=="Ambient", select = avg.area) #subset mean area for day 23 ambient condition
CG.norm.area.med <- subset(Init.Shell.size, Day=="Day23" & Treatment=="Medium", select = avg.area) #subset mean area for day 23 medium condition
CG.norm.area.low <- subset(Init.Shell.size, Day=="Day23" & Treatment=="Low", select = avg.area) #subset mean area for day 23 high condition

CG.norms <- function(x) { 
  if(x == "Ambient") y <- CG.norm.area.amb #if Treatment equals Ambient assign day 23 Ambient mean as normalization factor
  if(x == "Medium") y <- CG.norm.area.med #if Treatment equals Medium assign day 23 Medium mean as normalization factor
  if(x == "Low") y <- CG.norm.area.low #if Treatment equals Low assign day 23 Low mean as normalization factor
  return(y)
}

CG.data$A.norm <- as.numeric(sapply(CG.data$Treatment,CG.norms)) #add a column with the normalization values
CG.data$A.rel <- CG.data$Area/CG.data$A.norm #normalize the area to be relative size
CG.area <- do.call(data.frame,aggregate(Area ~ Day*Treatment, data = CG.data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
CG.Len <- do.call(data.frame,aggregate(Length ~ Day*Treatment, data = CG.data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
CG.Wid <- do.call(data.frame,aggregate(Width ~ Day*Treatment, data = CG.data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
CG.Ratio <- do.call(data.frame,aggregate(Ratio ~ Day*Treatment, data = CG.data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
CG.ARel <- do.call(data.frame,aggregate(A.rel ~ Day*Treatment, data = CG.data, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
CG.Shell.size <- cbind(CG.area, CG.Len$Length.mean, CG.Len$Length.se, CG.Wid$Width.mean, CG.Wid$Width.se, CG.Ratio$Ratio.mean, CG.Ratio$Ratio.se, CG.ARel$A.rel.mean, CG.ARel$A.rel.se) #bind columns into dataframe
colnames(CG.Shell.size) <- c("Day", "Treatment", "avg.area", "se.area","avg.len", "se.len", "avg.wid", "se.wid","avg.ratio", "se.ratio","avg.Arel", "se.Arel") #rename columns

# Significance testing for common garden
CG <- aov(log10(A.rel) ~ Treatment*Day, data=CG.data)
CG.res <- anova(CG)
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
hist(CG$residuals)
boxplot(CG$residuals)
plot(CG)

#PostHoc Tukey adjustment comparison of least-squares means
CG.lsm <- lsmeans(CG, ~ Treatment*Day, adjust="tukey") #compute least-squares means for Treatment*Day from ANOVA model
CG.pairs <- cld(CG.lsm, alpha=.05, Letters=letters) #list pairwise tests and letter display
CG.pairs #view results

#Common Garden Plotting
CG.Shell.size$Day.Num <- c(135,51,135,51,135,51) #set days as numeric for plotting

Fig.CG.size <- ggplot(CG.Shell.size, aes(x=Day.Num, y=avg.Arel, group=Treatment)) + 
  geom_errorbar(aes(ymin=CG.Shell.size$avg.Arel-CG.Shell.size$se.Arel, ymax=CG.Shell.size$avg.Arel+CG.Shell.size$se.Arel), colour="black", width=.1, position = position_dodge(width = 0.6)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.6)) +   
  geom_point(aes(shape=Treatment), size = 3, position = position_dodge(width = 0.6)) +
  scale_x_continuous(breaks=seq(0,160,10)) +
  annotate("text", x=43, y=1.85, label = "a", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=132, y=2.15, label = "b", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=c(43,43), y=c(2.45,2.2), label = "ab", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=c(141,141), y=c(3.15,3.4), label = "c", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=35, y=0.5, label = ".", size = 1) + #add text to the graphic for posthoc letters
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

#Exposure 2 Experimental Component
ReExp <- subset(seed.size, Timepoint=="Reexposure", select = Date:Area) #subset dataframe to exposure 2 data
ReExp$Ratio <- ReExp$Length/ReExp$Width #calculate the length:width of the shell

#calculating mean for Exposure 2 normalizations
E2.norm.area.amb <- subset(CG.Shell.size, Day=="Day135" & Treatment=="Ambient", select = avg.area) #subset mean area for day 135 ambient condition
E2.norm.area.med <- subset(CG.Shell.size, Day=="Day135" & Treatment=="Medium", select = avg.area) #subset mean area for day 135 medium condition
E2.norm.area.low<- subset(CG.Shell.size, Day=="Day135" & Treatment=="Low", select = avg.area) #subset mean area for day 135 low condition

E2.norms <- function(x) { 
  if(x == "Ambient") y <- E2.norm.area.amb #if Treatment equals Ambient assign day 135 Ambient mean as normalization factor
  if(x == "Medium") y <- E2.norm.area.med #if Treatment equals Ambient assign day 135 Medium mean as normalization factor
  if(x == "Low") y <- E2.norm.area.low #if Treatment equals Ambient assign day 135 Low mean as normalization factor
  return(y)
}

ReExp$A.norm <- as.numeric(sapply(ReExp$Treatment,E2.norms)) #add a column with the normalization values
ReExp$A.rel <- ReExp$Area/ReExp$A.norm #normalize the area to be relative size
ReExp.area <- do.call(data.frame,aggregate(Area ~ Treatment*Secondary*Day, data = ReExp, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
ReExp.Len <- do.call(data.frame,aggregate(Length ~ Treatment*Secondary*Day, data = ReExp, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
ReExp.Wid <- do.call(data.frame,aggregate(Width ~ Treatment*Secondary*Day, data = ReExp, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
ReExp.Ratio <- do.call(data.frame,aggregate(Ratio ~ Treatment*Secondary*Day, data = ReExp, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
ReExp.ARel <- do.call(data.frame,aggregate(A.rel ~ Treatment*Secondary*Day, data = ReExp, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
ReExp.Shell.size <- cbind(ReExp.area, ReExp.Len$Length.mean, ReExp.Len$Length.se, ReExp.Wid$Width.mean, ReExp.Wid$Width.se, ReExp.Ratio$Ratio.mean, ReExp.Ratio$Ratio.se, ReExp.ARel$A.rel.mean, ReExp.ARel$A.rel.se) #bind columns into dataframe
colnames(ReExp.Shell.size) <- c("Treatment", "Secondary", "Day", "avg.area", "se.area","avg.len", "se.len", "avg.wid", "se.wid","avg.ratio", "se.ratio","avg.Arel", "se.Arel") #rename columns

# Significance testing for secondary exposure
Rexps <- aov(log10(A.rel) ~ Treatment*Secondary*Day, data=ReExp)
Rexp.res <- anova(Rexps)
par(mfrow=c(3,2))
par(mar=c(1,1,1,1))
hist(Rexps$residuals)
boxplot(Rexps$residuals)
plot(Rexps)

#PostHoc
Day145 <- subset(ReExp.Shell.size, Day=="Day145", select = Treatment:se.Arel)
Day158 <- subset(ReExp.Shell.size, Day=="Day158", select = Treatment:se.Arel)

E2.lsm <- lsmeans(Rexps, ~ Treatment*Secondary*Day, adjust="tukey") #compute least-squares means for Treatment*Day from ANOVA model
E2.pairs <- cld(E2.lsm, alpha=.05, Letters=letters) #list pairwise tests and letter display
E2.pairs #view results

#Exposure2 Plotting
Fig.Exp2.D10.size <- ggplot(Day145, aes(x=Secondary, y=avg.Arel, group=Treatment)) + 
  geom_errorbar(aes(ymin=Day145$avg.Arel-Day145$se.Arel, ymax=Day145$avg.Arel+Day145$se.Arel), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 3, position = position_dodge(width = 0.05)) +
  annotate("text", x=0.85, y=1.3, label = "b", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=c(0.85,0.85,2.2,2.2,2.25), y=c(0.9,0.95,0.88, 0.92,0.97), label = "ab", size = 3) + #add text to the graphic for posthoc letters
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
  annotate("text", x=c(2.2,2.2), y=c(0.81, 0.85), label = "a", size = 3) + #add text to the graphic for posthoc letters
  annotate("text", x=c(0.85,0.85,0.85,2.25), y=c(0.96,1.04,1.2,0.95), label = "ab", size = 3) + #add text to the graphic for posthoc letters
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
All.Area <- do.call(data.frame,aggregate(Area ~ Day*Treatment*Secondary, data = seed.size, function(x) c(mean = mean(x), se = std.error(x)))) #calculate mean and sem
All.Area <- All.Area[with(All.Area, order(Day, Treatment)), ] #order by day and then treatment
write.table (All.Area, file="/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_oa/RAnalysis/Output/All.Shell.Area.csv", sep=",", row.names = FALSE) #save data to output file

##### CAPTURE STATISTICAL OUTPUT AND FIGURES TO FILE #####
setwd(file.path(mainDir, 'Output'))

capture.output(Flow1.Rate, FL1.res, Flow2.Rate, FL2.res, temp1.tank.res, temp1.trt.res, 
               TA1.tank.res, TA1.trt.res, temp2.tank.res, temp2.trt.res, 
               TA2.tank.res, TA2.trt.res, cells1.tank.res, cells1.trt.res,
               cells2.tank.res, cells2.trt.res,
               file="Geoduck_Juvenile_Descriptive_Statistics.txt")

capture.output(Init.res,E1.pairs, CG.res, CG.pairs, Rexp.res, E2.pairs, file="Geoduck_Juvenile_Statistical_Results.txt")

#Capture Figures to File
Figure1.Size <- arrangeGrob(Fig.Exp1.size,Fig.CG.size, Fig.Exp2.D10.size, Fig.Exp2.D23.size, ncol=2)
ggsave(file="Geoduck_Size.pdf", Figure1.Size, width = 6.5, height = 8, units = c("in"))

FigureS1 <- arrangeGrob(All.pH, All.temp, ncol=1)
ggsave(file="Fig.S1.pdf", FigureS1, width = 6, height = 8, units = c("in"))

setwd(file.path(mainDir, 'Data'))
