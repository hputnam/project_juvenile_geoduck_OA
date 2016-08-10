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
#pH_Calibration_Files/
#SW_Chem_Juvenile_Geoduck.csv
#Cell_Counts_Juvenile_Geoduck.csv
#Size_Juvenile_Geoduck.csv
#Avetch_Temperature_Juvenile_Geoduck.csv

#############################################################
setwd("/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/Data/") #set working directory
mainDir<-'/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/' #set main directory

##### CONTINUOUS HEADER TANK PH EXPOSURE 1 #####

##### CONTINUOUS HEADER TANK PH INDOOR COMMON GARDEN #####

##### CONTINUOUS HEADER TANK PH EXPOSURE 2 #####

##### CONTINUOUS TANK TEMPERATURE EXPOSURE 1 #####

##### CONTINUOUS TANK TEMPERATURE INDOOR COMMON GARDEN #####

##### CONTINUOUS TANK TEMPERATURE OUTDOOR COMMON GARDEN #####

##### CONTINUOUS WATERBATH TEMPERATURE EXPOSURE 2 #####

##### SEAWATER CHEMISTRY ANALYSIS FOR DISCRETE MEASUREMENT#####

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

##### SEACARB CALCULATIONS Seacarb Calculations #####

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
#Calcualte descriptive stats by Tank
ccarb <- melt(carb.output[,c(1,3,4,5:14)], id.vars=c("Exposure", "Treatment", "Tank")) #reshape data into long format

SWC.Tanks <- ddply(ccarb, c("Exposure", "Tank", "variable"), summarise,
      N = length(na.omit(value)), #count the sample size removing NA
      mean = mean(value), #calculate average 
      sem = sd(value)/sqrt(N)) #calcualte the standard error of the mean

#Calcualte descriptive stats by Treatment
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

#Plotting
par(mar = rep(2, 4))
par(mfrow=c(5,4))
Exposure1.Plot<-function(index) {plot(Exposure1.long[,index] ~ Exposure1.long$Treatment, main=names(Exposure1.long[index]), pch=16,xlab="Treatment",ylab=names(Exposure1.long)[index])}
lapply(2:21,FUN=Exposure1.Plot)
div.off()

par(mar = rep(2, 4))
par(mfrow=c(5,4))
Exposure2.Plot<-function(index) {plot(Exposure2.long[,index] ~ Exposure2.long$Treatment, main=names(Exposure2.long[index]), pch=16,xlab="Treatment",ylab=names(Exposure2.long)[index])}
lapply(2:21,FUN=Exposure2.Plot)

##### ALGAL FEED DENSITY COUNTS EXPOSURE 1 #####

cell.counts <- read.csv("Cell_Counts_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
cell.counts$Avg.Cells <- rowMeans(cell.counts[,c("Count1",  "Count2")], na.rm = TRUE) #calculate average of counts
cell.counts$cells.ml <- cell.counts$Avg.Cells/cell.counts$Volume.Counted #calculate density

#Tanks
mean_cells=tapply(cell.counts$cells.ml, cell.counts$Tank, mean, na.rm = TRUE)
se_cells=tapply(cell.counts$cells.ml, cell.counts$Tank, std.error, na.rm = TRUE)

#Treatments
gmean_cells <- tapply(cell.counts$cells.ml, cell.counts$Treatment, mean, na.rm = TRUE)
gse_cells <- tapply(cell.counts$cells.ml, cell.counts$Treatment, std.error, na.rm = TRUE)

##### ALGAL FEED DENSITY COUNTS EXPOSURE 2 #####

##### Plot Tank and Treatment mean ± se #####
pdf("/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/Output/running_carbonate_chemistry_tanks_Seed.pdf")
par(cex.axis=0.8, cex.lab=0.8, mar=c(5, 5, 4, 2),mgp=c(3.7, 0.8, 0),las=1, mfrow=c(3,3), oma=c(0,0,2,0))

#Tanks
plot(c(11,16),c(0,6000),type="n",ylab=expression(paste("pCO"["2"])), xlab=expression(paste("Tank")))
plotCI(x=c(11,12,13,14,15,16), y=mean_pCO2,uiw=se_pCO2, liw=se_pCO2,add=TRUE,gap=0.001)

plot(c(11,16),c(6.5,8.5),type="n",ylab=expression(paste("pH")), xlab=expression(paste("Tank")))
plotCI(x=c(11,12,13,14,15,16), y=mean_pH,uiw=se_pH, liw=se_pH,add=TRUE,gap=0.001)

plot(c(11,16),c(12,15),type="n",ylab=expression(paste("Temperature °C")), xlab=expression(paste("Tank")))
plotCI(x=c(11,12,13,14,15,16), y=mean_Temp,uiw=se_Temp, liw=se_Temp,add=TRUE,gap=0.001)

plot(c(11,16),c(25,29),type="n",ylab=expression(paste("Salinity")), xlab=expression(paste("Tank")))
plotCI(x=c(11,12,13,14,15,16), y=mean_Sal,uiw=se_Sal, liw=se_Sal,add=TRUE,gap=0.001)

plot(c(11,16),c(1800,2200),type="n",ylab=expression(paste("Total Alkalinity µmol kg"^"-1")), xlab=expression(paste("Tank")))
plotCI(x=c(11,12,13,14,15,16), y=mean_TA,uiw=se_TA, liw=se_TA,add=TRUE,gap=0.001)

plot(c(11,16),c(1850,2400),type="n",ylab=expression(paste("DIC µmol kg"^"-1")), xlab=expression(paste("Tank")))
plotCI(x=c(11,12,13,14,15,16), y=mean_DIC,uiw=se_DIC, liw=se_DIC,add=TRUE,gap=0.001)

plot(c(11,16),c(20000,60000),type="n", ylab=expression(paste("Algal Feed (Cells ml"^"-1",")")), xlab=expression(paste("Tank")))
plotCI(x=c(11,12,13,14,15,16), y=mean_cells,uiw=se_cells, liw=se_cells,add=TRUE,gap=0.001)

title("Tank Conditions Juvenile Experiment", outer=TRUE)
dev.off()


pdf("/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/Output/running_carbonate_chemistry_treatments_Seed.pdf")
par(cex.axis=0.8, cex.lab=0.8, mar=c(5, 5, 4, 2),mgp=c(3.7, 0.8, 0),las=1, mfrow=c(3,3), oma=c(0,0,2,0))

#Treatments
plot(c(1,2,3),c(0,6000,6000), xaxt = "n", type="n",ylab=expression(paste("pCO"["2"])), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=gmean_pCO2,uiw=gse_pCO2, liw=gse_pCO2,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

plot(c(1,2,3),c(6.8,8.5,8.5), xaxt = "n", type="n",ylab=expression(paste("pH")), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=gmean_pH,uiw=gse_pH, liw=gse_pH,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

plot(c(1,2,3),c(13,15,15), xaxt = "n", type="n",ylab=expression(paste("Temperature °C")), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=gmean_Temp,uiw=gse_Temp, liw=gse_Temp,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

plot(c(1,2,3),c(25,29,29), xaxt = "n", type="n",ylab=expression(paste("Salinity")), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=gmean_Sal,uiw=gse_Sal, liw=gse_Sal,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

plot(c(1,2,3),c(1800,2200,2200), xaxt = "n", type="n",ylab=expression(paste("Total Alkalinity µmol kg"^"-1")), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=gmean_TA,uiw=gse_TA, liw=gse_TA,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

plot(c(1,2,3),c(1800,2400,2400), xaxt = "n", type="n",ylab=expression(paste("DIC µmol kg"^"-1")), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=gmean_DIC,uiw=gse_DIC, liw=gse_DIC,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

plot(c(1,2,3),c(20000,90000,90000), xaxt = "n", type="n", ylab=expression(paste("Algal Feed (Cells ml"^"-1",")")), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=gmean_cells,uiw=gse_cells, liw=gse_cells,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

##### JUVENILE GEODUCK SHELL SIZE #####

seed.size <- read.csv("Size_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
Initial <- subset(seed.size, Timepoint=="Initial", select = Date:Area)
Init.avg.area <- aggregate(Area ~ Day*Treatment, data=Initial, mean)
Init.se.area <- aggregate(Area ~ Day*Treatment, data=Initial, std.error)
Init.Shell.size <- cbind(Init.avg.area, Init.se.area$Area)
colnames(Init.Shell.size) <- c("Day", "Treatment", "avg.area", "se.area")

ReExp <- subset(seed.size, Timepoint=="Reexposure", select = Date:Area)
ReExp.avg.area <- aggregate(Area ~ Day*Treatment*Secondary, data=ReExp, mean)
ReExp.se.area <- aggregate(Area ~ Day*Treatment*Secondary, data=ReExp, std.error)
ReExp.Shell.size <- cbind(ReExp.avg.area, ReExp.se.area$Area)
colnames(ReExp.Shell.size) <- c("Day", "Treatment", "Secondary", "avg.area", "se.area")
ReExp.Shell.size$Groups <- c("Amb-Amb", "Med-Amb", "High-Amb", "Amb-Med", "Med-Med", "High-Med")

# Significance testing for initial exposure and common garden
Init <- aov(log10(Area) ~ Treatment*Day, data=Initial)
anova(Init)
par(mfrow=c(3,3))
hist(Init$residuals)
boxplot(Init$residuals)
plot(Init)

#posthoc <- TukeyHSD(x=X,conf.level=0.95)
#posthoc

Day1 <- subset(Initial, Day==1, select = Date:Area)
Day10 <- subset(Initial, Day==10, select = Date:Area)
Day23 <- subset(Initial, Day==23, select = Date:Area)
Day51 <- subset(Initial, Day==51, select = Date:Area)
Day135 <- subset(Initial, Day==135, select = Date:Area)

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

rect <- data.frame(xmin=3.1, xmax=Inf, ymin=-Inf, ymax=Inf) #identify background shading placement

Fig1 <- ggplot(Init.Shell.size, aes(x=Day, y=avg.area, group=Treatment)) + 
  geom_errorbar(aes(ymin=Init.Shell.size$avg.area-Init.Shell.size$se.area, ymax=Init.Shell.size$avg.area+Init.Shell.size$se.area), colour="black", width=.1, position = position_dodge(width = 0.6)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.6)) +   
  geom_point(aes(shape=Treatment), size = 4, position = position_dodge(width = 0.6)) +
  xlab("Days") +
  ylab("Seed Shell Area mm^2") +
  ylim(20,130) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank()) + #Set the plot background
  geom_rect(data=rect, aes(xmin=25, xmax=137, ymin=ymin, ymax=ymax),
            color="gray",
            alpha=0.5,
            inherit.aes = FALSE)
Fig1

# Significance testing for secondary exposure
Rexp <- aov(log10(Area) ~ Treatment*Secondary, data=ReExp)
anova(Rexp)
#Check Assumptions
par(mfrow=c(3,3))
hist(Rexp$residuals)
boxplot(Rexp$residuals)
plot(Rexp)

Rexp.posthoc <- TukeyHSD(x=Rexp,conf.level=0.95)
Rexp.posthoc

#Day145 <- subset(seed.size, Day==145, select = Date:Area)

Fig2 <- ggplot(ReExp.Shell.size, aes(x=Secondary, y=avg.area, group=Treatment)) + 
  geom_errorbar(aes(ymin=ReExp.Shell.size$avg.area-ReExp.Shell.size$se.area, ymax=ReExp.Shell.size$avg.area+ReExp.Shell.size$se.area), colour="black", width=.1, position = position_dodge(width = 0.05)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.05)) +   
  geom_point(aes(shape=Treatment), size = 4, position = position_dodge(width = 0.05)) +
  xlab("Secondary Exposure") +
  ylab("Seed Shell Area mm^2") +
  ylim(20,130) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank()) #Set the plot background

Fig2

#CAPTURE ALL STATISTICAL OUTPUT TO FILE
setwd(file.path(mainDir, 'Output'))
capture.output(Init, Rexp, file="Geoduck_Juvenile_Statistical_Results.txt")

#CAPTURE ALL FIGURES TO FILE
FigureX <- arrangeGrob(Fig1, Fig2, ncol=1)
ggsave(file="Geoduck_Size.pdf", FigureX, width = 6.5, height = 8, units = c("in"))

setwd(file.path(mainDir, 'Data'))

##### Avtech Data #####
#Load WiSH data
wish.data <- read.csv("Wish_data_Seed.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
date.time <- sub("-",",", wish.data$Date...Time)
date.time <- strsplit(date.time, ",")
date.time <- data.frame(matrix(unlist(date.time), nrow=length(date.time), byrow=T),stringsAsFactors=FALSE)
temp.data <- wish.data[,grepl("Tank", colnames(wish.data))] #search for and subset columns containing the header name "Tank"
temp.data <- cbind(date.time, temp.data)
colnames(temp.data) <- c("Date", "Time", "Tank3", "Tank6", "Tank4", "Tank1", "Tank5", "Tank2")
pH.data <-cbind(date.time, wish.data$pH.Exp.Treat...Custom.Value, wish.data$ph.Exp.Control...Custom.Value, wish.data$Header.2...c...Custom.Value, wish.data$Header.1...tr...Custom.Value)
colnames(pH.data) <- c("Date", "Time", "TpH7.38", "TpH7.04", "HpH7.04","HpH7.38")

  
##plot temp data
plot(temp.data$Tank1,type="l", col="pink", ylab=expression(paste("Temperature °C")), xlab=expression(paste("Time")), ylim=c(10, 20))
lines(temp.data$Tank2, col="lightblue" )
lines(temp.data$Tank3, col="blue")
lines(temp.data$Tank4, col="red")
lines(temp.data$Tank5, col="darkred")
lines(temp.data$Tank6, col="darkblue")
legend("topleft", c("Tank1","Tank2", "Tank3","Tank4","Tank5", "Tank6" ), col=c("pink","lightblue", "blue", "red", "darkred", "darkblue"), bty="n", lwd=1, cex=0.6) 

#plot pH Data
plot(pH.data$TpH7.38,type="l", col="lightblue", ylab=expression(paste("pH")), xlab=expression(paste("Time")), ylim=c(6.8, 8.0))
lines(pH.data$TpH7.04, col="pink" )
lines(pH.data$HpH7.38, col="blue" )
lines(pH.data$HpH7.04, col="red" )
legend("topleft", c("Tank Low","Tank Super Low", "Header Low","Header Super Low"), col=c("lightblue","pink", "blue", "red"), bty="n", lwd=1, cex=0.6) 

##### CAPTURE STATISTICAL OUTPUT AND FIGURES TO FILE #####
setwd(file.path(mainDir, 'Output'))
capture.output(Init, Rexp, file="Geoduck_Juvenile_Statistical_Results.txt")

#CAPTURE ALL FIGURES TO FILE
FigureX <- arrangeGrob(Fig1, Fig2, ncol=1)
ggsave(file="Geoduck_Size.pdf", FigureX, width = 6.5, height = 8, units = c("in"))


