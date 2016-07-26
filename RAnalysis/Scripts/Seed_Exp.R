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
library("seacarb") #seawater carbonate chemistry
library("reshape") #reshape data
library("plotrix") #functions in tapply
library("ggplot2")

#Required Data files
#pH_Calibration_Files/
#SW_Chem_Juvenile_Geoduck.csv
#Cell_Counts_Juvenile_Geoduck.csv
#Size_Juvenile_Geoduck.csv
#Avetch_Temperature_Juvenile_Geoduck.csv

#############################################################
setwd("/Users/hputnam/MyProjects/Geoduck_Epi/project_juvenile_geoduck_OA/RAnalysis/Data/") #set working directory


#####SEAWATER CHEMISTRY ANALYSIS FOR DISCRETE MEASUREMENT#####

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

##### Seacarb Calculations #####

#Calculate CO2 parameters using seacarb
carb.output <- carb(flag=8, var1=SW.chem$pH.Total, var2=SW.chem$TA/1000000, S= SW.chem$Salinity, T=SW.chem$Temperature, P=0, Pt=0, Sit=0, pHscale="T", kf="pf", k1k2="l", ks="d") #calculate seawater chemistry parameters using seacarb

carb.output$ALK <- carb.output$ALK*1000000 #convert to µmol kg-1
carb.output$CO2 <- carb.output$CO2*1000000 #convert to µmol kg-1
carb.output$HCO3 <- carb.output$HCO3*1000000 #convert to µmol kg-1
carb.output$CO3 <- carb.output$CO3*1000000 #convert to µmol kg-1
carb.output$DIC <- carb.output$DIC*1000000 #convert to µmol kg-1

carb.output <- cbind(SW.chem$Measure.Date,  SW.chem$Tank,  SW.chem$Treatment, carb.output) #combine the sample information with the seacarb output
colnames(carb.output) <- c("Date",  "Tank",  "Treatment",	"flag",	"Salinity",	"Temperature",	"Pressure",	"pH",	"CO2",	"pCO2",	"fCO2",	"HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat", 	"Calcite.Sat") #Rename columns to describe contents
carb.output <- subset(carb.output, select= c("Date",  "Tank",  "Treatment",	"Salinity",	"Temperature",		"pH",	"CO2",	"pCO2",	"HCO3",	"CO3",	"DIC", "TA",	"Aragonite.Sat"))

##### Seawater Chemistry Descriptive Statistics #####
#Tanks
mean_pCO2=tapply(carb.output$pCO2, carb.output$Tank, mean)
se_pCO2=tapply(carb.output$pCO2, carb.output$Tank, std.error)
mean_Temp=tapply(carb.output$Temperature, carb.output$Tank, mean)
se_Temp=tapply(carb.output$Temperature, carb.output$Tank, std.error)
mean_Sal=tapply(carb.output$Salinity, carb.output$Tank, mean)
se_Sal=tapply(carb.output$Salinity, carb.output$Tank, std.error)
mean_TA=tapply(carb.output$TA, carb.output$Tank, mean)
se_TA=tapply(carb.output$TA, carb.output$Tank, std.error)
mean_pH=tapply(carb.output$pH, carb.output$Tank, mean)
se_pH=tapply(carb.output$pH, carb.output$Tank, std.error)
mean_DIC=tapply(carb.output$DIC, carb.output$Tank, mean)
se_DIC=tapply(carb.output$DIC, carb.output$Tank, std.error)

#Treatments
gmean_pCO2 <- tapply(carb.output$pCO2, carb.output$Treatment, mean)
gse_pCO2 <- tapply(carb.output$pCO2, carb.output$Treatment, std.error)
gmean_Temp <- tapply(carb.output$Temperature, carb.output$Treatment, mean)
gse_Temp <- tapply(carb.output$Temperature, carb.output$Treatment, std.error)
gmean_Sal <- tapply(carb.output$Salinity, carb.output$Treatment, mean)
gse_Sal <- tapply(carb.output$Salinity, carb.output$Treatment, std.error)
gmean_TA <- tapply(carb.output$TA, carb.output$Treatment, mean)
gse_TA <- tapply(carb.output$TA, carb.output$Treatment, std.error)
gmean_pH <- tapply(carb.output$pH, carb.output$Treatment, mean)
gse_pH <- tapply(carb.output$pH, carb.output$Treatment, std.error)
gmean_DIC <- tapply(carb.output$DIC, carb.output$Treatment, mean)
gse_DIC <- tapply(carb.output$DIC, carb.output$Treatment, std.error)

mean.carb.output <- rbind(gmean_pCO2, gse_pCO2, gmean_pH, gse_pH, gmean_Temp, gse_Temp, gmean_Sal, gse_Sal, gmean_TA, gse_TA, gmean_DIC, gse_DIC)
mean.carb.output <- as.data.frame(mean.carb.output)
row.names(mean.carb.output) <- c("mean pCO2", "SE pCO2", "mean pH", "SE pH", "mean Temperature", "SE Temperature", "mean Salinity", "SE Salinity", "mean Total Alkalinity", "SE Total Alkalinity", "mean DIC", "SE DIC")
mean.carb.output$Variables <- row.names(mean.carb.output)
write.table (mean.carb.output, file="/Users/hputnam/MyProjects/Geoduck_Epi/project-geoduck-oa/RAnalysis/Output/Seawater_chemistry_table_Output_Seed.csv", sep=",", row.names = FALSE)

##### Algal Feed Density Cell Counts #####

cell.counts <- read.csv("Cell_Counts_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
cell.counts$Avg.Cells <- rowMeans(cell.counts[,c("Count1",  "Count2")], na.rm = TRUE) #calculate average of counts
cell.counts$cells.ml <- cell.counts$Avg.Cells/cell.counts$Volume.Counted #calculate density

#Tanks
mean_cells=tapply(cell.counts$cells.ml, cell.counts$Tank, mean, na.rm = TRUE)
se_cells=tapply(cell.counts$cells.ml, cell.counts$Tank, std.error, na.rm = TRUE)

#Treatments
gmean_cells <- tapply(cell.counts$cells.ml, cell.counts$Treatment, mean, na.rm = TRUE)
gse_cells <- tapply(cell.counts$cells.ml, cell.counts$Treatment, std.error, na.rm = TRUE)

##### Seed Size #####

seed.size <- read.csv("Size_Juvenile_Geoduck.csv", header=TRUE, sep=",", na.strings="NA") #load data with a header, separated by commas, with NA as NA
Day.10.size <-subset(seed.size, Day=="Day10")
Day.23.size <-subset(seed.size, Day=="Day23")
avg.area <- aggregate(Area ~ Day*Treatment, data=seed.size, mean)
se.area <- aggregate(Area ~ Day*Treatment, data=seed.size, std.error)
Shell.size <- cbind(avg.area, se.area$Area)
colnames(Shell.size) <- c("Day", "Treatment", "avg.area", "se.area")
Day.10 <-subset(Shell.size, Day=="Day10")
Day.23 <-subset(Shell.size, Day=="Day23")
colnames(Day.10) <- c("Day", "Treatment", "avg.area", "se.area")
colnames(Day.23) <- c("Day", "Treatment", "avg.area", "se.area")

plot(c(1,3), c(20,40), type="n", xaxt = "n", ylab=expression(paste("Shell Area mm2")), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=Day.10$avg.area,uiw=Day.10$se.area, liw=Day.10$se.area,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

D.10 <- lm(Area ~ Treatment, data=Day.10.size)
anova(D.10)

plot(c(1,3), c(20,40), type="n", xaxt = "n", ylab=expression(paste("Shell Area mm2")), xlab=expression(paste("Treatment")))
axis(1, at=1:3, labels=c("pH 7.89", "pH 7.38", "pH 7.04"))
plotCI(x=c(1:3), y=Day.23$avg.area,uiw=Day.23$se.area, liw=Day.23$se.area,add=TRUE,gap=0.001, pch=20, col=c("blue", "pink", "red"))

D.23<- lm(Area ~ Treatment, data=Day.23.size)
anova(D.23)

X <- aov(Area ~ Treatment*Day, data=seed.size)
anova(X)
posthoc <- TukeyHSD(x=X, 'Treatment', conf.level=0.95)
posthoc
hist(X$residuals)

color= c("blue","pink","red")
interaction.plot(seed.size$Day, seed.size$Treatment, seed.size$Area, fun=mean, col = color, lty=1, legend=FALSE, xaxt = "n", frame.plot=TRUE, ylab=expression(paste("Seed Shell Area mm"^"2")))
axis(1, at=c(1,2,3,4), labels=c("Day0", "Day1", "Day10", "Day23", "Day51"))
legend("topleft", c("pH 7.89", "pH 7.38", "pH 7.04"), col=c("blue", "pink", "red"), lty=1, bty="n", cex=0.8) 

rect <- data.frame(xmin=3.1, xmax=Inf, ymin=-Inf, ymax=Inf)

FigX <- ggplot(Shell.size, aes(x=Day, y=avg.area, group=Treatment)) + 
  geom_errorbar(aes(ymin=Shell.size$avg.area-Shell.size$se.area, ymax=Shell.size$avg.area+Shell.size$se.area), colour="black", width=.1, position = position_dodge(width = 0.2)) +
  geom_line(aes(linetype=Treatment), size = 0.5, position = position_dodge(width = 0.2)) +   
  geom_point(aes(shape=Treatment), size = 4, position = position_dodge(width = 0.2)) +
  xlab("TimePoint") +
  ylab("Seed Shell Area mm^2") +
  ylim(20,90) +
  geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            color="gray",
            alpha=0.5,
            inherit.aes = FALSE) +
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank()) #Set the plot background
FigX


geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              color="grey20",
              alpha=0.5,
              inherit.aes = FALSE)


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

#Seed size
N<-aggregate(N~Treatment*Day,data=seed.size, sum) 
boxplot(Area~Treatment*Day,data=seed.size,col=c("blue", "pink", "red"), xaxt = "n", frame.plot=TRUE, ylab=expression(paste("Seed Shell Area mm"^"2")))
axis(1, at=c(5,8,11), labels=c("Day1", "Day10", "Day23"))
legend("topleft", c("pH 7.89", "pH 7.38", "pH 7.04"), fill=c("blue", "pink", "red"), bty="n", cex=0.6) 

# Add data points
# mylevels<-levels(seed.size$Treatment)
# for(i in 1:length(mylevels))
# {
#   thislevel<-mylevels[i]
#   thisvalues<-seed.size[seed.size$Treatment==thislevel, "Avg.Area"]
#   
#   # take the x-axis indices and add a jitter, proportional to the N in each level
#   myjitter<-jitter(rep(i, length(thisvalues)), amount=levelProportions[i]/10)
#   points(myjitter, thisvalues, pch=20, col=rgb(0,0,0,.2))   
# }

# 
title("Treatment Conditions Juvenile Experiment", outer=TRUE)
dev.off()

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



