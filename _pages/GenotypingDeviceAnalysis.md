---
title: "Alison Thompson - GenotypingDeviceAnalysis"
layout: textlay
excerpt: "GenotypingDeviceAnalysis"
sitemap: false
permalink: /GenotypingDeviceAnalysis/
---

Introduction
------------

This code was used for a project to examine the genetics of single cells. The experients use a microfluidic device to partition single cells into small wells where a reaction can be carried out. This is an picture of the device:

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/IMG_photo_of_SDChip_at_FH_crop.jpg"|absolute_url}}) **Figure 1. Photo of the device with three replicate arrays.**

This device has three replicate arrays, each will be loaded with a solution containing cells and the reagents for the reaction. Each replicate array has 1024 wells, and each well will hold about 5 nanoliters of solution. The device shown here is made on a 2"x3" microscope slide. As the device is loaded, the cells are distributed randomly, so we first image the arrays to locate wells with single cells.

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/CellImage1.png"|absolute_url}}) **Figure 2. Image of cells in wells.**

This image is aquired using a fluorescence microscope, and because the cells are stained with a fluorescent dye, we can see them as the tiny white spots in the image. The cells are about 15 micrometers in diameter, and the wells are 200 micrometers wide and 400 micrometers tall. 34 images are captured to cover all the device area.

After imaging, the reaction (PCR) takes about 2 hours. The reaction will only occur if a specific genetic sequence is present in the well. If the reaction does occur, we will observe an increased fluorescence signal in that well as a dye is released. We are looking for three different genetic sequences with three different dyes, so we image the device in three channels.

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Typhoon_composite.png"|absolute_url}}) **Figure 3. Image of fluorescence post-reaction.**

From these 37 images, we will have enough information to report single-cell genotypes in our three replicates! To do this, we need to know which wells contain single cells, and whether they are negative or positive for the control sequence, the wild-type sequence, and/or the mutant sequence. The cells can be either wild-type, mutant, or heterozygous.

Directory Setup and Data Cleaning
---------------------------------

The images are processed using ImageJ using a macro language to overlay a grid of ROI's on the wells and extract the critical data. These include the number of cell-sized particles per well and the integrated fluorescence value for each well. ImageJ exports these values as .csv's. We are going to merge the 37 csv's into a single data frame using R.

Load the necessary R packages

``` r
require(plyr); require(dplyr); require(tidyr); require(reshape2); require(ggplot2)
require(ggExtra); require(ggrepel); require(gridExtra);  require(RColorBrewer); require(knitr)
options(stringsAsFactors = FALSE)
```

Create a results folder and read in metadata about the experimental conditions

``` r
#Sets up the directory structure for where to look for data, and where to save results.  
#This assumes that the raw images have been processed using the ImageJ macros.
runDate <- format(Sys.time(), "%Y%m%d-%H%M")
output <- c(paste0("./OutputData-",runDate))
dir.create(output, showWarnings = FALSE)
metadataDir <- c("./Metadata")
fluorInput <- c("./TyphoonCSVs")
writeLines(capture.output(sessionInfo()),
           paste(output,"/TyphoonDR-sessionInfo.txt", sep = ""))

#Allows us to read in the probeScheme to determine which fluorophore is amp/wt/mut and also read in which template is in each array.
arrayList <- c("A1", "A2", "A3")
runMetadata <- read.delim("runinfo.txt")
probeScheme <- runMetadata[grepl( "^chip.*", runMetadata[,1])>0,]
#check these for new probe schemes, assumes ps1 is FAM-amp, HEX-WT, Cy5-MUT
ampFluor <- "FAM"; wtFluor <- "HEX"; mutFluor <- "Cy5"; pSch = "PS1"
if(grepl(".*ZygSwap", probeScheme)== T) {
    ampFluor <- "FAM"; wtFluor <- "Cy5"; mutFluor <- "HEX"; pSch = "ZygSwap"}
if(grepl(".*AmpSwap", probeScheme)== T) {
    ampFluor <- "Cy5"; wtFluor <- "HEX"; mutFluor <- "FAM"; pSch = "AmpSwap"}
fluorList <- c(ampFluor, wtFluor, mutFluor)
TemplateRows <- c(runMetadata[grepl( "control.*", runMetadata[,1])>0,], runMetadata[grepl( "array.*", runMetadata[,1])>0,])
TemplateList <- gsub("^.*\\s([A-z])", "\\1", TemplateRows)
```

**Experimental Conditions:**
The probe scheme used for this chip is **ZygSwap**.
The amplification probe is in the **FAM** channel, the wild-type probe is in the **Cy5** channel, and the mutant probe is in the **HEX** channel. The templates used were **OCI-AML3, OCI-AML3, OCI-AML3** in replicates A1, A2, and A3, respectively. The run outputs will be saved in the folder **./OutputData-20180324-1445**.

*note: in this example all replicates contain OCI-AML3, a cell line.*

Create a single dataframe from all Typhoon images. Label arrays, rows, and columns

``` r
#Read in the individual csv's for each fluorphore to make a single data frame.  
fileNames <- list.files(path = fluorInput, pattern = "*.csv", full.names = TRUE)
dataList <- lapply(fluorList, function(x)
    {read.csv(file = fileNames[grep(x, fileNames)], stringsAsFactors = FALSE)})
names(dataList) <- fluorList
#Give each line a "Well" number as an ID variable
dataList <- lapply(dataList, function(x) {colnames(x)[1] <- "Well";x})
allFluors <- melt(dataList, id=colnames(dataList[[1]]))
#The last column in the Fluorphore (dye).  Label it "Fluor"
colnames(allFluors)[length(colnames(allFluors))] <- "Fluor"
#ImageJ labeled the area is X.Area.  Change it to "Area"
names(allFluors)[names(allFluors) == "X.Area"] <- "Area"

rm(dataList)

#Name the replicates "A1 - A3".  This will allow us to compare their statistics.
allFluors$Array[allFluors$Well <= 1024] <- "A1"
allFluors$Array[allFluors$Well > 1024 & allFluors$Well <= 2048] <- "A2"
allFluors$Array[allFluors$Well > 2048 & allFluors$Well <= 3072] <- "A3"

#Add a template column, which tells us what sample was in the array.
ArraytoTemplate <- data_frame(Array = arrayList, Template = TemplateList)
allFluors <- left_join(allFluors, ArraytoTemplate, by = c("Array"))

#Generate an adusted integrated density (AID) equal to the mean pixel intensity in the well
allFluors$AID <- allFluors$IntDen/(allFluors$Area)
allFluors$Fluor <- factor(allFluors$Fluor, levels = fluorList)
allFluors$Array <- factor(allFluors$Array)

procData <- allFluors[c("Well", "Fluor", "Array", "AID")]
ampArea <- allFluors %>% filter(Fluor == ampFluor) %>% select(Well, Area, Array, Template)
```

    ## Warning: package 'bindrcpp' was built under R version 3.3.3

``` r
wideForm <- procData %>% group_by(Array) %>% spread(Fluor, AID)
colnames(wideForm)[which(grepl(ampFluor, colnames(wideForm)))] <- "ampFluor"
colnames(wideForm)[which(grepl(wtFluor, colnames(wideForm)))] <- "wtFluor"
colnames(wideForm)[which(grepl(mutFluor, colnames(wideForm)))] <- "mutFluor"

primeTime <- merge(ampArea, wideForm)
primeTime <- primeTime[order(primeTime$Well),]
primeTime$ArrayRow <- c(rep(1:16, each = 64, times = 3))
primeTime$ArrayCol <- seq(1:64)

rm(procData, ampArea, wideForm)
```

Read in cell images data

``` r
if(dir.exists("./OlympusCSVs")){

cellInput <- c("./OlympusCSVs")
#Find all the cell image csv's
filenamesOlympus <- list.files(path = cellInput, pattern = ".csv", full.names = T)
cellData <- lapply(filenamesOlympus, function(x) read.csv(x, header = T))
cellMelt <- melt(cellData, id = colnames(cellData[[1]]))
colnames(cellMelt) <- c("imageWell", "cellCount", "imageId")
cellMelt$ROIorder <- seq(1:3072)
#All 34 cell images have their own csv, and when we read them in the well order is different from the order in our main data frame.  We'll create a lookup table to translate the well ID's.  Then we can tack a cell count onto each well!  
lookup <- data.frame(ROIorder = 1:3072)
#most images have 11 columns x 8 rows of wells, but some only have 8 columns x 8 rows.  This lookup table only works if we image the cells in the same pattern every time!
lookup$ArrayCol <- rep(c(rep(1:11, len = 88), rep(12:22, len = 88),
                         rep(23:33, len = 88), rep(34:44, len = 88),
                         rep(45:55, len = 88), rep(56:64, len = 144),
                         rep(45:55, len = 88), rep(34:44, len = 88),
                         rep(23:33, len = 88), rep(12:22, len = 88),
                         rep(1:11, len = 88)), 3)
lookup$ArrayRow <- rep(c(rep(1:8, each = 11, len = 440),
                         rep(1:8, each = 9),rep(9:16, each = 9),
                         rep(9:16, each = 11, len = 440)),3)

lookup$Array <- rep(c("A1", "A2", "A3"), each = 1024)
#arrange the data so that the wells are in the same order as our main data frame
lookup <- arrange(lookup, Array, ArrayRow, ArrayCol)
#Add a Well column which will match the Well ID's in the main data frame
lookup$Well <- seq(from = 1, to = 3072)
fullMergeData <- join(lookup, cellMelt, by = "ROIorder", type = "left")
cellsReadyToRoll <- select(fullMergeData, Well, cellCount)
#Add a cell count to our data frame
primeTime <- left_join(primeTime, cellsReadyToRoll, by = "Well")
primeTime$cellCount <- as.numeric(primeTime$cellCount)
}
```

Plot the fluorescence intenisty in each fluorophore/dye channel. It looks pretty messy here, but we expect to see a bimodal distribution of negative and positive wells for each channel. We'll do some filtering and corrections so that we only analyze the wells that filled and imaged without error.

``` r
ggplot(data = allFluors, aes(x=Well, y=AID)) +
  geom_point(color = "skyblue") +
  facet_grid(Fluor ~ Array, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(panel.grid = element_blank())
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/scatter%20plot%20of%20AID%20by%20Well,%20where%20each%20panel%20is%20one%20channel-1.png"|absolute_url}})

**Figure 4. Scatter plot of the AIDs (fluorescence intensity) for each replicate array and fluorophore.**

``` r
ggplot(data = allFluors, aes(x=AID)) +
      geom_histogram(fill="skyblue",
          bins=30, color="grey50")+
  xlim(0, 120000) +
  xlab("AID") + ylab("Density") +
  stat_density(geom="line", color="red") +
  ggtitle("Custer Identification") +
  facet_grid(Fluor ~ Array) +
  theme_bw() +
  theme(panel.grid = element_blank())
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/histogram%20of%20AID,%20for%20each%20array%20and%20fluorophore-1.png"|absolute_url}})

**Figure 5. Histogram of the AIDs (fluorescence intensity) for each replicate array and fluorophore.**

Thresholding
------------

The first correction we perform is an area threshold. This way, we can eliminate wells that didn't fully fill from further analysis.

``` r
tArea <- 27
ggplot(data = primeTime %>% filter(Area > 5), aes(x=Area)) +
    geom_density(fill = "skyblue") + facet_grid(. ~ Array) +
    geom_vline(xintercept = tArea, color = "springgreen4") +
    theme_bw() +
    theme(panel.grid = element_blank())
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/area%20threshold-1.png"|absolute_url}})

``` r
#Create a Filled column that will be True if the well is filled
ruleThemAll <- primeTime %>% mutate(Filled = Area > tArea)
```

**Figure 6. Well Area Threshold. A density curve of well area for each array of this chip is shown with the area threshold drawn in green at 27. All wells below the threshold are empty or partially filled. These individual wells will be marked as "Filled = FALSE" and will not be evaluated in results. This threshold can be adjusted if necessary in the code chunk above.**

We also noticed some sample evaporation occured at the edges of each replicate. We can remove wells on the edges from further analysis.

``` r
ruleThemAll <- primeTime %>% group_by(Array) %>% mutate(Filled = ArrayCol < 63 & ArrayCol > 2 &
                                                          ArrayRow < 16 & ArrayRow > 1 &
                                                          Area > tArea)
```

Our dyes and filters are not perfectly separated in the visible spectrum, and some of our FAM fluorescence shows up in our HEX channel. The more FAM fluorescence in a well, the more artificat we see in the HEX channel for that well. We previously determined that a correction factor of 0.45 could account for this.

``` r
#Modify dfAllControls to use HEX correction, where m is avgHEXSlope, HEXFluor(adj) = HEXFluor(initial) - m*FAMFluor
avgHEXSlope = 0.45
if(pSch == "PS1") {ruleThemAll <- mutate(ruleThemAll, wtFluor = wtFluor - avgHEXSlope*ampFluor)}
if(pSch == "ZygSwap") {ruleThemAll <- mutate(ruleThemAll, mutFluor = mutFluor - avgHEXSlope*ampFluor)}
if(pSch == "AmpSwap") {ruleThemAll <- mutate(ruleThemAll, wtFluor = wtFluor - avgHEXSlope*mutFluor)}
```

After filtering and correcting, we can finally draw thresholds to determine whether a well is positive vs. negative for the genetic sequence!

I tried a few methods for finding the negative peak. It is *usually* the largest peak in a density curve of fluorescence intensity, but sometimes the second mode of the bimodal distribution is larger. I've found the kmeans clustering, when setting the number of clusters as 2, works well to consistently identify the two populations.

``` r
#this function assigns two clusters to the long form of the data frame.  Cluster 0 = low AID cluster and 1 = high AID cluster.  If centers are less that 10000 fluorescence units apart, data is considered a single cluster and all measurements are assigned to cluster 0.
AddCluster <- function(AID) {
  set.seed(112);
  km <- kmeans(AID, centers = 2, nstart = 25);
  low <- which.min(km[["centers"]]);
  sep <- max(km[["centers"]]) - min(km[["centers"]]);
  if (sep < 10000) {
    km$clust <- 0
  } else {
  km$clust[km$clust == low] <- 0;
  km$clust[km$clust >= 1] <- 1;
  }
  return(km$clust)
}

#call to function.  Creates a long data frame from a wide form and runs AddPeaks function to produce a "Cluster" column with the cluster assignment.
WithCluster <- ruleThemAll %>% gather(Fluor, AID, ampFluor:mutFluor) %>%
          filter(Filled == TRUE, !is.na(AID)) %>%
          group_by(Array, Fluor) %>%
          mutate(Cluster = AddCluster(AID))

#plots histograms colored by clusters for all runID's
ggplot(data = WithCluster, aes(x=AID)) +
      geom_histogram(aes(fill=as.factor(Cluster)),
          bins=30, color="grey50")+
  xlim(0, 120000) +
  xlab("AID") + ylab("Density") +
  stat_density(geom="line", color="red") +
  ggtitle("Custer Identification") +
  facet_grid(Fluor ~ Array) +
  theme_bw() +
  theme(panel.grid = element_blank())
```

    ## Warning: Removed 18 rows containing missing values (geom_bar).

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/calculate%20maximums%20via%20kmeans-1.png"|absolute_url}})

``` r
#function to find the full width half max standard deviation
FWHMSD <- function(X) {
  Y = density(X, na.rm = TRUE);
      Z = data.frame(Intensity = Y$x, Density = Y$y);
    Q = Z$Intensity[which.max(Z$Density)]
  S = filter(Z, Intensity > Q)
  W = filter(S, Density < (max(Z$Density)/2))
  return(2*(W$Intensity[1]-Q)/2.355)
}

#call to standard deviation function, make a table with SD's of AID's for each Fluor
snowFlake <- WithCluster %>% filter(Cluster == 0) %>%
          group_by(Array, Fluor) %>%
          do(data.frame(SD = FWHMSD(.$AID)))

#also calculate the center of cluster 0, while we test to see if we should use the center or the max of the density curve.  We can just calculate the mean of cluster 0, which would be equivalent to the center, since this is 1D data.
snowCenter <- WithCluster %>% filter(Cluster == 0) %>%
  group_by(Array, Fluor) %>%
  summarise(Center0 = mean(AID))

#also add the centers to the dataframe while we test using the center vs. the max.
snowPeaks <- left_join(snowFlake, snowCenter, by = c("Array", "Fluor"))
```

**Figure 7. Histogram of fluorescence AID for each channel and array, with color by cluster distinction**

``` r
#Max SD's x multiplier for each probe from negative templates (NTC or off-target plasmid).  
#See "Z:\grp\SDChipImageAnalysis\Multi-FWHMAnalysis-20170919-1439\20170919-1439SDstatsNoFilter_totals.xlsx" for calculations.
#change multiplier (xsigamp/mut/wt) to the x sigma you want to use for thresholding.  maxAID + xsig*sigma.
xsigamp <- 3
xsigmut <- 2
xsigwt <- 6

#Add thresholds to snowPeaksAdj
snowDrift <- snowPeaks %>%
  mutate(Threshold = if_else(Fluor == "ampFluor", Center0 + SD*xsigamp,
                     if_else(Fluor == "mutFluor", Center0 + SD*xsigmut,
                     if_else(Fluor == "wtFluor", Center0 + SD*xsigwt, 0))))

#Add this threshold to the df with all the AID's
snowDrift1 <- snowDrift %>% select(Array, Fluor, Threshold) %>% spread(Fluor, Threshold) %>% select(Array, ampThresh = ampFluor, mutThresh = mutFluor, wtThresh = wtFluor)

snowCone1 <- left_join(ruleThemAll, snowDrift1, by = c("Array"))
```

``` r
ampThreshPlot2 <- ggplot(snowCone1 %>% filter(Filled == TRUE), aes(x=ampFluor)) +
                            geom_histogram(aes(x = ampFluor), bins = 30, fill = "skyblue2") +
                            geom_vline(aes(xintercept = ampThresh), color = "orchid4") +
                            xlim(0, 50000) +
                            xlab("AID") + ylab("Density") +
                            facet_grid(. ~ Array) +
                            ggtitle("Amplification Threshold") +
                            theme_bw() +
                            theme(strip.text.y = element_blank(), axis.text.x = element_text(angle=90, hjust=1), panel.grid = element_blank())
ampThreshPlot2
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/view%20amp%20thresh%20without%20density%20curve-1.png"|absolute_url}}) **Figure 8. Threshold for amplification fluorophore.**

``` r
mut_vs_wt <- ggplot(snowCone1 %>% filter(Filled == TRUE),
            aes(x=wtFluor, y=mutFluor)) +
            geom_point(alpha = 0.6, stroke = 0, aes(color = Array)) +
            xlim(0,100000) +
            ylim(0,100000) +
            coord_fixed() +
            labs(x = "wtFluor", y = "mutFluor") +
            facet_wrap("Array", nrow = 2) +
            ggtitle("Wild-Type and Mutant Thresholds") +
            geom_vline(aes(xintercept = wtThresh), color = "orchid4") +
            geom_hline(aes(yintercept = mutThresh), color = "orchid4") +
            theme_bw() +
            theme(strip.text.y = element_blank(), axis.text.x = element_text(angle=90, hjust=1), panel.grid = element_blank())
mut_vs_wt
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Preliminary%20View%20of%20wt%20and%20mut%20Thresholds-1.png"|absolute_url}})

**Figure 9. Thresholds for mutant and wild-type fluorophore. A scatter plot of mutant vs. wild-type AIDs for each filled well is shown with the generated mutant and wild-type thresholds. These quadrants will be used to determine zygosity in each well.**

Results
-------

Using the calculated thresholds, we can assign a zygosity to each well where a reaction occured: wild-type, mutant, or heterozygous

``` r
#Define Amp positivity
snowCone1 <- mutate(snowCone1, AmpPos = FALSE)
snowCone1 <- snowCone1 %>% mutate(AmpPos = ampFluor >= ampThresh)
snowCone1$AmpPos <- as.factor(snowCone1$AmpPos)

#this turns any NA's into 0's
snowCone1$wtFluor[is.na(snowCone1$wtFluor)] <- 0
snowCone1$mutFluor[is.na(snowCone1$mutFluor)] <- 0
snowCone1$ampFluor[is.na(snowCone1$ampFluor)] <- 0

#Define zygosity
snowCone1 <- snowCone1 %>%
  mutate(Zygosity = if_else(Filled == TRUE & AmpPos == TRUE & wtFluor >= wtThresh & mutFluor < mutThresh, "WT",
                    if_else(Filled == TRUE & AmpPos == TRUE & wtFluor < wtThresh & mutFluor >= mutThresh, "MUT",
                    if_else(Filled == TRUE & AmpPos == TRUE & wtFluor >= wtThresh & mutFluor >= mutThresh, "HET", "UNCALLED"))))

#Define FluorQC.  This labels the filled wells that did not fit any of the zygosity definitions.
snowCone1 <- snowCone1 %>%
  mutate(FluorQC =
           if_else(Filled == TRUE & AmpPos == FALSE & wtFluor >= wtThresh & mutFluor < mutThresh, "wtFluor",
          if_else(Filled == TRUE & AmpPos == FALSE & wtFluor < wtThresh & mutFluor >= mutThresh, "mutFluor",
          if_else(Filled == TRUE & AmpPos == FALSE & wtFluor >= wtThresh & mutFluor >= mutThresh, "Both",
          if_else(Filled == TRUE & AmpPos == TRUE & wtFluor < wtThresh & mutFluor < mutThresh, "AmpOnly",
          if_else(Filled == TRUE & AmpPos == FALSE & wtFluor < wtThresh & mutFluor < mutThresh, "Negative", "ThisIsFine"))))))

snowCone1$Zygosity <- factor(snowCone1$Zygosity, levels = c("WT", "HET", "MUT","UNCALLED"))
snowCone1$FluorQC <- factor(snowCone1$FluorQC, levels = c("wtFluor", "mutFluor", "Both", "AmpOnly", "Negative", "ThisIsFine"))
```

``` r
FinalZygosity <- ggplot(snowCone1 %>% filter(Filled ==TRUE, AmpPos == TRUE | FluorQC == "Negative"),
            aes(x=wtFluor, y=mutFluor)) +
            geom_point(aes(color = Zygosity),alpha = 0.4, stroke = 0) +
            geom_vline(aes(xintercept = wtThresh), color = "springgreen4") +
            geom_hline(aes(yintercept = mutThresh), color = "springgreen4") +
        scale_color_manual(values = c("dodgerblue4","magenta4", "red3", "black"), limits = levels(snowCone1$Zygosity)) +
            xlim(0,100000) +
            ylim(0,100000) +
            coord_fixed() +
            labs(x = "wtFluor", y = "mutFluor") +
            facet_wrap("Array", nrow = 2) +
            ggtitle("Zygosity") +
            theme_bw() +
            theme(axis.text.x = element_text(angle=90, hjust=1), strip.text.y = element_text(size = 7, angle = 0))
FinalZygosity
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/View%20Results%20of%20wt%20and%20mut%20Thresholding-1.png"|absolute_url}})

**Figure 10. Zygosity of each well that passed quality control filtering. Color indicates zygosity.**

We can also take a peak at wells that did not get assigned a zygosity.

``` r
FluorQC <- ggplot(snowCone1 %>% filter(Filled ==TRUE, FluorQC != "Negative", FluorQC != "ThisIsFine"),
            aes(x=wtFluor, y=mutFluor)) +
            geom_point(aes(color = FluorQC),alpha = 0.6, stroke = 0) +
            geom_vline(aes(xintercept = wtThresh), color = "springgreen4") +
            geom_hline(aes(yintercept = mutThresh), color = "springgreen4") +
        scale_color_manual(values = c("dodgerblue4","magenta4", "red3", "green4", "orange4", "black"), limits = levels(snowCone1$FluorQC)) +
            xlim(0,100000) +
            ylim(0,100000) +
            coord_fixed() +
            labs(x = "wtFluor", y = "mutFluor") +
            facet_wrap("Array", nrow = 2) +
            ggtitle("Fluorescence QC") +
            theme_bw() +
            theme(axis.text.x = element_text(angle=90, hjust=1), strip.text.y = element_text(size = 7, angle = 0))
FluorQC
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/View%20Results%20Fluor%20QC-1.png"|absolute_url}})

**Figure 11. Condition of any wells that was filled but not assigned a zygosity. Since we required wells to have both amplification control positivity and one other positive reaction to be labeled with a zygosity, this removed some questionable wells.**

**Table showing well counts per array. Every filled well is either Amp pos or neg, FilledAmped = WT + MUT + HET + ampOnly; and FilledNoAmp = NOAMPwt + NOAMPmut + NOAMPboth + TripleNeg.**

``` r
snowCone1$AmpPos <- as.logical(snowCone1$AmpPos)

statsFilledWells <- snowCone1 %>% group_by(Array, Template) %>%
    summarise(TotalFilled = sum(Filled),
              FilledAmped = sum(Filled & AmpPos),
              WT = n_distinct(which(Zygosity == "WT")),
              MUT = n_distinct(which(Zygosity == "MUT")),
              HET = n_distinct(which(Zygosity == "HET")),
              ampOnly = n_distinct(which(FluorQC == "AmpOnly")),
              FilledNoAmp = n_distinct(which(Filled == TRUE & AmpPos == FALSE)),
              NOAMPwt = n_distinct(which(FluorQC == "wtFluor")),
              NOAMPmut = n_distinct(which(FluorQC == "mutFluor")),
              NOAMPboth = n_distinct(which(FluorQC == "Both")),
              TripleNeg = n_distinct(which(FluorQC == "Negative"))
              )
kable(statsFilledWells)
```

| Array | Template |  TotalFilled|  FilledAmped|   WT|  MUT|  HET|  ampOnly|  FilledNoAmp|  NOAMPwt|  NOAMPmut|  NOAMPboth|  TripleNeg|
|:------|:---------|------------:|------------:|----:|----:|----:|--------:|------------:|--------:|---------:|----------:|----------:|
| A1    | OCI-AML3 |          763|          212|   75|   81|   56|        0|          551|       15|        15|          0|        521|
| A2    | OCI-AML3 |          558|          180|   74|   67|   38|        1|          378|        9|         3|          0|        366|
| A3    | OCI-AML3 |          428|          115|   39|   46|   30|        0|          313|        2|         4|          1|        306|

``` r
write.csv(statsFilledWells, file = paste0(output,"/", runDate, ".statsFilledWells.csv"), row.names = F)
```

### Cell Specific Results

Now we can combine our zygosity measurements with our known cell locations, and find the zygosity of wells that contained single cells!

``` r
if(dir.exists("./OlympusCSVs")){
set1 <- brewer.pal(n = 9, "Set1")
set3 <- brewer.pal(n = 12, "Set3")
zygcolors <- c(set3[5], set3[10], set3[4], "black", "white")

zygBarPlotCells <- ggplot(data = snowCone1 %>% filter(Zygosity != "UNCALLED" & cellCount == 1),
                     aes(x = Zygosity, fill = Zygosity, group = Array), color = "black") +
                      geom_bar(aes(y = ..prop.., fill = factor(..x..))) +
                      ylim(0, 1.2) +
                      scale_fill_manual(values = zygcolors) + theme_bw() +
                      geom_text(aes(label = ..count.., y = ..prop..), stat = "count", vjust = -0.5) +
                      labs(x = "Zygosity", y = "Frequency") +
                      scale_x_discrete(limits=c("WT", "HET", "MUT")) +
                      theme(legend.position = "none") +
                      facet_grid(. ~ Array) +
                      theme(strip.text.y = element_text(size = 7, angle = 0))
print(zygBarPlotCells)

ggsave(paste0(output, "/zygBarPlot-singleCells.jpg"),
       zygBarPlotCells,
       width = 6, height = 6)
}
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Plot%20zygosities%20of%20cells-1.png"|absolute_url}})

**Figure 12. Bar plot of single-cell zygosities.**

We can also calculate some performance metrics for our device. For instance, we can calculate:
false positives: no cell was present but the reaction occured false negative: a cell was present but no reaction did not occur true positive: a cell was present and the reaction occured true negative: no cell was present and no reaction occured

``` r
if(dir.exists("./OlympusCSVs")){
failureData <- snowCone1 %>%
  filter(cellCount >= 0) %>%
  group_by(Array, Template) %>%
    summarise(FalsePositives = n_distinct(which(Filled == TRUE & AmpPos == TRUE & cellCount == 0)),
              FalseNegatives = n_distinct(which(Filled == TRUE & AmpPos == FALSE & cellCount == 1)),
              TruePositives = n_distinct(which(Filled == TRUE & AmpPos == TRUE & cellCount == 1)),
              TrueNegatives = n_distinct(which(Filled == TRUE & AmpPos == FALSE & cellCount == 0)),
              DoubletPlus = n_distinct(which(Filled == TRUE & cellCount > 1)),
              EmptyWells = n_distinct(which(Filled == FALSE)),
              FPR = FalsePositives/(FalsePositives + TrueNegatives),
              FNR = FalseNegatives/(FalseNegatives + TruePositives)
              )

failColors <- brewer.pal(6, "Accent")
failurePlotCounts <- melt(failureData %>% select(-c(EmptyWells, FPR, FNR)), variable.name = "Mode", value.name = "WellCount")
failurePlot <- ggplot(data = failurePlotCounts) +
  geom_bar(aes(x = Array, y = WellCount, fill = Mode), stat = "identity") +
  scale_fill_manual(values = failColors) +
  scale_y_continuous(expand = c(0, 10)) +
  theme_bw()
print(failurePlot)

ggsave(paste0(output, "/failurePlot.jpg"),
       failurePlot,
       width = 6, height = 6)

write.csv(failureData, file = paste0(output,"/", runDate, ".failureRateResults.csv"), row.names = F)
kable(failureData)
}
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Failure%20Rates-1.png"|absolute_url}})

| Array | Template |  FalsePositives|  FalseNegatives|  TruePositives|  TrueNegatives|  DoubletPlus|  EmptyWells|        FPR|        FNR|
|:------|:---------|---------------:|---------------:|--------------:|--------------:|------------:|-----------:|----------:|----------:|
| A1    | OCI-AML3 |             118|              44|             71|            496|           34|         261|  0.1921824|  0.3826087|
| A2    | OCI-AML3 |             118|              29|             56|            347|            8|         466|  0.2537634|  0.3411765|
| A3    | OCI-AML3 |              77|              12|             33|            298|            8|         596|  0.2053333|  0.2666667|

**Figure 13. Bar plot and table of device performance metrics.**

We can also see if our distribution of cells in the replicates was random. If it were random, we would expect the counts of empty, single cell, and multiple-cell wells to fit the Poisson density function.

``` r
if(dir.exists("./OlympusCSVs")){
experimentalLambda <- snowCone1 %>%
  filter(cellCount >= 0) %>% group_by(Array, Template) %>%
    summarise(TotalCounts =   n_distinct(which(Filled == TRUE)),
              ZeroCells =     n_distinct(which(Filled == TRUE & cellCount == 0)),
              OneCell =       n_distinct(which(Filled == TRUE & cellCount == 1)),
              DoubletPlus =   n_distinct(which(Filled == TRUE & cellCount > 1)),
              ExpLambda =    -log(ZeroCells/TotalCounts),
              PredSingles =  TotalCounts*ExpLambda*exp(-ExpLambda)/1,
              PredDoubPlus = TotalCounts - ZeroCells - PredSingles
              )

PoissonTest <- experimentalLambda %>% gather("Occupancy", "Counts", ZeroCells:DoubletPlus) %>% select(Array:ExpLambda, Occupancy, Counts) %>% mutate(Origin = "Actual")
PoissonTest2 <- experimentalLambda %>% mutate(PredZero = ZeroCells) %>% select(Array:TotalCounts, ExpLambda, "OneCell"= PredSingles, "DoubletPlus" = PredDoubPlus, "ZeroCells" = PredZero) %>% gather("Occupancy", "Counts", OneCell:ZeroCells) %>% mutate(Origin = "Predicted")
PoissonFull <- PoissonTest %>% bind_rows(PoissonTest2)                

PoissonFit <- ggplot(data = PoissonFull, aes(x = Occupancy, y = Counts, fill = Origin)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(limits=c("ZeroCells", "OneCell", "DoubletPlus")) +
  facet_grid(. ~ Array) +
  ggtitle("Actual and Poisson-Predicted Cell Counts") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1), strip.text.y = element_text(size = 7, angle = 0))
print(PoissonFit)
}
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Poisson%20Calculations-1.png"|absolute_url}})

**Figure 14. Actual counts of cell distributions versus those predicted by the Poisson density function.**

``` r
if(dir.exists("./OlympusCSVs")){
PoissonSummary <- PoissonFull %>% group_by(Origin, Occupancy) %>% summarise(Counts = sum(Counts))

PoissonSummaryPlot <- ggplot(data = PoissonSummary, aes(x = Occupancy, y = Counts, fill = Origin)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_x_discrete(limits=c("ZeroCells", "OneCell", "DoubletPlus")) +
  ggtitle("Combined Actual and Poisson-Predicted Cell Counts") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))
print(PoissonSummaryPlot)
}
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Poisson%20Summary-1.png"|absolute_url}})

**Figure 15. For sum of all replicates, actual counts of cell distributions versus those predicted by the Poisson density function.**

### Array Maps

Finally, we can veiw a representative map of the original images, to find locations of the zygosities and cells.

``` r
plotcolors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Zygositymapdf <- snowCone1 %>%
                            mutate(MapVal = if_else(AmpPos == FALSE & Filled == TRUE, "NonAmp",
                                            if_else(Filled == FALSE, "Empty",
                                                as.character(Zygosity))))

Zygositymapdf$MapVal <- factor(Zygositymapdf$MapVal, levels = c("WT", "HET", "MUT", "UNCALLED", "Empty", "NonAmp"))
Zygositymapdf$MapVal <- mapvalues(Zygositymapdf$MapVal,
                                  from = c("WT", "HET", "MUT", "UNCALLED", "Empty", "NonAmp"),
                                  to = c("Wild-Type", "Heterozygous", "Mutant", "Uncalled", "Empty", "No Amplification"))
zygosityColors <- c(plotcolors[6], plotcolors[4], plotcolors[7], "black", "gray95", "gray80")

Zygositymap <- ggplot() +
  geom_point(data = Zygositymapdf, aes(x = ArrayCol, y = ArrayRow, color = MapVal), size = 2) +
  labs(color = "Zygosity") +
  scale_color_manual(values = zygosityColors, limits = levels(Zygositymapdf$MapVal)) +
  scale_y_continuous(expand = c(0,1), trans = "reverse", name = "") +
  scale_x_continuous(name="")+
  facet_grid(Array ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.size = unit(0.5, "line"), legend.margin = margin(-20, 1, 0, 1),
        panel.grid = element_blank(), axis.ticks=element_blank(),
        axis.text = element_blank(), text = element_text(size = 8))
Zygositymap
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Zygosity%20Array%20Maps-1.png"|absolute_url}})

**Figure 16. Map of zygosity locations.**

``` r
if(dir.exists("./OlympusCSVs")){
Cellmapdf <- snowCone1 %>% filter(cellCount >= 0) %>%
  mutate(MapVal = if_else(cellCount == 0, "No Cells",
                                            if_else(cellCount == 1, "Single Cell",
                                            if_else(cellCount > 1, "Doublet Plus",
                                                    "unknown"))))

Cellmapdf$MapVal <- factor(Cellmapdf$MapVal, levels = c("No Cells", "Single Cell", "Doublet Plus"))
CellColors <- c("gray90", plotcolors[4], "black")

Cellmap <- ggplot() +
  geom_point(data = Cellmapdf, aes(x = ArrayCol, y = ArrayRow, color = MapVal), size = 2) +
  labs(color = "Cell Status") +
  scale_color_manual(values = c(CellColors, brewer.pal(9, "Blues")[9], "black"), limits = levels(Cellmapdf$MapVal)) +
  scale_y_continuous(expand = c(0,1), trans = "reverse", name = "") +
  scale_x_continuous(name="")+
  facet_grid(Array~.) +
  theme_bw() +
  theme(legend.position = "bottom", legend.margin = margin(-20, 0, 0, 0),
        panel.grid = element_blank(), axis.ticks=element_blank(),
        axis.text = element_blank(), text = element_text(size = 8))
print(Cellmap)
}
```

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Cell%20Array%20Maps-1.png"|absolute_url}})

**Figure 17. Map of single cell locations.**
