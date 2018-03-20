---
title: "Alison Thompson - TyphoonDataReduction"
layout: textlay
excerpt: "TyphoonDataReduction"
sitemap: false
permalink: /TyphoonDataReduction/
---

Instructions
------------

1.  Click "Knit" icon in menu bar to execute all code
2.  Change variables in code if necessary

Directory Setup and Data Cleaning
---------------------------------

Initialize libraries

``` r
require(plyr); require(dplyr); require(tidyr); require(reshape2); require(ggplot2)
require(ggExtra); require(ggrepel); require(gridExtra);  require(RColorBrewer); require(knitr)
options(stringsAsFactors = FALSE)
```

Create an output folder in the working directory and read the Typhoon metadata file

``` r
#Sets up the directory structure for where to look for data, and where to save results.  
#This assumes that the raw Typhoon and Olympus (if cells) data has been processed using the ImageJ macros.  This .rmd must run from a processed data folder within grp/SDGenotypingPhase2/ImageAnalysis
runDate <- format(Sys.time(), "%Y%m%d-%H%M")
output <- c(paste0("./OutputData-",runDate))
dir.create(output, showWarnings = FALSE)
metadataDir <- c("./Metadata")
fluorInput <- c("./TyphoonCSVs")
writeLines(capture.output(sessionInfo()),
           paste(output,"/TyphoonDR-sessionInfo.txt", sep = ""))

#Allows us to read in the probeScheme to determine which fluor is amp/wt/mut and also read in which template is in each array.
arrayList <- c("R01", "R02", "A1", "A2", "A3")
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

The probe scheme used for this chip is **PS1**.
The amplification probe is in the **FAM** channel, the wild-type probe is in the **HEX** channel, and the mutant probe is in the **Cy5** channel. The templates used were **BS, BS, HET, WT, WT**. The run outputs will be saved in the folder **./OutputData-20180319-2103**.

Create a single dataframe from all Typhoon images. Label arrays and well rows and columns

``` r
#Read in the individual csv's for each fluorphore to make a single data frame.  
fileNames <- list.files(path = fluorInput, pattern = "*.csv", full.names = TRUE)
dataList <- lapply(fluorList, function(x)
    {read.csv(file = fileNames[grep(x, fileNames)], stringsAsFactors = FALSE)})
names(dataList) <- fluorList
dataList <- lapply(dataList, function(x) {colnames(x)[1] <- "Well";x})
allFluors <- melt(dataList, id=colnames(dataList[[1]]))
colnames(allFluors)[length(colnames(allFluors))] <- "Fluor"
names(allFluors)[names(allFluors) == "X.Area"] <- "Area"

rm(dataList)

#This is for the design with two control arrays of 256 wells and three arrays of 1024 wells, appearing in the Typhoon image top-to-bottom in that order.  Change the numbers if using a new design!
allFluors$Array[allFluors$Well <= 256] <- "R01"
allFluors$Array[allFluors$Well > 256 & allFluors$Well <= 512] <- "R02"
allFluors$Array[allFluors$Well > 512 & allFluors$Well <= 1536] <- "A1"
allFluors$Array[allFluors$Well > 1536 & allFluors$Well <= 2560] <- "A2"
allFluors$Array[allFluors$Well > 2560 & allFluors$Well <= 3584] <- "A3"

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
primeTime$ArrayRow <- c(rep(1:4, each = 64, times = 2), rep(1:16, each = 64, times = 3))
primeTime$ArrayCol <- seq(1:64)

rm(procData, ampArea, wideForm)
```

``` r
ggplot(data = allFluors, aes(x=Well, y=AID)) + geom_point() + facet_grid( ~ Fluor)
```

![](TyphoonDataReduction_files/figure-markdown_github/histogram%20of%20AID%20by%20Array%20and%20Color-1.png)

**Figure. Scatter plot of the AIDs vs. well \# for each fluorophore.**

If there is Olympus data in the working directory, read in cell data

``` r
if(dir.exists("./OlympusCSVs")){

cellInput <- c("./OlympusCSVs")
filenamesOlympus <- list.files(path = cellInput, pattern = ".csv", full.names = T)
#images <- as.numeric(str_sub(gsub(".csv", "", filenamesOlympus), -2,-1))
cellData <- lapply(filenamesOlympus, function(x) read.csv(x, header = T))
#names(cellData) <- images
cellMelt <- melt(cellData, id = colnames(cellData[[1]]))
colnames(cellMelt) <- c("imageWell", "cellCount", "imageId")
cellMelt$ROIorder <- seq(1:3072)

lookup <- data.frame(ROIorder = 1:3072)
lookup$ArrayCol <- rep(c(rep(1:11, len = 88), rep(12:22, len = 88),
                         rep(23:33, len = 88), rep(34:44, len = 88),
                         rep(45:55, len = 88), rep(56:64, len = 144),
                         rep(45:55, len = 88), rep(34:44, len = 88),
                         rep(23:33, len = 88), rep(12:22, len = 88),
                         rep(1:11, len = 88)), 3)
lookup$ArrayRow <- rep(c(rep(1:8, each = 11, len = 440),
                         rep(1:8, each = 9),rep(9:16, each = 9),
                         rep(9:16, each = 11, len = 440)),3)
#lookup$estImageId <- c(rep(1:5, each = 88), rep(6:7, each = 72),
#                       rep(8:17, each = 88), rep(18:19, each = 72),
#                       rep(20:29, each = 88), rep(30:31, each = 72),
#                       rep(32:36, each = 88))
lookup$Array <- rep(c("A1", "A2", "A3"), each = 1024)
lookup <- arrange(lookup, Array, ArrayRow, ArrayCol)
#since array1 starts at well 513 (two 256 well arrays R01 and R02 come first in Typhoon image)
lookup$Well <- seq(from = 513, to = 3584)
fullMergeData <- join(lookup, cellMelt, by = "ROIorder", type = "left")
cellsReadyToRoll <- select(fullMergeData, Well, cellCount)

primeTime <- left_join(primeTime, cellsReadyToRoll, by = "Well")
primeTime$cellCount <- as.numeric(primeTime$cellCount)
}
```

Thresholding
------------

``` r
tArea <- 27
ggplot(data = primeTime %>% filter(Area > 5), aes(x=Area)) +
    geom_density() + facet_grid(. ~ Array) +
    geom_vline(xintercept = tArea, color = "springgreen4") +
    theme_bw()
```

![](TyphoonDataReduction_files/figure-markdown_github/area%20threshold-1.png)

``` r
ruleThemAll <- primeTime %>% mutate(Filled = Area > tArea)
```

**Figure. Well Area Threshold. A density curve of well area for each array of this chip is shown with the area threshold drawn in green at 27. All wells below the threshold are empty or partially filled. These individual wells will be marked as "Filled = FALSE" and will not be evaluated in results. This threshold can be adjusted if necessary in the code chunk above.**

``` r
ggplot(data = ruleThemAll %>% filter(Filled == TRUE), aes(x=ArrayCol, y=ampFluor)) +
  geom_point() +
  facet_grid( ~ Array)
```

![](TyphoonDataReduction_files/figure-markdown_github/histogram%20of%20AID%20by%20Array%20and%20Color%202-1.png)

``` r
ruleThemAll <- primeTime %>% group_by(Array) %>% mutate(Filled = ArrayCol < 63 & ArrayCol > 2 &
                                                          ArrayRow < 16 & ArrayRow > 1 &
                                                          Area > tArea)
```

**Figure. Adjusted integrated denisity (AID) of filled wells. A scatter plot of AID vs. column number is shown to determine effects of evaporation. A pattern of high values for the first and last columns indicates evaporation during thermalcycling. By default, the first and last column are marked as "unfilled" to remove them from analysis. Columns can be selected for removal in the above code chunk.**

Should any arrays be excluded from the run due to experimental error? Examples include air entering the array during fill or cycling, or cross-talk during thermalcycling. Arrays can be selected for removal in this code chunk. Alternatively, templates can be labeled "BS" in the runInfo.txt file in the folder containing the original Typhoon images.

``` r
#Are any of the arrays "BS"?  AKA filling malfunction, had lots of cross-talk during PCR.  Un-# BS arrays below
#ruleThemAll$Template[ruleThemAll$Array == "R01"] <- "BS"
#ruleThemAll$Template[ruleThemAll$Array == "R02"] <- "BS"
#ruleThemAll$Template[ruleThemAll$Array == "A1"] <- "BS"
#ruleThemAll$Template[ruleThemAll$Array == "A2"] <- "BS"
#ruleThemAll$Template[ruleThemAll$Array == "A3"] <- "BS"

ruleThemAll <- ruleThemAll %>% filter(Template != "BS")
```

``` r
chipID <- substring(gsub(".*amthomps_", "", gsub("-Results.*", "", getwd())), 1, 100)
tagname <- paste0("-NoThresh-", chipID, "-", pSch, ".csv")
write.csv(ruleThemAll, file = paste0(output, "/", runDate, tagname), row.names = F)
```

``` r
avgHEXSlope = NA
if(pSch == "PS1" & "MUT" %in% ruleThemAll$Template){
  print(ggplot(ruleThemAll %>% filter(Filled == TRUE, Template == "MUT"),
       aes(x = ampFluor, y = wtFluor, color = Array)) +
         geom_point() +
         ggtitle("HEX Bleed-Through") +
          geom_smooth(method = "lm"))
  HEXSlopes <- ruleThemAll  %>% filter(Filled == TRUE, Template == "MUT") %>%
    group_by(Array)%>% do(model = lm(wtFluor ~ ampFluor, data=.)) %>%
    mutate(Slope = coef(model)[2]) %>% select(Array, Slope)
  avgHEXSlope <-mean(HEXSlopes$Slope)
}

if(pSch == "ZygSwap" & "WT" %in% ruleThemAll$Template){
  print(ggplot(ruleThemAll %>% filter(Filled == TRUE, Template == "WT"),
        aes(x = ampFluor, y = mutFluor, color = Array)) +
         geom_point() +
          ggtitle("HEX Bleed-Through") +
          geom_smooth(method = "lm"))

  HEXSlopes <- ruleThemAll  %>% filter(Filled == TRUE, Template == "WT") %>%
    group_by(Array) %>% do(model = lm(mutFluor ~ ampFluor, data=.)) %>%
    mutate(Slope = coef(model)[2]) %>% select(Array, Slope)
  avgHEXSlope <-mean(HEXSlopes$Slope)
}

if(pSch == "AmpSwap" & "MUT" %in% ruleThemAll$Template){
  print(ggplot(ruleThemAll %>% filter(Filled == TRUE, Template == "MUT"),
        aes(x = mutFluor, y = wtFluor, color = Array)) +
         geom_point() +
          ggtitle("HEX Bleed-Through") +
          geom_smooth(method = "lm"))
  HEXSlopes <- ruleThemAll  %>% filter(Filled == TRUE, Template == "MUT") %>%
    group_by(Array)%>% do(model = lm(wtFluor ~ mutFluor, data=.)) %>%
    mutate(Slope = coef(model)[2]) %>% select(Array, Slope)
  avgHEXSlope <-mean(HEXSlopes$Slope)
}
```

**Figure. For templates off-target for the HEX fluorophore were used, scatterplot of the HEX fluorophore vs. the FAM fluorophore to determine if fluorescence bleed-through is occurring.**

The average HEX slope between the arrays analyzed is:
**NA**.

We currently are not correcting for this HEX bleed-through. If you want to correct the data, edit the header of this chunk to `eval=TRUE`.

``` r
#Modify dfAllControls to use HEX correction, where m is avgHEXSlope, HEXFluor(adj) = HEXFluor(initial) - m*FAMFluor
if(pSch == "PS1") {ruleThemAll <- mutate(ruleThemAll, wtFluor = wtFluor - avgHEXSlope*ampFluor)}
if(pSch == "ZygSwap") {ruleThemAll <- mutate(ruleThemAll, mutFluor = mutFluor - avgHEXSlope*ampFluor)}
if(pSch == "AmpSwap") {ruleThemAll <- mutate(ruleThemAll, wtFluor = wtFluor - avgHEXSlope*mutFluor)}

#View results post-HEX correction
if(pSch == "PS1") {ggplot(ruleThemAll %>% filter(Filled == TRUE, Template == "MUT"),
       aes(x = ampFluor, y = wtFluor, color = Array)) +
         geom_point() +
          geom_smooth(method = "lm")
}

if(pSch == "ZygSwap") {ggplot(ruleThemAll %>% filter(Filled == TRUE, Template == "WT"),
        aes(x = ampFluor, y = mutFluor, color = Array)) +
         geom_point() +
          geom_smooth(method = "lm")
}

if(pSch == "AmpSwap") {ggplot(ruleThemAll %>% filter(Filled == TRUE, Template == "MUT"),
        aes(x = ampFluor, y = mutFluor, color = Array)) +
         geom_point() +
          geom_smooth(method = "lm")
}
```

kmeans, max density curve for first cluster (currently not running this chunk, testing the next kmeans chunk as of 3/14/18)

kmeans, cluster centers

``` r
#function to find cluster value of first of two clusters for each array and fluorophore
NegPeak <- function(A, B) {
kmeansdf <- ruleThemAll %>% gather(Fluor, AID, ampFluor:mutFluor) %>%
  filter(Filled == TRUE, Array == A, Fluor == B);
return(min(kmeans(kmeansdf$AID, centers=2, nstart = 25)[["centers"]]));
}

#call to function
VArrays <- levels(droplevels(ruleThemAll$Array))
VFluors <- c("ampFluor", "wtFluor", "mutFluor")
#parameters to iterate over
ParamA <- rep(VArrays, each = 3)
ParamB <- rep(VFluors, times = length(VArrays))
snowPeaks <- mapply(NegPeak, ParamA, ParamB)
snowPeaks <- data.frame(Array=names(snowPeaks), maxAID=snowPeaks, row.names=NULL)
snowPeaks$Fluor <- ParamB
```

``` r
#using previously determined maxSD of negative controls, draw thresholds at 6 std deviations above the mean of the negative wells.  Edit this chunk if thresholds should be adjusted.
#Max SD's x multiplier for each probe from negative templates (NTC or off-target plasmid).  
#See "Z:\grp\SDChipImageAnalysis\Multi-FWHMAnalysis-20170919-1439\20170919-1439SDstatsNoFilter_totals.xlsx" for calculations.
#change multiplier (xsigamp/mut/wt) to the x sigma you want to use for thresholding.  maxAID + xsig*sigma.
xsigamp <- 6
xsigmut <- 6
xsigwt <- 6
PS1ampSD      <- 1098*xsigamp
PS1mutSD      <- 1280*xsigmut
ZygSwapampSD    <- 1098*xsigamp
ZygSwapwtSD     <- 788*xsigwt
PS1wtSD       <- 1059*xsigwt
ZygSwapmutSD    <- 2086*xsigmut
AmpSwapampSD  <- 1098*xsigamp
AmpSwapmutSD  <- 1278*xsigmut
AmpSwapwtSD  <- 1098*xsigwt


#Add a threshold to the maxAID table as a simple "+ constant" for everything
#snowDrift <- snowPeaksAdj %>% mutate(Threshold = maxAID + 400000)

#Add thresholds to snowPeaksAdj
snowDrift <- snowPeaks %>% mutate(Scheme = pSch) %>%
  mutate(Threshold = if_else(Scheme == "PS1" & Fluor == "ampFluor", maxAID + PS1ampSD,
                    if_else(Scheme == "PS1" & Fluor == "wtFluor", maxAID + PS1wtSD,
                    if_else(Scheme == "PS1" & Fluor == "mutFluor", maxAID + PS1mutSD,
                    if_else(Scheme == "ZygSwap" & Fluor == "ampFluor", maxAID + ZygSwapampSD,
                    if_else(Scheme == "ZygSwap" & Fluor == "wtFluor", maxAID + ZygSwapwtSD,
                    if_else(Scheme == "ZygSwap" & Fluor == "mutFluor", maxAID + ZygSwapmutSD,
                    if_else(Scheme == "AmpSwap" & Fluor == "ampFluor", maxAID + AmpSwapampSD,
                    if_else(Scheme == "AmpSwap" & Fluor == "wtFluor", maxAID + AmpSwapwtSD,
                    if_else(Scheme == "AmpSwap" & Fluor == "mutFluor", maxAID + AmpSwapmutSD,
                            0))))))))))

#Add this threshold to the df with all the AID's
snowDrift1 <- snowDrift %>% select(Array, Fluor, Threshold) %>% spread(Fluor, Threshold) %>% select(Array, ampThresh = ampFluor, mutThresh = mutFluor, wtThresh = wtFluor)

snowCone1 <- left_join(ruleThemAll, snowDrift1, by = c("Array"))
```

    ## Warning: Column `Array` joining factor and character vector, coercing into
    ## character vector

Calculate fluorescence thresholds (currently not running, testing new kmeans method as of 3/14)

``` r
ampThreshPlot <- ggplot(snowCone1 %>% filter(Filled == TRUE), aes(x=ampFluor)) +
                            geom_histogram(aes(x = ampFluor, y = ..density..), bins = 30) +
                            geom_point(stat = "density", alpha = 0.6, stroke = 0) +
                            geom_vline(aes(xintercept = ampThresh), color = "orchid4") +
                            xlim(0, 120000) +
                            xlab("AID") + ylab("Density") +
                            facet_grid(. ~ Array) +
                            ggtitle("Amplification Threshold") +
                            theme_bw() +
                            theme(strip.text.y = element_blank(), axis.text.x = element_text(angle=90, hjust=1))
ampThreshPlot
```

![](TyphoonDataReduction_files/figure-markdown_github/view%20amp%20thresh-1.png)

**Figure. Threshold for amplification fluorophore. A histogram and density curve of ampFluor AID is shown. Wells with ampFluor AID left of the threshold will be marked as amplification negative for results.**

``` r
ampThreshPlot2 <- ggplot(snowCone1 %>% filter(Filled == TRUE), aes(x=ampFluor)) +
                            geom_histogram(aes(x = ampFluor), bins = 30) +
                            geom_vline(aes(xintercept = ampThresh), color = "orchid4") +
                            xlim(0, 120000) +
                            xlab("AID") + ylab("Density") +
                            facet_grid(. ~ Array) +
                            ggtitle("Amplification Threshold") +
                            theme_bw() +
                            theme(strip.text.y = element_blank(), axis.text.x = element_text(angle=90, hjust=1))
ampThreshPlot2
```

![](TyphoonDataReduction_files/figure-markdown_github/view%20amp%20thresh%20without%20density%20curve-1.png)

**Figure. Threshold for amplification fluorophore without denisty curve. A histogram-only version of the previous plot.**

``` r
mut_vs_wt <- ggplot(snowCone1 %>% filter(Filled == TRUE),
            aes(x=wtFluor, y=mutFluor)) +
            geom_point(alpha = 0.6, stroke = 0, aes(color = Array)) +
            xlim(0,120000) +
            ylim(0,120000) +
            coord_fixed() +
            labs(x = "wtFluor", y = "mutFluor") +
            facet_wrap("Array", nrow = 2) +
            ggtitle("Wild-Type and Mutant Thresholds") +
            geom_vline(aes(xintercept = wtThresh), color = "orchid4") +
            geom_hline(aes(yintercept = mutThresh), color = "orchid4") +
            theme_bw() +
            theme(strip.text.y = element_blank(), axis.text.x = element_text(angle=90, hjust=1))
mut_vs_wt
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](TyphoonDataReduction_files/figure-markdown_github/Preliminary%20View%20of%20wt%20and%20mut%20Thresholds-1.png)

**Figure. Thresholds for mutant and wild-type fluorophore. A scatter plot of mutant vs. wild-type AIDs for each filled well is shown with the generated mutant and wild-type thresholds. These quadrants will be used to determine zygosity in each well.**

Results
-------

Determine Zygosities in the individual arrays

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

    ## Warning: Removed 114 rows containing missing values (geom_point).

![](TyphoonDataReduction_files/figure-markdown_github/View%20Results%20of%20wt%20and%20mut%20Thresholding-1.png)

**Figure. Zygosity of each well that passed quality control filtering. Color indicates zygosity.**

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

![](TyphoonDataReduction_files/figure-markdown_github/View%20Results%20Fluor%20QC-1.png)

**Figure. Failure assignment of any wells that was filled but not assigned a zygosity.**

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
| A1    | HET      |          725|          319|   37|   37|  242|        3|          406|        0|         0|          0|        406|
| A2    | WT       |          773|          430|  394|    0|   36|        0|          343|        0|         0|          0|        343|
| A3    | WT       |          667|          402|  334|    0|   67|        1|          265|        4|         0|          0|        261|

``` r
write.csv(statsFilledWells, file = paste0(output,"/", runDate, ".statsFilledWells.csv"), row.names = F)
```

### HET Specific Results

``` r
if("HET" %in% snowCone1$Template){

set1 <- brewer.pal(n = 9, "Set1")
set3 <- brewer.pal(n = 12, "Set3")
zygcolors <- c(set3[5], set3[10], set3[4], "black")

zygBarPlot <- ggplot(data = filter(snowCone1, Zygosity != "UNCALLED", Template == "HET"),
                     aes(x = Zygosity, fill = Zygosity, group = Array), color = "black") +
    geom_bar(aes(y = ..prop.., fill = factor(..x..))) +
    scale_fill_manual(values = zygcolors) + theme_bw() +
    geom_text(aes(label = ..count.., y = ..prop..), stat = "count", vjust = -0.5) +
    labs(x = "Zygosity", y = "Frequency") + theme(legend.position = "none") #+
#    facet_grid(. ~ Array)
print(zygBarPlot)


totalSampleSize <- snowCone1 %>% group_by(Array, Filled) %>% filter(Filled == TRUE, Template == "HET") %>%
    summarize(TotalCounts = n()) %>% ungroup() %>% select(Array, TotalCounts)

zeroK <- snowCone1 %>% group_by(Array, Filled, AmpPos) %>% filter(Filled == TRUE, AmpPos == FALSE, Template == "HET") %>% summarize(ZeroCounts = n()) %>% ungroup() %>%select(Array, ZeroCounts)

dropoutData <- left_join(totalSampleSize, zeroK, by = c("Array"))

#Use Poisson probability mass function pmf = (lamba^k*e^-lamba)/k! to calculate experimental lamba
dropoutData <- mutate(dropoutData, expZeroFreq = ZeroCounts/TotalCounts)
dropoutData <- mutate(dropoutData, expLambda = -log(expZeroFreq))
dropoutData <- mutate(dropoutData, predOneFreq = expLambda*exp(-expLambda)/1)
dropoutData <- mutate(dropoutData, predTwoPlusFreq = 1-expZeroFreq-predOneFreq)


zygCounts <- snowCone1  %>% group_by(Array, Zygosity) %>% filter(Template == "HET") %>%
    summarize(Counts = length(Zygosity)) %>% filter(Zygosity != "UNCALLED")
zygCounts[is.na(zygCounts)] <- 0

dropoutData <- left_join(dropoutData, spread(zygCounts, Zygosity, Counts), by = c("Array"))
dropoutData[is.na(dropoutData)] <- 0

dropoutData <- mutate(dropoutData, singleHets = HET-round(predTwoPlusFreq*TotalCounts))
dropoutData <- mutate(dropoutData, mutDrop = WT/(WT+MUT+singleHets))
dropoutData <- mutate(dropoutData, wtDrop = MUT/(WT+MUT+singleHets))

ggsave(paste0(output, "/zygBarPlot-dropoutTable.jpg"),
       zygBarPlot,
       width = 6, height = 6)
write.csv(dropoutData, file = paste0(output,"/", runDate, ".dropoutData.csv"), row.names = F)
}
```

![](TyphoonDataReduction_files/figure-markdown_github/ZygBarPlot%20and%20dropoutTable-1.png)

**Figure. Bar chart of zygosities of HET wells containing one copy of plasmid.**

**Table. Drobout Data.**

``` r
if("HET" %in% snowCone1$Template){
kable(dropoutData)
}
```

| Array |  TotalCounts|  ZeroCounts|  expZeroFreq|  expLambda|  predOneFreq|  predTwoPlusFreq|   WT|  HET|  MUT|  singleHets|    mutDrop|     wtDrop|
|:------|------------:|-----------:|------------:|----------:|------------:|----------------:|----:|----:|----:|-----------:|----------:|----------:|
| A1    |          725|         406|         0.56|  0.5798185|    0.3246984|        0.1153016|   37|  242|   37|         158|  0.1594828|  0.1594828|

**Table. Dropout Summary.**

``` r
if("HET" %in% snowCone1$Template){
#only useful if you have more than one array with Template = HET.  If only one HET array, will give NaN for standard deviation.
dropoutSummary <- dropoutData %>% group_by(Array) %>% summarise(AvgWtADO = mean(wtDrop),
                                                                 SDWtADO = sd(wtDrop),
                                                                 AvgMutADO = mean(mutDrop),
                                                                 SdMutADO = sd(mutDrop),
                                                                 n = n())

write.csv(dropoutSummary, file = paste0(output,"/", runDate, ".dropoutSummary.csv"), row.names = F)
kable(dropoutSummary)
}
```

| Array |   AvgWtADO|  SDWtADO|  AvgMutADO|  SdMutADO|    n|
|:------|----------:|--------:|----------:|---------:|----:|
| A1    |  0.1594828|      NaN|  0.1594828|       NaN|    1|

### Cell Specific Results

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
       zygBarPlot,
       width = 6, height = 6)
}
```

**Figure. Bar plot of single-cell zygosities.**

Failure rate calulations

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
              FDR = FalsePositives/(FalsePositives + TruePositives),
              FPR = FalsePositives/(FalsePositives + TrueNegatives),
              FNR = FalseNegatives/(FalseNegatives + TruePositives)
              )

failColors <- brewer.pal(6, "Accent")
failurePlotCounts <- melt(failureData %>% select(-c(FDR, FPR, FNR)), variable.name = "Mode", value.name = "WellCount")
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

**Figure. Bar plot of failure counts.**

Poisson Distribution calculations

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

cellNumPlot <- melt(experimentalLambda %>% select(Array, ZeroCells, OneCell, DoubletPlus), variable.name = "Mode", value.name = "CellCount")
CellCountPlot <- ggplot(data = cellNumPlot) +
  geom_bar(aes(x = Array, y = CellCount, fill = Mode), stat = "identity") +
  theme_bw()
print(CellCountPlot)
}
```

Compare estimated and experimental zero/single/doublet+ counts, per array

``` r
if(dir.exists("./OlympusCSVs")){

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

Compare estimated and experimental zero/single/doublet+ counts, all data added together

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

Save a .jpeg of the Poisson plots

``` r
if(dir.exists("./OlympusCSVs")){
ggsave(PoissonFit, filename = paste0(output,"/", runDate, "PoissonFit.jpeg"))
ggsave(PoissonSummaryPlot, filename = paste0(output,"/", runDate, "PoissonSummaryPlot.jpeg"))
ggsave(CellCountPlot, filename = paste0(output,"/", runDate, "CellCountPlot.jpeg"))

write.csv(experimentalLambda, file = paste0(output,"/", runDate, ".poissionEstimateResults.csv"), row.names = F)
write.csv(PoissonFull, file = paste0(output,"/", runDate, ".PoissonFull.csv"), row.names = F)
write.csv(PoissonSummary, file = paste0(output,"/", runDate, ".PoissonSummary.csv"), row.names = F)
}
```

``` r
write.csv(snowCone1, file = paste0(output,"/", runDate, ".cellsnowCone.csv"), row.names = F)
```

### Array Maps

``` r
#Colorblind color pallet
#1=black, 2=orange, 3=sky-blue, 4=bluish-green, 5=yellow, 6=blue, 7=vermillion, 8=reddish-purple
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

![](TyphoonDataReduction_files/figure-markdown_github/Zygosity%20Array%20Maps-1.png)

**Figure. Map of array zygosities.**

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

**Figure. Map of single cell zygosities.**

**If you made changes to the rmd, save the file.**
