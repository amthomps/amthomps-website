---
title: "Alison Thompson - AnalyzeParticles"
layout: textlay
excerpt: "AnalyzeParticles"
sitemap: false
permalink: /AnalyzeParticles/
---

Introduction
------------

The purpose of this experiment was to determine the effect of Triton concentration on the fluorescence intensity of cells in PCR buffer containing 0.5X EvaGreen. EvaGreen stains DNA and can only penetrate cells with compromised membranes, so more Triton could lead to more permeabilization and better cell staining.

The cells were imaged on a hemocytometer at 1.25X zoom (2X objective, 0.63X c-mount). This is an example cell image:

``` r
require(dplyr); require(tidyr); require(ggplot2); require(knitr)
```

![]({{"/images/AnalyzeParticles_files/figure-markdown_github/Example_image_3_FAM_50ms.png"|absolute_url}})

**Figure 1. Raw fluorescence image of cells stained with EvaGreen on an Olympus MVX10 microscope.**

Images were run through the imageJ macro CellIntensityAnalysis which performs flat-fielding and particle analysis to report the size and intensity of each cell above a certain threshold. Here is an example of where the macro found particles meeting the selection criteria:

![]({{"/images/AnalyzeParticles_files/figure-markdown_github/Example_map_3_FAM_50ms.png"|absolute_url}})

**Figure 2. Map of particles measured by the ImageJ macro.**

ImageJ Analysis
---------------

The following code is writen using ImageJ's macro language. If all the macro code chunks are saved as a single .txt file, it can be run in ImageJ by first opening the image containing particles to be analyzed, then navigating on the imageJ toolbar to plugins -&gt; macros -&gt; run.

The first few lines reset some parameters to default for consistancy, and set the image name and directory as variables so that we can eventually deposit a .csv file with the same name as the image.

``` r
#first remove any global scale so that all size measurements are in units of pixels
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel global");
dir = File.directory;
name = getTitle();
namebase = replace(name, ".tif", "");
```

Next, we perform flat-fielding. This image was aquired on an Olympus MVX-10 at 1.6X magnification. This is a macro zoom microscope for imaging wide areas. The distorted intensity profile across the image is very noticable in these images - the edges are much dimmer than the center. I used a protocol similar to the one on this website for flat-field acquisition: <http://nic.ucsf.edu/resources/how-to-acquire-flat-field-correction-images/>

To perform the correction, we will first multiply every pixel intensity by a factor, then divide our image by the flat-field image to correct for the distortion. The average pixel intensity value in my flat-field image is about 50, so I chose to multiply the raw image by 50 before the correction. This allows me to maintain the intensity range in the raw image.

``` r
#mulitply the image by a factor of 50.  
#To flat field, we need to divide the image by the flat field image.
#We multiply all pixels by a factor first to maintain our dynamic range.
run("Multiply...", "value=50.000");
open(dir + "20180118-flat-field-zoom-125-test4.tif");
#Flat-fielding operation
imageCalculator("Divide", name,"20180118-flat-field-zoom-125-test4.tif");
```

We can now use ImageJ's built-in Analyze Particles function. This function requires binary images to operate. Pixels with intensity = 1 are particles, and intensity = 0 is background. We make this binary image by making a copy our cells image, applying a threshold, and converting everything above the threshold to 1 and below to 0. We set the size (pixels^2) and circularity (0-1, where 1 is most circular), and Analyze Particles defines an ROI for each particle that fits this criteria. We can then redirect these ROIs to our cell image and measure the size, integrated density, and mean intensity for each particle. The resulting csv is exported to the directory. We now have a .csv file with the same name as it's unprocessed .tif file.

``` r
selectWindow(name);
#Duplicate the image, find the location of particles in this copy, and measure values at the mapped locations in the original image.
run("Duplicate...", " ");
#a threshold of 805 seems appropriate for these images -
#above 805 is cell material, below 805 is background after flat-fielding.  
#4095 is the max value for these images, which were aquired with a 12-bit camera.
setThreshold(805, 4095);
setOption("BlackBackground", false);
#convert the image copy to binary - particle finder's requirement.
run("Convert to Mask");
#set measurements to particle area, particle integrated density, and particle median.
run("Set Measurements...", "area integrated median redirect=" + name +" decimal=3");
#Set the particle parameters.  
#Currently we are looking for particles between 2-200 square pixels,
#with 0.8-1.00 circularity.
run("Analyze Particles...", "size=2-200 circularity=0.80-1.00 show=Outlines display");
saveAs("results",  dir + namebase + ".csv");
```

Cell Intensity Analysis in R
----------------------------

In this experiment, I wanted to see if the intensity of the cells changed when the Triton concentration was varied between 0.02%, 0.03%, and 0.05%. Triton should allow the dye to penetrate the cell membrane more easily, so one might expect an increase in cell intensity with increasing amounts of triton to a certain limit. One image containing a few thousand cells was captured and analyzed at each condition. Once all images were processed with the ImageJ macro above, the following R code was used to pull in all csv's in the folder for analysis and visualization.

``` r
#find all csv's in the current directory
fileNames <- list.files(path = ".", pattern = "*.csv")
DataList <- lapply(fileNames, read.csv, header = TRUE)
ImageTitle <- substring(gsub(".csv", "", fileNames), 1, 100)
names(DataList) <- ImageTitle
#make a data frame containing one row for every particle measured
AllImages <- bind_rows(DataList, .id = "ImageTitle")
```

The conditions and exposure times were part of the filenames, so I pull those out and include them in the data frame.

``` r
#separate the image title column into the experimental condition "condition", and the camera exposure time "exposure".
AllImages <- separate(AllImages, ImageTitle, sep = "_FAM_", into = c("Condition", "Exposure"))
#keep only the measurement columns for particle area, integrated density "IntDen", and the median intensity.
AllImages <- select(AllImages, Condition, Exposure, Area, IntDen, Median)
AllImages <- mutate(AllImages, Mean = IntDen/Area)
#Decode the condition identity - in this experiment the condition # corresponds to the Triton concentration (as % w/v Triton)
AllImages <- mutate(AllImages, Triton = case_when(
  Condition == 1 ~ 0.02,
  Condition == 2 ~ 0.03,
  Condition == 3 ~ 0.05
))
```

``` r
HistPart50 <- ggplot(AllImages %>% filter(Exposure == "50ms"), aes(x = Mean)) +
  geom_histogram() +
  xlab("Cell Intensity") +
  ylab("Cell Count") +
  facet_grid(Triton ~ .) +
  theme_bw() +
  theme(panel.grid = element_blank())
HistPart50
```

![]({{"/images/AnalyzeParticles_files/figure-markdown_github/plot%20particle%20histograms-1.png"|absolute_url}})

**Figure 3. Histogram of particle intensity at each condtion.** Cells were imaged in 0.5X EvaGreen in various concentrations of Triton (0.02%, 0.03%, 0.05% w/v).

``` r
ViolinPart50 <- ggplot(AllImages %>% filter(Exposure == "50ms"),
                     aes(y = Mean, x = Triton, fill = Condition)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = TRUE) +
  geom_hline(yintercept = 805, color = "orchid4") +
  ylim(0, 1700) +
  xlim(0, 0.06) +
  ylab("Cell Intensity") +
  xlab("% Titon Concentration") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank())
ViolinPart50
```

![]({{"/images/AnalyzeParticles_files/figure-markdown_github/plot%20violin%20distributions-1.png"|absolute_url}})

**Figure 4. Violin plot of particle intensity at each condtion.** Same data from the histogram in Figure 3. Cells were imaged in 0.5X EvaGreen in various concentrations of Triton (w/v%). Purple line is at the threshold for cell detection, so any cells with intensity below this threshold were not measured by ImageJ. Lines within the violin plot represent the 1st, 2nd, and 3rd quartiles. In a violin plot, the width of the colored space correlates to the number of cells with that intensity.

**Table 1. Median particle intensity.**

``` r
#compute the median particle intensity
ImgSummary <- AllImages %>% filter(Exposure == "50ms") %>%
  group_by(Condition, Triton) %>%
  summarize(Med_Intensity = median(Mean),
            N = n())
kable(ImgSummary)
```

| Condition |  Triton|  Med\_Intensity|     N|
|:----------|-------:|---------------:|-----:|
| 1         |    0.02|        1094.333|  2013|
| 2         |    0.03|        1173.786|  3274|
| 3         |    0.05|        1340.354|  2840|

Remarks
-------

For the project this is a part of, we want the cells to be as bright as possible without increasing the dye concentration. The results of this experiment seem to indicate that increasing Triton does raise the average fluorescence intensity of our cells. It might be interesting to find the limits to this, but other factors limit the amount of Triton we can use. Nevertheless, an interesting result from one very quick experiment is a good thing, and more experiments will follow to further optimize cell brightness in our system.
