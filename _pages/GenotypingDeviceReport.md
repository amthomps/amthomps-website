---
title: "Alison Thompson - GenotypingDeviceReport"
layout: textlay
excerpt: "GenotypingDeviceReport"
sitemap: false
permalink: /GenotypingDeviceReport/
---

Introduction
------------

This code was used for a project to examine the genetics of single cells. The experiments use a microfluidic device to partition single cells into small wells where a reaction can be carried out. This is an picture of the device:

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/IMG_photo_of_SDChip_at_FH_crop.png"|absolute_url}})  
**Figure 1. Photo of the device with three replicate arrays.**

This device has three replicate arrays, each will be loaded with a solution containing cells and the reagents for the reaction. Each replicate array has 1024 wells, and each well will hold about 5 nanoliters of solution. The device shown here is made on a 2"x3" microscope slide. As the device is loaded, the cells are distributed randomly, so we first image the arrays to locate wells with single cells.

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/CellImage1.png"|absolute_url}})  
**Figure 2. Image of cells in wells.**

This image is acquired using a fluorescence microscope, and because the cells are stained with a fluorescent dye, we can see them as the tiny white spots in the image. The cells are about 15 micrometers in diameter, and the wells are 200 micrometers wide and 400 micrometers tall. 34 images are captured to cover all the device area.

After imaging, the reaction (PCR) takes about 2 hours. The reaction will only occur if a specific genetic sequence is present in the well. If the reaction does occur, we will observe an increased fluorescence signal in that well as a dye is released. We are looking for three different genetic sequences with three different dyes, so we image the device in three channels.

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Typhoon_composite.png"|absolute_url}})  
**Figure 3. Image of fluorescence post-reaction.**

From these 37 images, we will have enough information to report single-cell genotypes in our three replicates! To do this, we need to know which wells contain single cells, and whether they are negative or positive for the control sequence, the wild-type sequence, and/or the mutant sequence. The cells can be either wild-type, mutant, or heterozygous.

Experimental Details
--------------------

The images are processed using ImageJ using a macro language to overlay a grid of ROI's on the wells and extract the critical data. These include the number of cell-sized particles per well and the integrated fluorescence value for each well. ImageJ exports these values as .csv's. We are going to merge the 37 csv's into a single data frame using R.

**Experimental Conditions:**
The probe scheme used for this chip is **ZygSwap**.
The amplification probe is in the **FAM** channel, the wild-type probe is in the **Cy5** channel, and the mutant probe is in the **HEX** channel. The templates used were **OCI-AML3, OCI-AML3, OCI-AML3** in replicates A1, A2, and A3, respectively. The run outputs will be saved in the folder **./OutputData-20180407-1406**.

*note: in this example all replicates contain OCI-AML3, a cell line.*

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/scatter%20plot%20of%20AID%20by%20Well,%20where%20each%20panel%20is%20one%20channel-1.png"|absolute_url}})

**Figure 4. Scatter plot of the AIDs (fluorescence intensity) for each replicate array and fluorophore before filtering and corrections.**

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/histogram%20of%20AID,%20for%20each%20array%20and%20fluorophore-1.png"|absolute_url}})

**Figure 5. Histogram of the AIDs (fluorescence intensity) for each replicate array and fluorophore.**

Thresholding
------------

The first correction we perform is an area threshold to eliminate wells that didn't fully fill from further analysis.

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/area%20threshold-1.png"|absolute_url}})

**Figure 6. Well Area Threshold. A density curve of well area for each array of this chip is shown with the area threshold drawn in green at 27. All wells below the threshold are empty or partially filled. These individual wells will be marked as "Filled = FALSE" and will not be evaluated in results. This threshold can be adjusted if necessary in the code chunk above.**

We also noticed some sample evaporation occurred at the edges of each replicate. We can remove wells on the edges from further analysis.

Our dyes and filters are not perfectly separated in the visible spectrum, and some of our FAM fluorescence shows up in our HEX channel. The more FAM fluorescence in a well, the more artifact we see in the HEX channel for that well. We previously determined that a correction factor of 0.45 could account for this.

After filtering and correcting, we can finally draw thresholds to determine whether a well is positive vs. negative for the genetic sequence!

I tried a few methods for finding the negative peak. It is *usually* the largest peak in a density curve of fluorescence intensity, but sometimes the second mode of the bimodal distribution is larger. I've found the kmeans clustering, when setting the number of clusters as 2, works well to consistently identify the two populations.

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/calculate%20maximums%20via%20kmeans-1.png"|absolute_url}})

**Figure 7. Histogram of fluorescence AID for each channel and array, with color by cluster distinction**

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/view%20amp%20thresh%20without%20density%20curve-1.png"|absolute_url}})    

**Figure 8. Threshold for amplification fluorophore.**

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Preliminary%20View%20of%20wt%20and%20mut%20Thresholds-1.png"|absolute_url}})

**Figure 9. Thresholds for mutant and wild-type fluorophore. A scatter plot of mutant vs. wild-type AIDs for each filled well is shown with the generated mutant and wild-type thresholds. These quadrants will be used to determine zygosity in each well.**

Results
-------

Using the calculated thresholds, we can assign a zygosity to each well where a reaction occurred: wild-type, mutant, or heterozygous

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/View%20Results%20of%20wt%20and%20mut%20Thresholding-1.png"|absolute_url}})

**Figure 10. Zygosity of each well that passed quality control filtering. Color indicates zygosity.**

We can also take a peak at wells that did not get assigned a zygosity.  
![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/View%20Results%20Fluor%20QC-1.png"|absolute_url}})

**Figure 11. Condition of any wells that was filled but not assigned a zygosity. Since we required wells to have both amplification control positivity and one other positive reaction to be labeled with a zygosity, this removed some questionable wells.**

**Table showing well counts per array. Every filled well is either Amp pos or neg, FilledAmped = WT + MUT + HET + ampOnly; and FilledNoAmp = NOAMPwt + NOAMPmut + NOAMPboth + TripleNeg.**

| Array | Template |  TotalFilled|  FilledAmped|   WT|  MUT|  HET|  ampOnly|  FilledNoAmp|  NOAMPwt|  NOAMPmut|  NOAMPboth|  TripleNeg|
|:------|:---------|------------:|------------:|----:|----:|----:|--------:|------------:|--------:|---------:|----------:|----------:|
| A1    | OCI-AML3 |          763|          212|   75|   81|   56|        0|          551|       15|        15|          0|        521|
| A2    | OCI-AML3 |          558|          180|   74|   67|   38|        1|          378|        9|         3|          0|        366|
| A3    | OCI-AML3 |          428|          115|   39|   46|   30|        0|          313|        2|         4|          1|        306|

### Cell Specific Results

Now we can combine our zygosity measurements with our known cell locations, and find the zygosity of wells that contained single cells!  
![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Plot%20zygosities%20of%20cells-1.png"|absolute_url}})

**Figure 12. Bar plot of single-cell zygosities.**

We can also calculate some performance metrics for our device. For instance, we can calculate:
false positives: no cell was present but the reaction occurred false negative: a cell was present but no reaction did not occur true positive: a cell was present and the reaction occurred true negative: no cell was present and no reaction occurred
![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Failure%20Rates-1.png"|absolute_url}})

**Figure 13. Bar plot of device error rates.**

**Table showing error rates for each array.**
| Array | Template |  FalsePositives|  FalseNegatives|  TruePositives|  TrueNegatives|  DoubletPlus|  EmptyWells|        FPR|        FNR|
|:------|:---------|---------------:|---------------:|--------------:|--------------:|------------:|-----------:|----------:|----------:|
| A1    | OCI-AML3 |             118|              44|             71|            496|           34|         261|  0.1921824|  0.3826087|
| A2    | OCI-AML3 |             118|              29|             56|            347|            8|         466|  0.2537634|  0.3411765|
| A3    | OCI-AML3 |              77|              12|             33|            298|            8|         596|  0.2053333|  0.2666667|


We can also see if our distribution of cells in the replicates was random. If it were random, we would expect the counts of empty, single cell, and multiple-cell wells to fit the Poisson density function.  
![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Poisson%20Calculations-1.png"|absolute_url}})

**Figure 14. Actual counts of cell distributions versus those predicted by the Poisson density function.**

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Poisson%20Summary-1.png"|absolute_url}})

**Figure 15. For sum of all replicates, actual counts of cell distributions versus those predicted by the Poisson density function.**

### Array Maps

Finally, we can veiw a representative map of the original images, to find locations of the zygosities and cells.  
![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Zygosity%20Array%20Maps-1.png"|absolute_url}})

**Figure 16. Map of zygosity locations.**

![]({{"/images/GenotypingDeviceAnalysis_files/figure-markdown_github/Cell%20Array%20Maps-1.png"|absolute_url}})

**Figure 17. Map of single cell locations.**
