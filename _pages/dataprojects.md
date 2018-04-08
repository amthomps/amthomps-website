---
title: "Alison Thompson - Data Projects"
layout: textlay
excerpt: "Data Projects"
sitemap: false
permalink: /dataprojects.html
---

## Data projects  

Some of my projects for data analysis in R and ImageJ are presented here.  

### Genotyping Device Analysis  

I recently submitted a manuscript describing a device I developed to perform single-cell genotyping of single cells.  As part of this project, I built a script to produce reproducible summary reports for every sample I run.  These reports include summary statistics, error rates, and descriptive visualizations.  These reports are built in R using the knitr package.  

[View the report version with code.](GenotypingDeviceAnalysis.md)   

[View the report version without code.](GenotypingDeviceReport.md)   

### Particle Analysis in ImageJ  

[ImageJ](https://imagej.nih.gov/ij/) is an open source image processing program.  Custom analysis pipelines can be written in ImageJ's macro language.  In my research, one of the most useful built-in ImageJ functions has been *Analyze Particles*.  In this example, I use a custom macro to locate single cells in an image captured using a fluorescence microscope.  Using the *Analyze Particles* function, I can find the size and fluorescence intensity of each cell in the image.  I can then summarize data from multiple experimental conditions in R to find optimal conditions for cell staining.

[View the Analyze Particles walkthrough.](AnalyzeParticles)
