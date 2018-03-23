## Welcome to the research website for Alison Thompson!

## Introduction  
This website was created to present some of my research and projects.  It is currently a work in progress, but I hope to continue to add content as I expand my data science portfolio and learn more about web publishing.

This website was created with Jekyll and was modified from a website of another GitHub user.  If you are building a website and like this style, I recommend that you check out this [GitHub repo](https://github.com/allanlab/allanlab) and read this [about page](http://www.allanlab.org/aboutwebsite.html).  These were a huge help to me.  

## Instructions  
#### To Publish a New page
New content written in markdown should be added to the \_pages folder.  A corresponding file with instructions on how to convert the markdown to html should be added to the \_layouts folder.  All changes should be pushed to the master branch.

#### To Edit Tabs in the Header
To add or modify tabs, navigate to \_includes/header.html and edit the navbar lines.  Push to the master.

#### To Update CV  
Replace CV.pdf in the \_pages folder.  If adding or removing pages, adjust the iframe size in \_pages/cv.md. Push to the master.

#### To Add Projects Produced in RStudio
This process is still a bit clunky, but until I find a better way, follow these instructions.
1. Use knitr to knit the document to a github_document
2. Move the .md document to the \_pages folder
3. Move the folder containing plots and images to the images folder
4. Update the links in the .md document to link to the corresponding images
5. Add html page in the \_layouts
6. Push to the master.
