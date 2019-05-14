# SNPs 'n Yarn
Generative art image combining a correlation plot and a Manhattan plot

## Description

This script generates two plots, a correlation matrix with the nodes positioned on a circle and the correlation between each node as a curve.
Strength of the correlation is indicated by a color on a color scale. The Manhattan plot is plotted using polar coordinates, with the chromosomes going clockwise.

## Required packages

This script requires three packages, two of which are as of now not publicly available. The `tidyverse` package is used to modulate and plot the data using the `ggplot` function.
The two packages that are not yet publicly available are used to generate sample GWAS data (`normentR`( and to set the color scales (`danielR`).
The color settings can be easily circumvented. The scales used in the example images are all from the `scico` package, using the `vik` palette. The `danielR` functions are nothing more than a wrapper for some of these color settings.
To generate sample GWAS data, a function we wrote for the `normentR` package is used (`simulateGWAS`). This function outputs a data frame that mimics the output from the `--qassoc` flag in PLINK.
If you have such a file available, then use that, otherwise there are probably other packages available to simulate GWAS data, or it's possible to simulate some simple data.
The only data used to plot the GWAS data are in the `CHR`, `BPcum`, and `P` columns.

## Creation of the sample images

The sample images as included in this repository can be created by overlaying the transparent plots on a grey surface and scaling them to make them fit the required format (50x50cm 300dpi canvas).
The `ggsave` function is optional. I used the saved .png's to generate the final image that was used for submission in Photoshop.

## Development

This code will be updated from time to time when I learn more about the possibilities of the `tidyverse` package, and to filter out bugs.
Thanks goes to my colleagues for feedback and input on the whole process of writing this script and visualization of the data!
