sacc-tDCS: Frontal Eye Field coordinates
================
Leon Reteig

-   [Setup](#setup)
-   [Coordinates in native space](#coordinates-in-native-space)
-   [Transform to MNI space](#transform-to-mni-space)
    -   [Run shell script](#run-shell-script)
    -   [MNI coordinates](#mni-coordinates)
-   [Create spherical ROI for each coordinate](#create-spherical-roi-for-each-coordinate)
    -   [Exclude subjects](#exclude-subjects)
    -   [Run script](#run-script)

R notebook for analysis of frontal eye field coordinates in the `sacc-tDCS` dataset.

Setup
=====

``` r
# Load some libraries
library(here) # file paths
```

    ## here() starts at /Volumes/research$/reteig/sacc-tDCS

``` r
library(tidyverse) # importing, transforming, and visualizing data frames
```

    ## Loading tidyverse: ggplot2
    ## Loading tidyverse: tibble
    ## Loading tidyverse: tidyr
    ## Loading tidyverse: readr
    ## Loading tidyverse: purrr
    ## Loading tidyverse: dplyr

    ## Conflicts with tidy packages ----------------------------------------------

    ## filter(): dplyr, stats
    ## lag():    dplyr, stats

``` r
library(knitr) # R markdown output (html, pdf, etc.)
# set default output and figure options
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 7, fig.asp = 0.618, out.width = "75%", fig.align = "center")

sessionInfo()
```

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: OS X El Capitan 10.11.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] C
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] knitr_1.15.1    dplyr_0.5.0     purrr_0.2.2     readr_1.1.0    
    ## [5] tidyr_0.6.1     tibble_1.3.0    ggplot2_2.2.1   tidyverse_1.1.1
    ## [9] here_0.1       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.10     cellranger_1.1.0 compiler_3.4.0   plyr_1.8.4      
    ##  [5] forcats_0.2.0    tools_3.4.0      digest_0.6.12    lubridate_1.6.0 
    ##  [9] jsonlite_1.4     evaluate_0.10    nlme_3.1-131     gtable_0.2.0    
    ## [13] lattice_0.20-35  psych_1.7.3.21   DBI_0.6-1        yaml_2.1.14     
    ## [17] parallel_3.4.0   haven_1.0.0      xml2_1.1.1       stringr_1.2.0   
    ## [21] httr_1.2.1       hms_0.3          rprojroot_1.2    grid_3.4.0      
    ## [25] R6_2.2.0         readxl_1.0.0     foreign_0.8-67   rmarkdown_1.5   
    ## [29] modelr_0.1.0     reshape2_1.4.2   magrittr_1.5     backports_1.0.5 
    ## [33] scales_0.4.1     htmltools_0.3.6  rvest_0.3.2      assertthat_0.2.0
    ## [37] mnormt_1.5-5     colorspace_1.3-2 stringi_1.1.5    lazyeval_0.2.0  
    ## [41] munsell_0.4.3    broom_0.4.2

List the subjects that were excluded in the eye data analyses.

``` r
subs2exclude <- c("S21","S25","S16","S22","S28")
```

Coordinates in native space
===========================

These were determined for each subject's scan; see `neuronav_notes.md` for further info.

``` r
dataFile <- here("data", "FEF_coords_native.csv")
nativeCoords <- read_csv2(dataFile)
nativeCoords %>% 
  select(-folder, -scan) %>% # drop columns with folder and scan names
  filter(!(subject %in% subs2exclude)) %>% # drop rows with excluded subjects
  kable(.)
```

| subject |  native.X|  native.Y|  native.Z|
|:--------|---------:|---------:|---------:|
| S01     |        96|       132|       137|
| S02     |        89|       118|       132|
| S03     |       153|       124|       177|
| S04     |        96|       127|       131|
| S05     |        93|       116|       138|
| S06     |       101|       131|       182|
| S07     |        92|       127|       147|
| S08     |       101|       123|       159|
| S09     |        99|       125|       158|
| S10     |        94|       124|       125|
| S11     |       147|       123|       149|
| S12     |        97|       121|       176|
| S13     |       144|       125|       140|
| S14     |        87|       122|       129|
| S15     |       118|       111|       125|
| S17     |       161|       135|       175|
| S18     |       153|       147|       173|
| S19     |       146|       155|       172|
| S20     |        78|       124|        99|
| S24     |       160|       134|       141|
| S26     |       156|       141|       150|
| S27     |       148|       122|       156|
| S29     |        90|       131|       153|
| S30     |        95|       125|       149|
| S32     |       147|       140|       146|
| S33     |       159|       149|       168|

Transform to MNI space
======================

Run shell script
----------------

`mni_coords.sh` reads in the native coordinates, converts them to MNI space, and writes out a `.csv` of the same format as for the native coordinates.

Under `# FSL setup`, replace `FSLDIR` with the directory path you installed FSL to.

Arguments:

1.  **Input:** csv file with coordinates in native space
2.  **Output:** csv file coordinates in MNI space
3.  Passing `all` as a 3rd argument first (re)runs skull stripping and registration to the MNI template.

``` bash
# FSL Setup
FSLDIR=/usr/local/fsl
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
. ${FSLDIR}/etc/fslconf/fsl.sh

# Run shell script to compute MNI coordinates
bash src/func/mni_coords.sh data/FEF_coords_native.csv data/FEF_coords_MNI.csv all
```

MNI coordinates
---------------

Load in the MNI coordinates that were written to a `.csv` file by the shell script.

``` r
dataFile <- here("data", "FEF_coords_MNI.csv")
mniCoords <- read_delim(dataFile, ";")
mniCoords <- filter(mniCoords, !(subject %in% subs2exclude)) # exclude subjects
mniCoords %>%
  select(-folder, -scan) %>% # drop columns with folder and scan names
  kable(.)
```

| subject |   MNI\_X|      MNI\_Y|   MNI\_Z|
|:--------|--------:|-----------:|--------:|
| S01     |  29.4081|    1.069700|  54.8724|
| S02     |  33.0338|   -2.245090|  50.4239|
| S03     |  30.5837|   -1.480360|  50.5881|
| S04     |  25.7061|   -3.761730|  56.3648|
| S05     |  29.7780|   -5.201110|  55.8046|
| S06     |  29.7783|   -1.120820|  58.2622|
| S07     |  38.1233|    2.975730|  45.9613|
| S08     |  31.5014|    0.526375|  45.6200|
| S09     |  28.5103|    3.632940|  51.2804|
| S10     |  28.1080|   -1.933630|  50.7210|
| S11     |  30.5811|   -3.787490|  51.9535|
| S12     |  36.5129|   -0.386049|  46.7630|
| S13     |  26.2263|   -1.069220|  54.7079|
| S14     |  37.5005|   -1.588150|  52.5889|
| S15     |  31.8188|   -8.357000|  58.9671|
| S17     |  31.0229|   -5.121250|  54.3116|
| S18     |  34.9669|    8.390570|  49.8188|
| S19     |  28.0703|   -3.834910|  52.7607|
| S20     |  41.2127|   -1.745790|  47.5514|
| S24     |  37.3474|   -0.862753|  43.4371|
| S26     |  34.3288|   -2.861210|  49.1761|
| S27     |  27.6566|  -10.061100|  50.9743|
| S29     |  30.3161|   -5.272910|  55.2661|
| S30     |  26.8320|   -3.868670|  54.5619|
| S32     |  29.0049|    4.888810|  49.1404|
| S33     |  30.3099|   -3.897600|  50.9083|

Calculate statistics over subjects:

``` r
mniCoords %>%
  gather(dimension, coord, MNI_X:MNI_Z) %>%
  group_by(dimension) %>%
  summarise(average = mean(coord), standard.deviation = sd(coord), minimum = min(coord), maximum = max(coord)) %>%
  kable(.)
```

| dimension |    average|  standard.deviation|   minimum|   maximum|
|:----------|----------:|-------------------:|---------:|---------:|
| MNI\_X    |  31.470735|            4.036810|   25.7061|  41.21270|
| MNI\_Y    |  -1.806643|            3.931623|  -10.0611|   8.39057|
| MNI\_Z    |  51.645608|            3.922225|   43.4371|  58.96710|

Create spherical ROI for each coordinate
========================================

Exclude subjects
----------------

Write a new \`.csv' file with MNI coordinates, without the coordinates from the excluded subjects.

``` r
exclMNI <- here("data", "FEF_coords_MNI_excl.csv")
write_delim(mniCoords, exclMNI, delim = ";")
Sys.setenv(exclMNI = exclMNI) # port variable with csv file name to bash
```

Run script
----------

The coords can be very nicely visualized with [Surf Ice](https://www.nitrc.org/projects/surfice/) by running the `plot_coords.m` script (requires MATLAB)

Alternatively, `fef_rois.sh` takes the MNI coordinates and adds spherical masks, centered around each subjects FEF coordinate, to the MNI template.

Under `# FSL setup`, replace `FSLDIR` with the directory path you installed FSL to.

Arguments:

1.  **Input:** csv file with coordinates in MNI space
2.  **Output:** NIfTI file of `MNI_152_1mm` template with spherical ROI around each subject's frontal eye field
3.  Radius of the spherical ROI in mm

``` bash

# FSL Setup
FSLDIR=/usr/local/fsl
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
. ${FSLDIR}/etc/fslconf/fsl.sh

# Run shell script to make spherical ROIs
bash src/func/fef_rois.sh $exclMNI neuronav/FEF_ROIs_bigspheres.nii.gz 2
```
