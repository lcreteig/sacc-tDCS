sacc-tDCS: Data inspection
================
Leon Reteig

-   [Load data](#load-data)
-   [Inspect distributions](#inspect-distributions)
    -   [Histograms for each subject](#histograms-for-each-subject)
    -   [Stimulation effects across subjects](#stimulation-effects-across-subjects)
    -   [Session effects in each subject](#session-effects-in-each-subject)
-   [Outliers](#outliers)
    -   [Outlier trials](#outlier-trials)
-   [Outliers in median latency](#outliers-in-median-latency)
-   [Drift correction](#drift-correction)

R notebook for inspection of eye tracking data in the `sacc-tDCS` dataset. Previous processing:

-   Raw data were parsed into events (saccades, fixations, etc.) by the EyeLink data were collected on.
-   Events were extracted and saccade measures were computed with a MATLAB script.

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
library(formatR)
# set default output and figure options
knitr::opts_chunk$set(message = FALSE, warning = FALSE, tidy.opts=list(width.cutoff=80), tidy=TRUE, fig.width = 7, fig.asp = 0.618, out.width = "75%", fig.align = "center")

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
    ##  [1] formatR_1.5     knitr_1.15.1    dplyr_0.5.0     purrr_0.2.2    
    ##  [5] readr_1.1.0     tidyr_0.6.1     tibble_1.3.0    ggplot2_2.2.1  
    ##  [9] tidyverse_1.1.1 here_0.1       
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

Load data
=========

The .csv file with the eye tracking data was created in MATLAB.

``` r
# Load the data frame
dataFile <- here("data", "sacc-tDCS_data.csv")
groupData <- read_csv(dataFile, col_names = TRUE, na = "NaN", progress = FALSE, col_types = cols(stimulation = col_factor(c("anodal", 
    "cathodal")), leg = col_factor(c("pre", "tDCS", "post")), type = col_factor(c("lateral", 
    "center")), direction = col_factor(c("left", "right"))))
```

``` r
kable(head(groupData))
```

| subject | stimulation | leg |  block|  trial| type    | direction |  deviation.start|  deviation.end.x|  deviation.end.y|  amplitude|  latency|    drift.x|   drift.y|
|:--------|:------------|:----|------:|------:|:--------|:----------|----------------:|----------------:|----------------:|----------:|--------:|----------:|---------:|
| S01     | anodal      | pre |      1|      1| lateral | right     |         0.462897|         0.170455|       -0.0080638|    8.02463|      433|  0.0953736|  0.102814|
| S01     | anodal      | pre |      1|      1| center  | left      |         0.459092|         1.032850|        0.0665268|    7.16262|      439|  0.0953736|  0.102814|
| S01     | anodal      | pre |      1|      2| lateral | right     |         0.344561|        -0.344967|        0.2197400|    7.74873|      291|  0.0953736|  0.102814|
| S01     | anodal      | pre |      1|      2| center  | left      |         0.550230|         0.361201|        0.3507760|    7.43233|      198|  0.0953736|  0.102814|
| S01     | anodal      | pre |      1|      3| lateral | right     |         0.514736|        -0.588470|        0.1673250|    7.61080|      281|  0.0953736|  0.102814|
| S01     | anodal      | pre |      1|      3| center  | left      |         0.620728|         1.576610|        0.4031910|    6.35043|      376|  0.0953736|  0.102814|

-   **subject**: subject ID
-   **stimulation**: Whether data are from the `anodal` or `cathodal` session
-   **leg**: Whether data are before (`pre`), during (`tDCS`), or after (`post`) tDCS
-   **block**: After each block participant had a brief break and tracker was recalibrated
-   **trial**: trial number within a block
-   **type**:
    -   `lateral` - fixation in center of display, saccade made towards the periphery
    -   `center` - fixation in periphery, saccade made back towards the center of the display
-   **direction**: `left` for saccades towards the left of current fixation position; `right` for saccades to the right
-   **deviation.start** : distance (in visual angle) from saccade start point to fixation
-   **deviation.end.x**: distance (in visual angle) from x-coordinate of saccade end point to x-coordinate of target location
-   **deviation.end.y**: same for y-coordinate
-   **amplitude**: distance (in visual angle) between saccade start and end point
-   **latency**: time (in ms) from target onset to start of saccade
-   **drift.x**: distance (in visual angle) between x-coordinate of average fixation position during the break to x-coordinate of fixation stimulus. This stimulus was displayed at each break in the task, so this data can be used as an estimate of offsets to do drift correction.
-   **drift.y**: same for y-coordinate

Inspect distributions
=====================

### Histograms for each subject

``` r
histType <- ggplot(groupData, aes(latency, fill = type)) + facet_wrap(~subject, ncol = 5, 
    scales = "free_y") + geom_histogram(binwidth = 5, color = "grey50", size = 0.2) + 
    xlim(-50, 300)
histType
```

<img src="inspection_files/figure-markdown_github/Histogram per subject-1.png" width="75%" style="display: block; margin: auto;" />

**Stray observations:**

-   Center saccades are much faster than lateral, though not to the same degree in all subjects
-   Some have a fat short latency tail - too fast saccades that are virtually all towards the center: S10, S11, S22, S28
-   Some appear bimodal (S28), but when split for type this is generally because center saccades are faster: S05, S06, S11
-   Some look almost normally distributed (S10); others are very heavily right-skewed (S32)
-   Some are super sharp (S08); others really broad (S01)

### Stimulation effects across subjects

``` r
dens <- ggplot(groupData, aes(latency, color = stimulation, linetype = leg)) + facet_grid(type ~ 
    direction) + geom_density() + xlim(0, 250) + scale_color_brewer(palette = "Set1")
dens
```

<img src="inspection_files/figure-markdown_github/Density: stimulation across subjects-1.png" width="75%" style="display: block; margin: auto;" />

### Session effects in each subject

``` r
denstDCS <- ggplot(groupData[groupData$leg == "pre" & groupData$type == "lateral", 
    ], aes(latency, color = stimulation)) + facet_wrap(~subject, ncol = 5, scales = "free_y") + 
    geom_density() + xlim(0, 250) + scale_color_brewer(palette = "Set1") + ggtitle("Lateral saccades, baseline block")
denstDCS
```

<img src="inspection_files/figure-markdown_github/Density: anodal vs. cathodal baseline per subject-1.png" width="75%" style="display: block; margin: auto;" />

For most subjects, the latency distributions in both sessions are reasonably similar. Note that this is the baseline block, so we also wouldn't expect any differences. In that light, it is a little worrisome that for some subjects the distributions differ markedly (S01, S02, S09, S12, S17, S21, S26, S29, S32).

Outliers
========

Outlier trials
--------------

``` r
tooFast <- 50
tooSlow <- 400
badFix <- 1.8
badSacc <- 8
```

Criteria for outliers:

-   Discard fast saccades, with a latency of 50 ms or less
-   Discard slow saccades, saccades with a latency of 400 ms or more
-   Discard inaccurate fixations, with saccade starting point more than 1.8 degrees or more away from fixation
-   Discard faulty saccades, with x-coordinate of saccade end point 8 degree or more away from the target

In [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045), this was:

-   Fast saccades: 50 ms
-   Slow saccades: 400 ms
-   Bad fixations: 1.8 degrees
-   Faulty saccades: opposite hemifield of target (here, that would be 8 degrees as targets were that eccentric)

``` r
# Mark outliers
groupData <- mutate(groupData, outlier = "non.outlier", # fill vector for all trials
                    outlier = ifelse(latency < tooFast, "fast", outlier), # mark too fast trials as "fast"
                    outlier = ifelse(latency > tooSlow, "slow", outlier), # mark too slow trials as "slow"
                    outlier = ifelse(deviation.start > badFix, "fixation", outlier), # mark bad fixations as "fixation"
                    outlier = ifelse(deviation.end.x > badSacc, "saccade", outlier), # mark inaccurate saccades as "saccade"
                    outlier = ifelse(is.na(latency), "none", outlier) # mark absence of saccade as "none"
                    )
```

### Plot latencies of outlier trials per subject

``` r
outlierPlotTrials <- ggplot(groupData[!(groupData$outlier %in% c("none", "non.outlier")), ], aes(interaction(stimulation,leg,block,trial), latency, color = outlier, shape = type)) +
  facet_wrap(~subject, ncol = 5, scales = "free_y") +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = tooFast, linetype = "dashed") +
  geom_hline(yintercept = tooSlow, linetype = "dashed") +
  xlab('Trial') +
  theme(axis.text.x = element_blank(), # remove x-axis (just trial count)
  axis.ticks.x = element_blank())
outlierPlotTrials
```

<img src="inspection_files/figure-markdown_github/Plot outlier trials per subject-1.png" width="75%" style="display: block; margin: auto;" />

### Tally outlier types

``` r
outlierCount <- groupData %>% group_by(subject, stimulation, leg, direction, type, 
    outlier) %>% summarize(outlier_count = n())  # for each condition and subject, count how many (non)outliers there are
```

``` r
outlierTable <- outlierCount %>%
  group_by(subject,outlier) %>%
  summarize(outlier_count = sum(outlier_count)) %>% # create column with sum across all conditions per subject
  mutate(total = sum(outlier_count[outlier != "non.outlier"])) %>% # create column with sum across all outlier types
  spread(outlier, outlier_count) %>% # one column per outlier type
  select(subject, non.outlier, total, none, fast, slow, fixation, saccade) #reorder columns
kable(outlierTable, caption = "Number of outlier saccades per subject")
```

| subject |  non.outlier|  total|  none|  fast|  slow|  fixation|  saccade|
|:--------|------------:|------:|-----:|-----:|-----:|---------:|--------:|
| S01     |         5414|    346|   166|    49|    70|        57|        4|
| S02     |         5711|     49|    14|    13|    NA|        21|        1|
| S03     |         5389|    371|    89|    89|     4|       181|        8|
| S04     |         5563|    197|    44|    47|     2|        96|        8|
| S05     |         5009|    751|   215|   211|     9|       296|       20|
| S06     |         5075|    685|   255|   213|    NA|       202|       15|
| S07     |         5641|    119|    71|    31|    NA|        13|        4|
| S08     |         5367|    393|   214|   135|     1|        36|        7|
| S09     |         5586|    174|    59|    54|     9|        49|        3|
| S10     |         4578|   1182|   454|   352|     2|       363|       11|
| S11     |         5177|    583|   280|   164|    13|       115|       11|
| S12     |         5593|    167|    66|    54|     2|        40|        5|
| S13     |         5191|    569|   394|    87|     2|        80|        6|
| S14     |         5036|    724|   299|   119|     6|       279|       21|
| S15     |         5279|    481|   329|    75|     3|        56|       18|
| S16     |         3879|   1881|   858|   432|    10|       512|       69|
| S17     |         5210|    550|   160|   140|     5|       233|       12|
| S18     |         5717|     43|    15|    16|    NA|         8|        4|
| S19     |         5357|    403|   114|   141|     1|       141|        6|
| S20     |         5110|    650|   272|   188|     9|       162|       19|
| S21     |         5467|    293|    66|    41|    13|       168|        5|
| S22     |         2898|   2862|   847|   297|     7|      1707|        4|
| S24     |         4568|   1192|   571|   179|     4|       422|       16|
| S25     |         4185|   1575|   514|   384|    NA|       647|       30|
| S26     |         5628|    132|    30|    18|     5|        77|        2|
| S27     |         5624|    136|    26|    41|    NA|        63|        6|
| S28     |         2781|   2979|  1913|   329|     4|       694|       39|
| S29     |         5223|    537|   148|   119|    12|       245|       13|
| S30     |         5094|    666|   297|   230|     2|       137|       NA|
| S32     |         5078|    682|   106|   116|     5|       444|       11|
| S33     |         5387|    373|   200|   108|     8|        55|        2|

### Plot outlier counts per subject and type

``` r
max_n <- nrow(filter(groupData, subject == "S01", type == "center"))  # max amount of saccades in experiment

ggplot(filter(outlierCount, outlier != "non.outlier"), aes(subject, outlier_count, 
    fill = outlier)) + geom_col() + scale_y_continuous("number of saccades", limits = c(0, 
    max_n), sec.axis = sec_axis(~./max_n * 100, name = "percent of all saccades")) + 
    coord_flip() + facet_wrap(~type)
```

<img src="inspection_files/figure-markdown_github/Plot outlier counts per subject and type-1.png" width="75%" style="display: block; margin: auto;" />

**Stray observations:**

Differences between subjects:

-   Some subjects have extremely clen data with barely any outliers (S02, S07, S18)
-   Most subjects have quite a few outliers, especially S10, S16, S22, S25 and S28
-   Only S01 has a sizable amount of slow saccades
-   Only S16 has a sizable amount of saccades to the opposite hemifield
-   Those subjects with many fast saccades also tend to have many bad fixations

General patterns:

-   Occurence of outliers seems stable throughout the session: there aren't more/less in the beginning / end
-   There are very few inaccurate saccades (makes sense, because task is easy and criterion is not strict)
-   Most outliers are too fast saccades and bad fixations
-   Slow saccades are lateral, fast saccades are to the center, because only the latter are predictable
-   Bad fixations appear to be mostly center saccades (but this varies a lot). Perhaps the eyes already drift back towards the center, before executing the saccade?
-   Or is the source of bad fixations simply poor quality of the eye tracker data? e.g. people are actually fixating, but due to drift it appears they are not.

### Tally number of non-outlier saccades

Importantly, we should also see how many saccades are left per condition after excluding outliers

``` r
trialCount <- groupData %>%
  filter(outlier == "non.outlier") %>% # keep only non-outlier saccades
  # Equally cut leg factor into 15-minute intervals, because the "post" leg is currently 30 minutes 
  mutate(leg = as.character(leg), # cannot edit leg if it's still a factor
         leg = replace(leg, leg == "post" & block <= 3, "post.1"),
         leg = replace(leg, block > 3, "post.2"),
         leg = factor(leg, levels = c("pre", "tDCS", "post.1", "post.2")) # refactor and order levels
         ) %>%
  group_by(subject,stimulation,leg,direction,type) %>%
  summarise(saccades = n()) %>% # count how many there are per subject per condition
  arrange(saccades) # sort ascending
kable(head(trialCount))
```

| subject | stimulation | leg    | direction | type   |  saccades|
|:--------|:------------|:-------|:----------|:-------|---------:|
| S28     | anodal      | post.1 | right     | center |         2|
| S28     | anodal      | post.2 | right     | center |         2|
| S28     | cathodal    | post.2 | right     | center |         3|
| S28     | anodal      | post.1 | left      | center |         5|
| S28     | anodal      | post.2 | left      | center |         6|
| S28     | cathodal    | post.1 | right     | center |         6|

### List candidates for rejection

S28 should definitely be rejected, as certain conditions have only 2 useable saccades! S16 and S22 (together with S28) also have quite few saccade counts: less than 50 in some conditions. This is mostly because there are many missing saccades (`none`). Inspection of the data shows that this is not due to poor data quality, but because these subjects move their eyes too soon (i.e. before the stimulus even). That also explains why these low saccade counts occur almost exclusively in the *center* saccade condition (as there the location of the target was predictable).

``` r
subs2exclude <- c("S28", "S16", "S22", "S21", "S25")
```

S21 and S25 should also be excluded. Their anodal and cathodal sessions were separated by less than 48 hours, which is in violation of the protocol.

### Saccade counts for included participants

For all participants that will be included in the data analysis, compute descriptives of how many valid saccades remain for each type per cell (i.e. each stimulation, direction and leg combination used for statistical analysis)

``` r
trialCount %>%
  filter(!(subject %in% subs2exclude)) %>% # for all included participants
  group_by(type) %>%
  summarise(average = mean(saccades), standard.deviation = sd(saccades), minimum = min(saccades), maximum = max(saccades)) %>%
  kable(.)
```

| type    |   average|  standard.deviation|  minimum|  maximum|
|:--------|---------:|-------------------:|--------:|--------:|
| lateral |  175.2284|            6.244285|      142|      180|
| center  |  155.5529|           20.816516|       74|      180|

Let's also look at proportion of saccades for each outlier type:

``` r
n_total <- nrow(filter(groupData, !(subject %in% subs2exclude)))  # total amount of saccades across all included sessions/subjects
outlierCount %>% filter(!(subject %in% subs2exclude), !(outlier %in% c("non.outlier", 
    "none"))) %>% group_by(outlier) %>% summarise(percentage = sum(outlier_count)/n_total * 
    100) %>% kable(.)
```

| outlier  |  percentage|
|:---------|-----------:|
| fast     |   1.9958600|
| fixation |   2.5848024|
| saccade  |   0.1555823|
| slow     |   0.1161859|

Outliers in median latency
==========================

Next to rejecting outlier trials, we could also consider rejecting outlier subjects or certain conditions from the statistical tests. One way to detect outliers (that is itself robust to outliers, unlike the standard deviation) is the MAD-median rule (see [this blogpost](https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/) and [this preprint](http://dx.doi.org/10.1101/151811)). The MAD is the *median absolute devian from the median*.

``` r
outliersPerCondition <- groupData %>%
  filter(outlier == "non.outlier") %>% # drop all outlier trials
  group_by(subject,stimulation,leg,type,direction) %>%
  summarise(latency = median(latency)) %>% # compute median latency
  group_by(stimulation,leg,type,direction) %>%
  mutate(mad.median.rule = (abs(latency - median(latency)) / mad(latency))) %>% # deviation from the median, standardized by the MAD
  filter(mad.median.rule > 2.24) # MAD-median rule
kable(head(outliersPerCondition))
```

| subject | stimulation | leg  | type    | direction |  latency|  mad.median.rule|
|:--------|:------------|:-----|:--------|:----------|--------:|----------------:|
| S01     | anodal      | pre  | lateral | left      |    219.5|         4.890058|
| S01     | anodal      | pre  | lateral | right     |    195.0|         2.810378|
| S01     | anodal      | pre  | center  | left      |    214.0|         7.194568|
| S01     | anodal      | pre  | center  | right     |    211.0|         3.943177|
| S01     | anodal      | tDCS | lateral | left      |    210.5|         4.143300|
| S01     | anodal      | tDCS | lateral | right     |    192.0|         3.035208|

These are all the subject-condition combinations for which the MAD-median rule is violated. If we were to reject subjects with one more more violation, we would have to remove 11 subjects.

There is some overlap with the analysis of outlier trials. For instance, there is reason to remove S28 in both analyses. But this is not always the case: S26 has very clean trial-data, but is still flagged as an outlier here because their median latencies are quite slow.

Drift correction
================

Calibration isn't perfect, so there are always small offsets between the measurements and what people are actually looking at. Further, these measurement errors can increase with time away from calibration, which is known as drift.

After every 20 trials (40 saccades) there was a break in the task, in which we asked subjects to fixate the center of the screen before continueing. The offsets recorded here should thus be a good estimate of drift, since here you can trust that subjects were actually looking at fixation spot on.

To do drift correction, we simply subtract the offsets recorded in the break from the x- and y- coordinates of the eye data of interest.

``` r
# Add columns with drift-corrected values
groupData <- mutate(groupData, deviation.end.x.corr = deviation.end.x - drift.x, 
    deviation.end.y.corr = deviation.end.y - drift.y)
```

If this does indeed work, then most trials should now have a smaller deviation.

``` r
# Calculate percentage of trials with smaller deviations after drift correction
groupData %>%
  group_by(subject, stimulation) %>% # split per subject and session
  summarize(
    trials = n(),
    smaller.drift.x = sum(abs(deviation.end.x.corr) < abs(deviation.end.x), na.rm = TRUE) / trials * 100,
    smaller.drift.y = sum(abs(deviation.end.y.corr) < abs(deviation.end.y), na.rm = TRUE) / trials * 100
  ) %>%
  
  # Plot
  ggplot(aes(smaller.drift.x, smaller.drift.y)) +
  facet_wrap(~stimulation) +
  geom_hline(yintercept = 50, linetype = "dashed") + geom_vline(xintercept = 50, linetype = "dashed") +
  geom_point(aes(color = subject)) +
  coord_fixed(xlim = c(0,100), ylim = c(0,100)) +
  #geom_text_repel(aes(label = subject, color = subject)) +
  theme(legend.position = "none")
```

<img src="inspection_files/figure-markdown_github/Check correction-1.png" width="75%" style="display: block; margin: auto;" /> Seems that it doesn't really work at all, because for most subjects the errors are decreased on less than 50% of saccades! Of course, this is only a rough analysis, but this does match with SR Research's advice in the EyeLink manual, which states that drift correction may actually deteriorate the calibration maps.
