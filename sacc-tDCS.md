No evidence that frontal eye field tDCS affects latency or accuracy of prosaccades
================
Leon Reteig, Tomas Knapen, Richard Ridderinkhof, Heleen Slagter

-   [Group analyses](#group-analyses)
    -   [Load data](#load-data)
        -   [Eye tracking data](#eye-tracking-data)
        -   [Metadata](#metadata)
    -   [Preprocess eye data](#preprocess-eye-data)
        -   [Reject outlier saccades](#reject-outlier-saccades)
        -   [Exclude participants](#exclude-participants)
        -   [Collapse data across 15 minute intervals](#collapse-data-across-15-minute-intervals)
    -   [Subject and session descriptives](#subject-and-session-descriptives)
        -   [Sample demographics](#sample-demographics)
        -   [tDCS session order](#tdcs-session-order)
        -   [Saccade counts](#saccade-counts)
    -   [Saccade latency](#saccade-latency)
        -   [Figure 4](#figure-4)
        -   [Statistics](#statistics)
    -   [Saccade endpoint deviation](#saccade-endpoint-deviation)
        -   [Figure 6](#figure-6)
        -   [Statistics](#statistics-1)
    -   [Saccade endpoint variability](#saccade-endpoint-variability)
        -   [Figure 7](#figure-7)
        -   [Statistics](#statistics-2)
-   [Quantile analysis](#quantile-analysis)
    -   [Figure 5](#figure-5)
-   [Supplementary results](#supplementary-results)
    -   [tDCS adverse effect questionnaire](#tdcs-adverse-effect-questionnaire)
        -   [Figure S1](#figure-s1)
        -   [Statistics](#statistics-3)
    -   [Frontal eye field coordinates](#frontal-eye-field-coordinates)
        -   [Native space](#native-space)
        -   [MNI space (Table S1)](#mni-space-table-s1)

``` r
# Load some libraries
library(here) # file paths
```

    ## here() starts at /Volumes/psychology$/Researchers/reteig/sacc-tDCS

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
library(stringr) # manipulating strings
library(broom) # # transform model output into a data frame
library(cowplot) # formatting plots
```

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

``` r
library(knitr) # R notebook output
library(ez) # classical ANOVA
library(BayesFactor) # Bayesian statistics
```

    ## Loading required package: coda

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## ************
    ## Welcome to BayesFactor 0.9.12-2. If you have questions, please contact Richard Morey (richarddmorey@gmail.com).
    ## 
    ## Type BFManual() to open the manual.
    ## ************

``` r
library(rogme) # shift functions

# set default output and figure options
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 7, fig.asp = 0.618, out.width = "75%", fig.align = "center")

source(here("src", "lib", "InclusionBF.R"))
print(sessionInfo())
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
    ##  [1] rogme_0.1.0.9000     BayesFactor_0.9.12-2 Matrix_1.2-9        
    ##  [4] coda_0.19-1          ez_4.4-0             knitr_1.15.1        
    ##  [7] cowplot_0.7.0        broom_0.4.2          stringr_1.2.0       
    ## [10] dplyr_0.5.0          purrr_0.2.2          readr_1.1.0         
    ## [13] tidyr_0.6.1          tibble_1.3.0         ggplot2_2.2.1       
    ## [16] tidyverse_1.1.1      here_0.1            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtools_3.5.0       pbapply_1.3-2      reshape2_1.4.2    
    ##  [4] splines_3.4.0      haven_1.0.0        lattice_0.20-35   
    ##  [7] colorspace_1.3-2   htmltools_0.3.6    yaml_2.1.14       
    ## [10] mgcv_1.8-17        nloptr_1.0.4       foreign_0.8-67    
    ## [13] DBI_0.6-1          modelr_0.1.0       readxl_1.0.0      
    ## [16] plyr_1.8.4         MatrixModels_0.4-1 munsell_0.4.3     
    ## [19] gtable_0.2.0       cellranger_1.1.0   rvest_0.3.2       
    ## [22] mvtnorm_1.0-6      psych_1.7.3.21     evaluate_0.10     
    ## [25] forcats_0.2.0      SparseM_1.77       quantreg_5.33     
    ## [28] pbkrtest_0.4-7     parallel_3.4.0     Rcpp_0.12.10      
    ## [31] scales_0.4.1       backports_1.0.5    jsonlite_1.4      
    ## [34] lme4_1.1-13        mnormt_1.5-5       hms_0.3           
    ## [37] digest_0.6.12      stringi_1.1.5      grid_3.4.0        
    ## [40] rprojroot_1.2      tools_3.4.0        magrittr_1.5      
    ## [43] lazyeval_0.2.0     car_2.1-4          MASS_7.3-47       
    ## [46] xml2_1.1.1         lubridate_1.6.0    minqa_1.2.4       
    ## [49] assertthat_0.2.0   rmarkdown_1.5      httr_1.2.1        
    ## [52] R6_2.2.0           nnet_7.3-12        nlme_3.1-131      
    ## [55] compiler_3.4.0

``` r
base_font_size <- 8
mm_to_pt <- 7 * 0.35 # so geom_text size is same as axis text
base_font_family <- "Helvetica"
base_line_size <- .25
theme_custom <- theme_cowplot(font_size = base_font_size, font_family = base_font_family, line_size = base_line_size)
theme_set(theme_custom)
```

Group analyses
==============

Load data
---------

All (meta)data are stored as `.csv` files in the `/data` folder.

### Eye tracking data

Here we load the `.csv` file with the processed eye tracking data, which was created in MATLAB. To recreate it from the raw data, run the `analysis.m` script. This scripts calls the `processEDF.m` function to parse the raw eye tracking data.

``` r
# Load eye tracking data into data frame
dataFile <- here("data", "sacc-tDCS_data.csv")
groupData <- read_csv(dataFile, col_names = TRUE, na = "NaN", progress = FALSE, col_types = cols(
  stimulation = col_factor(c("anodal","cathodal")),
  leg = col_factor(c("pre","tDCS","post")),
  type = col_factor(c("lateral","center")),
  direction = col_factor(c("left","right")) 
))
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

### Metadata

#### Session info

``` r
# Load eye tracking data into data frame
dataFile <- here("data", "session_info.csv")
sessionData <- read_csv2(dataFile, col_names = TRUE, progress = FALSE, col_types = cols(
  session = col_factor(c("first","second")),
  stimulation = col_factor(c("anodal","cathodal"))
))
```

``` r
kable(head(sessionData))
```

| subject | session | stimulation | date     | day       | time     |
|:--------|:--------|:------------|:---------|:----------|:---------|
| S01     | first   | cathodal    | 10/07/15 | Friday    | 11:30:00 |
| S01     | second  | anodal      | 16/07/15 | Thursday  | 09:00:00 |
| S02     | first   | anodal      | 10/07/15 | Friday    | 09:00:00 |
| S02     | second  | cathodal    | 13/07/15 | Monday    | 09:00:00 |
| S03     | first   | cathodal    | 13/07/15 | Monday    | 11:30:00 |
| S03     | second  | anodal      | 15/07/15 | Wednesday | 13:00:00 |

-   **subject**: subject ID
-   **session**: Whether data are from the `first`) or `second` session
-   **stimulation**: Whether data are from the `anodal` or `cathodal` session
-   **date**: YY/MM/DD date the session took place
-   **day**: day of the week the session took place
-   **time**: time of day the session took place

#### Subject info

``` r
# Load eye tracking data into data frame
dataFile <- here("data", "subject_info.csv")
subjectData <- read_csv2(dataFile, col_names = TRUE, progress = FALSE, col_types = cols(
  session.order = col_factor(c("first.anodal", "first.cathodal"))
))
```

``` r
kable(head(subjectData))
```

| subject | session.order  | gender |  age| dominant.eye |
|:--------|:---------------|:-------|----:|:-------------|
| S01     | first.cathodal | female |   23| right        |
| S02     | first.anodal   | male   |   27| left         |
| S03     | first.cathodal | male   |   34| right        |
| S04     | first.anodal   | male   |   NA| right        |
| S05     | first.cathodal | male   |   24| right        |
| S06     | first.cathodal | female |   29| left         |

-   **subject**: subject ID
-   **session.order**: Whether subject had anodal stimulation in the first session (`first.anodal`) or cathodal stimulation in the first session (`first.cathodal`)
-   **gender**
-   **age**: in years
-   **dominant.eye**: result of eye dominance test

Preprocess eye data
-------------------

*See the `inspection.nb.html` notebook for more explanation and exploration of preprocessing, outlier saccades and subject exclusion*

### Reject outlier saccades

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

``` r
# Mark outliers
groupData <- mutate(groupData, outlier = "non.outlier", # fill vector for all saccades
                    outlier = ifelse(latency < tooFast, "fast", outlier), # mark too fast saccades as "fast"
                    outlier = ifelse(latency > tooSlow, "slow", outlier), # mark too slow saccades as "slow"
                    outlier = ifelse(deviation.start > badFix, "fixation", outlier), # mark bad fixations as "fixation"
                    outlier = ifelse(deviation.end.x > badSacc, "saccade", outlier), # mark inaccurate saccades as "saccade"
                    outlier = ifelse(is.na(latency), "none", outlier) # mark absence of saccade as "none"
                    )
```

``` r
preproc <- filter(groupData, outlier == "non.outlier")
```

### Exclude participants

-   S21 and S25 were tested &lt; 48h apart
-   S16, S22 and S28 have fewer than 50 saccades per condition after trial rejection

``` r
subs2exclude <- c("S21","S25","S16","S22","S28")
```

``` r
preproc <- filter(preproc, !(subject %in% subs2exclude))
```

### Collapse data across 15 minute intervals

Cut the post-block into two so we have four 15-minute intervals: one before, one during, and two after stimulation.

``` r
preproc <- preproc %>%
  mutate(leg = as.character(leg), # cannot edit leg if it's still a factor
         leg = replace(leg, leg == "post" & block <= 3, "post.1"),
         leg = replace(leg, block > 3, "post.2"),
         leg = factor(leg, levels = c("pre", "tDCS", "post.1", "post.2")) # refactor and order levels
         )
```

Subject and session descriptives
--------------------------------

### Sample demographics

For S23, no eye data was collected because eye tracker could not be calibrated, so here we don't include their subject/session data.

**Gender of remaining participants:**

``` r
subjectData %>%
  filter(!(subject %in% c(subs2exclude, "S23"))) %>% # remove rows with these subject
  group_by(gender) %>% # for each gender
  summarise(count = n_distinct(subject)) %>%
  kable(.)# count number of subjects
```

| gender |  count|
|:-------|------:|
| female |     14|
| male   |     12|

**Age:**

``` r
subjectData %>%
  filter(!(subject %in% c(subs2exclude, "S23"))) %>%
  summarise_at(vars(age), funs(mean, min, max, sd), na.rm = TRUE) %>% # apply summary functions to age column
  kable(.)
```

|      mean|  min|  max|        sd|
|---------:|----:|----:|---------:|
|  25.90476|   21|   34|  3.419134|

### tDCS session order

``` r
sessionData %>%
  filter(!(subject %in% c(subs2exclude, "S23"))) %>%
  group_by(session, stimulation) %>%
  summarise(count = n_distinct(subject)) %>% 
  kable(.)
```

| session | stimulation |  count|
|:--------|:------------|------:|
| first   | anodal      |     14|
| first   | cathodal    |     12|
| second  | anodal      |     12|
| second  | cathodal    |     14|

### Saccade counts

Average number of saccades per type after rejection:

``` r
preproc %>%
  group_by(subject,stimulation,leg,direction,type) %>% # for each cell
  summarise(saccades = n()) %>% # count how many saccacdes there are left
  group_by(type) %>% 
  summarise_at(vars(saccades), funs(mean, min, max, sd)) %>% # compute summary statistics per type
  kable(.)
```

| type    |      mean|  min|  max|         sd|
|:--------|---------:|----:|----:|----------:|
| lateral |  175.2284|  142|  180|   6.244285|
| center  |  155.5529|   74|  180|  20.816516|

Percentage of outlier saccades per type out of all sacccades:

``` r
n_total <- nrow(filter(groupData, !(subject %in% subs2exclude)))  # total amount of saccades across all included sessions/subjects

groupData %>%
  filter(!(subject %in% subs2exclude), !(outlier %in% c("non.outlier", "none"))) %>%
  group_by(subject,stimulation,leg,direction,type,outlier) %>%
  summarize(outlier_count = n()) %>% # for each condition and subject, count how many (non)outliers there are
  group_by(outlier) %>% # for each outlier type
  summarise(percentage = sum(outlier_count) / n_total *100) %>% # calculate percentage
  kable(.)
```

| outlier  |  percentage|
|:---------|-----------:|
| fast     |   1.9958600|
| fixation |   2.5848024|
| saccade  |   0.1555823|
| slow     |   0.1161859|

Saccade latency
---------------

*See the `median_latency.nb.html` notebook for more explanation and exploration of the median saccade latency data*

Compute the median for each 3 blocks

``` r
latencyMedianLeg <- preproc %>%
  group_by(subject,stimulation,direction,type) %>% 
  summarise(baseline = median(latency[leg == "pre"]), # take average of 3 blocks, make new column
            tDCS = median(latency[leg == "tDCS"]),
            post.1 = median(latency[leg == "post.1"]),
            post.2 = median(latency[leg == "post.2"])) %>%
  gather(leg,latency,baseline,tDCS,post.1,post.2) %>% # gather new columns to use as factor
  mutate(leg = factor(leg, levels = c("baseline", "tDCS", "post.1", "post.2"))) # reorder factor levels
```

Make a separate data frame where the median baseline latency is subtracted from subsequent blocks

``` r
latencyMedianBaseline <- latencyMedianLeg %>%
  group_by(subject,stimulation,direction,type) %>% # for each condition, subtract baseline scores and make new columns
  summarise(baseline = latency[leg == "baseline"] - latency[leg == "baseline"],
           tDCS = latency[leg == "tDCS"] - latency[leg == "baseline"],
           post.1 = latency[leg == "post.1"] - latency[leg == "baseline"],
           post.2 = latency[leg == "post.2"] - latency[leg == "baseline"]) %>%
  gather(leg, latency.baseline, baseline, tDCS, post.1, post.2)  %>% # gather new columns to use as factor 
  mutate(leg = factor(leg, levels = c("baseline", "tDCS", "post.1", "post.2"))) # reorder factor levels
```

### Figure 4

``` r
baselineLabelLatency <- latencyMedianLeg %>%
  group_by(stimulation,direction,type,leg) %>%
  filter(leg == "baseline") %>%
  summarise(latency.baseline = round(mean(latency)))
```

Plot median *lateral saccade* latency per participant

``` r
plot_sub_latency_lateral <- ggplot(filter(latencyMedianLeg, type == "lateral"), aes(leg, latency)) +
  facet_grid(stimulation ~ direction) +
  geom_line(aes(group = subject, color = subject), alpha = 0.5) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), linetype = direction), geom = "line") +
  stat_summary(fun.y = mean, geom = "point", aes(shape = stimulation)) +
  guides(colour = FALSE) +
  scale_fill_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous("median latency (ms)", breaks = seq(100,200,20), limits = c(90,220)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank(), legend.position = "none") +
  ggtitle("Lateral saccades")
```

Plot median *center saccade* latency per participant

``` r
plot_sub_latency_center <- ggplot(filter(latencyMedianLeg, type == "center"), aes(leg, latency)) +
  facet_grid(stimulation ~ direction) +
  geom_line(aes(group = subject, color = subject), alpha = 0.5) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), linetype = direction), geom = "line") +
  stat_summary(fun.y = mean, geom = "point", aes(shape = stimulation)) +
  scale_fill_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous("median latency (ms)", breaks = seq(100,200,20), limits = c(90,220)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank(), legend.position = "none") +
  ggtitle("Center saccades")
```

Plot group mean *lateral saccade* latency as change from baseline

``` r
plot_group_latency_lateral <- ggplot(filter(latencyMedianBaseline, type == "lateral"), aes(leg, latency.baseline)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = base_line_size) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "line", position = position_dodge(width = 0.5), size = .75) + 
  stat_summary(fun.data = mean_cl_normal, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "linerange", position = position_dodge(width = 0.5), show.legend = FALSE) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, shape = stimulation), geom = "point", position = position_dodge(width = 0.5), size = 2) +
  geom_text(data = subset(baselineLabelLatency, type=="lateral"), aes(y = c(2,1,4,3), label=latency.baseline, color = stimulation), position = position_dodge(width = 0.5), size = mm_to_pt) +
  scale_colour_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous("latency change from baseline (ms)", breaks = seq(-8,8,2)) +
  coord_cartesian(ylim = c(-8,8)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank())
```

Plot group mean *center saccade* latency as change from baseline

``` r
plot_group_latency_center <- ggplot(filter(latencyMedianBaseline, type == "center"), aes(leg, latency.baseline)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = base_line_size) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "line", position = position_dodge(width = 0.5), size = .75) + 
  stat_summary(fun.data = mean_cl_normal, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "linerange", position = position_dodge(width = 0.5), show.legend = FALSE) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, shape = stimulation), geom = "point", position = position_dodge(width = 0.5), size = 2) +
  geom_text(data = subset(baselineLabelLatency, type=="center"), aes(y = c(2,1,3,4), label=latency.baseline, color = stimulation), position = position_dodge(width = 0.5), size = mm_to_pt) +
  scale_colour_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_discrete(labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous("latency change from baseline (ms)", breaks = seq(-8,8,2)) +
  coord_cartesian(ylim = c(-8,8)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank())
```

Combine all the plots:

``` r
legend_fig4 <- get_legend(plot_group_latency_lateral)
figure_4 <- plot_grid(
  plot_sub_latency_lateral, plot_sub_latency_center,
  plot_group_latency_lateral + theme(legend.position = "none"),
  plot_group_latency_center + theme(legend.position = "none"), rel_heights = c(1,.75))
figure_4 <- plot_grid(figure_4,legend_fig4,rel_widths = c(1,1/10))
figure_4
```

<img src="sacc-tDCS_outline_files/figure-markdown_github/Figure 4-1.png" width="75%" style="display: block; margin: auto;" />

Save the plot:

``` r
ggsave("fig/figure_4.pdf", plot = figure_4, width = 180, height = 112.5, units = "mm")
ggsave("fig/figure_4.png", plot = figure_4, width = 180, height = 112.5, units = "mm")
```

### Statistics

#### Contralateral saccades in the anodal session

Means and SDs for effect of *anodal* tDCS on the latency changes in *contralateral* (here: left) *lateral* saccades, as this is what [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045) specifically found:

``` r
latencyMedianBaseline %>%
  filter(type == "lateral", stimulation == "anodal", direction == "left", leg != "baseline") %>%
  group_by(leg) %>%
  summarise_at(vars(latency.baseline), funs(mean, sd)) %>%
  kable(.)
```

| leg    |        mean|        sd|
|:-------|-----------:|---------:|
| tDCS   |  -0.1730769|  5.342176|
| post.1 |  -0.6153846|  7.134855|
| post.2 |   0.9615385|  9.651863|

#### Baseline tests

``` r
latencyMedianLeg %>%
  filter(leg == "baseline") %>%
  group_by(direction,type) %>% # for each of these 4 pairs
  nest() %>% # create a separate data frame of all remaining columns. These data frames are stored in a list-column named "data"
  mutate(stats = map(data, ~t.test(formula = latency~stimulation, paired = TRUE, data =.))) %>% # run t-test on the data frames
  mutate(tidy_model = map(stats, tidy)) %>% # force the four test outputs to also be data frames in a list column ("tidy_model")
  unnest(tidy_model, .drop = TRUE) %>% # unpack the list-column with results, drop the list-column we used to organize the data
  kable(.)
```

| direction | type    |  estimate|  statistic|    p.value|  parameter|    conf.low|  conf.high| method        | alternative |
|:----------|:--------|---------:|----------:|----------:|----------:|-----------:|----------:|:--------------|:------------|
| left      | lateral |  2.403846|   1.024429|  0.3154437|         25|  -2.4289077|   7.236600| Paired t-test | two.sided   |
| left      | center  |  5.557692|   1.910225|  0.0676402|         25|  -0.4344207|  11.549805| Paired t-test | two.sided   |
| right     | lateral |  2.769231|   1.060619|  0.2990028|         25|  -2.6081375|   8.146599| Paired t-test | two.sided   |
| right     | center  |  4.211538|   1.657484|  0.1099194|         25|  -1.0215893|   9.444666| Paired t-test | two.sided   |

#### Lateral saccades

##### Classical ANOVA

Repeated measures ANOVA matching [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045).

**Data**: Baseline subtracted, lateral saccades

**Dependent measure**: saccade latency

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

Prepare data:

``` r
medianLatencyLateralStats <- latencyMedianBaseline %>%
  ungroup() %>%
  filter(leg != "baseline", type == "lateral") %>%
  select(-type) %>%
  mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) %>%
  mutate(subject = factor(subject))
```

``` r
aovMedianLatencyLateral <- ezANOVA(data = data.frame(medianLatencyLateralStats), dv = latency.baseline, wid = subject, within = .(stimulation,leg,direction), type = 3)

kable(aovMedianLatencyLateral$ANOVA)
```

|     | Effect                    |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  0.8044339|  0.3783271|          |  0.0089853|
| 3   | leg                       |    2|   50|  3.4589912|  0.0391753| \*       |  0.0164093|
| 4   | direction                 |    1|   25|  0.0864110|  0.7712172|          |  0.0003663|
| 5   | stimulation:leg           |    2|   50|  0.4014143|  0.6715104|          |  0.0007505|
| 6   | stimulation:direction     |    1|   25|  0.5221926|  0.4766141|          |  0.0012979|
| 7   | leg:direction             |    2|   50|  3.6626747|  0.0327780| \*       |  0.0023536|
| 8   | stimulation:leg:direction |    2|   50|  2.5875456|  0.0852459|          |  0.0029246|

``` r
kable(aovMedianLatencyLateral$`Mauchly's Test for Sphericity`)
```

|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.5213492|  0.0004032| \*       |
| 5   | stimulation:leg           |  0.6360361|  0.0043831| \*       |
| 7   | leg:direction             |  0.7573179|  0.0355909| \*       |
| 8   | stimulation:leg:direction |  0.9123864|  0.3327710|          |

``` r
kable(aovMedianLatencyLateral$`Sphericity Corrections`)
```

|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.6762922|  0.0597258|                |  0.7012858|  0.0578252|                |
| 5   | stimulation:leg           |  0.7331572|  0.6087290|                |  0.7674993|  0.6179798|                |
| 7   | leg:direction             |  0.8047111|  0.0434681| \*             |  0.8517313|  0.0406111| \*             |
| 8   | stimulation:leg:direction |  0.9194442|  0.0903598|                |  0.9889658|  0.0859306|                |

##### Bayesian ANOVA

Same design as the classical ANOVA

Bayes factors agains the null model:

``` r
bfMedianLatencyLateral = anovaBF(latency.baseline~stimulation*leg*direction+subject, data = data.frame(medianLatencyLateralStats), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) 
bfMedianLatencyLateral = sort(bfMedianLatencyLateral, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfMedianLatencyLateral), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |         bf|
|-------------------------------------------------------------------------------------------------------------------------------|----------:|
| leg + subject                                                                                                                 |  0.8748482|
| stimulation + subject                                                                                                         |  0.7443607|
| stimulation + leg + subject                                                                                                   |  0.6528811|
| direction + subject                                                                                                           |  0.1330391|
| direction + leg + subject                                                                                                     |  0.1145413|
| stimulation + direction + subject                                                                                             |  0.1013966|
| stimulation + direction + leg + subject                                                                                       |  0.0863803|
| stimulation + leg + stimulation:leg + subject                                                                                 |  0.0471885|
| stimulation + direction + stimulation:direction + subject                                                                     |  0.0221681|
| stimulation + direction + stimulation:direction + leg + subject                                                               |  0.0192023|
| direction + leg + direction:leg + subject                                                                                     |  0.0105721|
| stimulation + direction + leg + direction:leg + subject                                                                       |  0.0088664|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |  0.0062142|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |  0.0019030|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |  0.0013428|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |  0.0006393|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |  0.0001362|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |  0.0000247|

Inclusion Bayes factor across matched models:

``` r
kable(inclusionBF(bfMedianLatencyLateral, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |     0.7472502|
| direction                 |     0.1330323|
| stimulation:direction     |     0.2199172|
| leg                       |     0.8735054|
| stimulation:leg           |     0.0721770|
| direction:leg             |     0.0971402|
| stimulation:direction:leg |     0.1813377|

#### Center saccades

##### Classical ANOVA

Repeated measures ANOVA matching [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045).

**Data**: Baseline subtracted, center saccades

**Dependent measure**: saccade latency

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

Prepare data:

``` r
medianLatencyCenterStats <- latencyMedianBaseline %>%
  ungroup() %>%
  filter(leg != "baseline", type == "center") %>%
  select(-type) %>%
  mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) %>%
  mutate(subject = factor(subject))
```

``` r
aovMedianLatencyCenter <- ezANOVA(data = data.frame(medianLatencyCenterStats), dv = latency.baseline, wid = subject, within = .(stimulation,leg,direction), type = 3)

kable(aovMedianLatencyCenter$ANOVA)
```

|     | Effect                    |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  3.0937191|  0.0908311|          |  0.0230246|
| 3   | leg                       |    2|   50|  2.3789819|  0.1030547|          |  0.0062189|
| 4   | direction                 |    1|   25|  4.7775590|  0.0384115| \*       |  0.0158632|
| 5   | stimulation:leg           |    2|   50|  0.1055936|  0.8999904|          |  0.0001735|
| 6   | stimulation:direction     |    1|   25|  1.8802459|  0.1824903|          |  0.0055307|
| 7   | leg:direction             |    2|   50|  3.1018543|  0.0537187|          |  0.0016853|
| 8   | stimulation:leg:direction |    2|   50|  1.9636026|  0.1510262|          |  0.0012036|

``` r
kable(aovMedianLatencyCenter$`Mauchly's Test for Sphericity`)
```

|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.5173954|  0.0003680| \*       |
| 5   | stimulation:leg           |  0.8690648|  0.1856204|          |
| 7   | leg:direction             |  0.9230394|  0.3825099|          |
| 8   | stimulation:leg:direction |  0.7823535|  0.0525819|          |

``` r
kable(aovMedianLatencyCenter$`Sphericity Corrections`)
```

|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.6744887|  0.1241893|                |  0.6991962|  0.1225297|                |
| 5   | stimulation:leg           |  0.8842239|  0.8774896|                |  0.9465498|  0.8902341|                |
| 7   | leg:direction             |  0.9285390|  0.0579604|                |  0.9999608|  0.0537210|                |
| 8   | stimulation:leg:direction |  0.8212564|  0.1600417|                |  0.8713552|  0.1575428|                |

##### Bayesian ANOVA

Same design as the classical ANOVA

Bayes factors agains the null model:

``` r
bfMedianLatencyCenter = anovaBF(latency.baseline~stimulation*leg*direction+subject, data = data.frame(medianLatencyCenterStats), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) 
bfMedianLatencyCenter = sort(bfMedianLatencyCenter, decreasing = TRUE) # sort such that winning model is at the top

kable(select(extractBF(bfMedianLatencyCenter), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |           bf|
|-------------------------------------------------------------------------------------------------------------------------------|------------:|
| stimulation + direction + subject                                                                                             |  545.2420547|
| stimulation + direction + stimulation:direction + subject                                                                     |  406.7143803|
| stimulation + direction + leg + subject                                                                                       |   97.7977620|
| stimulation + direction + stimulation:direction + leg + subject                                                               |   75.6169133|
| stimulation + subject                                                                                                         |   56.6669318|
| stimulation + leg + subject                                                                                                   |    9.6277038|
| stimulation + direction + leg + direction:leg + subject                                                                       |    9.2344705|
| direction + subject                                                                                                           |    7.8800046|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |    7.2861399|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |    6.2990737|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |    4.7452145|
| direction + leg + subject                                                                                                     |    1.3079424|
| stimulation + leg + stimulation:leg + subject                                                                                 |    0.6308269|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |    0.6002592|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |    0.4743342|
| leg + subject                                                                                                                 |    0.1579005|
| direction + leg + direction:leg + subject                                                                                     |    0.1251356|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |    0.0714910|

``` r
kable(select(extractBF(bfMedianLatencyCenter), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |           bf|
|-------------------------------------------------------------------------------------------------------------------------------|------------:|
| stimulation + direction + subject                                                                                             |  545.2420547|
| stimulation + direction + stimulation:direction + subject                                                                     |  406.7143803|
| stimulation + direction + leg + subject                                                                                       |   97.7977620|
| stimulation + direction + stimulation:direction + leg + subject                                                               |   75.6169133|
| stimulation + subject                                                                                                         |   56.6669318|
| stimulation + leg + subject                                                                                                   |    9.6277038|
| stimulation + direction + leg + direction:leg + subject                                                                       |    9.2344705|
| direction + subject                                                                                                           |    7.8800046|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |    7.2861399|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |    6.2990737|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |    4.7452145|
| direction + leg + subject                                                                                                     |    1.3079424|
| stimulation + leg + stimulation:leg + subject                                                                                 |    0.6308269|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |    0.6002592|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |    0.4743342|
| leg + subject                                                                                                                 |    0.1579005|
| direction + leg + direction:leg + subject                                                                                     |    0.1251356|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |    0.0714910|

Inclusion Bayes factor across matched models:

``` r
kable(inclusionBF(bfMedianLatencyCenter, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |    68.6247813|
| direction                 |     9.6723606|
| stimulation:direction     |     0.7506929|
| leg                       |     0.1813343|
| stimulation:leg           |     0.0638881|
| direction:leg             |     0.0953902|
| stimulation:direction:leg |     0.1507187|

###### Without S01

``` r
medianLatencyCenterNoS01 <- medianLatencyCenterStats %>%
  filter(subject != "S01") %>%
  mutate(subject = factor(subject))

bfMedianLatencyCenterNoS01 = anovaBF(latency.baseline~stimulation*leg*direction+subject, data = data.frame(medianLatencyCenterNoS01), whichModels="withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiCenterNoS01 <- sort(bfMedianLatencyCenterNoS01, decreasing = TRUE) # sort such that winning model is at the top
kable(select(extractBF(bfMedianLatencyCenterNoS01), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |            bf|
|-------------------------------------------------------------------------------------------------------------------------------|-------------:|
| stimulation + subject                                                                                                         |     1.9799111|
| direction + subject                                                                                                           |   599.1336816|
| stimulation + direction + subject                                                                                             |  1428.4811363|
| stimulation + direction + stimulation:direction + subject                                                                     |   347.6358459|
| leg + subject                                                                                                                 |     0.3008156|
| stimulation + leg + subject                                                                                                   |     0.6178388|
| direction + leg + subject                                                                                                     |   207.6795625|
| stimulation + direction + leg + subject                                                                                       |   532.0914361|
| stimulation + direction + stimulation:direction + leg + subject                                                               |   129.8412058|
| stimulation + leg + stimulation:leg + subject                                                                                 |     0.0411408|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |    33.4057211|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |     8.4260419|
| direction + leg + direction:leg + subject                                                                                     |    24.1946777|
| stimulation + direction + leg + direction:leg + subject                                                                       |    62.6160130|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |    15.1277324|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |     4.2072188|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |     1.2480044|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |     0.1690676|

``` r
kable(inclusionBF(bfMedianLatencyCenterNoS01, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |     2.4339362|
| direction                 |   710.9138133|
| stimulation:direction     |     0.2437298|
| leg                       |     0.3660414|
| stimulation:leg           |     0.0639315|
| direction:leg             |     0.1178280|
| stimulation:direction:leg |     0.1354704|

##### Follow-up tests

###### Main effect of stimulation

Follow-up One-sample t-test:

``` r
medianLatencyCenterStats %>%
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(latency = mean(latency.baseline)) %>% # average over all other variables (df is now still grouped per stimulation)
  summarise_if(is.numeric, funs(list(tidy(t.test(.))))) %>%  # run one-sample t-test for each stimulation condition, return tidy data frames
  unnest() %>% # unpack the list-column with data frame for each test
  kable(.)
```

| stimulation |   estimate|   statistic|    p.value|  parameter|   conf.low|  conf.high| method            | alternative |
|:------------|----------:|-----------:|----------:|----------:|----------:|----------:|:------------------|:------------|
| anodal      |  -1.608974|  -0.9875463|  0.3328370|         25|  -4.964508|   1.746559| One Sample t-test | two.sided   |
| cathodal    |   1.096154|   0.8335773|  0.4124131|         25|  -1.612139|   3.804446| One Sample t-test | two.sided   |

Follow-up Bayesian one-sample t-test:

``` r
medianLatencyCenterStats %>%
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(latency = mean(latency.baseline)) %>% # average over all other variables
  spread(stimulation, latency) %>% # make separate columns with test data
  summarise_if(is.numeric, funs(extractBF(ttestBF(.), onlybf = TRUE))) %>% # run Bayesian t-test on each column, keeping only the BF
  gather(stimulation,BF,anodal,cathodal) %>% # make row for each stimulation condition
  kable(.)
```

| stimulation |         BF|
|:------------|----------:|
| anodal      |  0.3218926|
| cathodal    |  0.2840670|

Saccade endpoint deviation
--------------------------

*See the `accuracy.nb.html` notebook for more explanation and exploration of the saccade endpoint deviation data*

Calculate endpoint deviation for each trial (Euclidian distance from target, defined by x- and y-coordinate deviations)

``` r
# Calculate end point deviation
accData <- preproc %>%
  group_by(subject,stimulation,leg,block,trial,direction,type) %>%
  summarise(deviation.end = sqrt(deviation.end.x^2 + deviation.end.y^2))
```

Compute the mean for each 3 blocks

``` r
accMeanLeg <- accData %>%
  group_by(subject,stimulation,direction,type) %>% 
  summarise(baseline = mean(deviation.end[leg == "pre"]), # take average of 3 blocks, make new column
            tDCS = mean(deviation.end[leg == "tDCS"]),
            post.1 = mean(deviation.end[leg == "post.1"]),
            post.2 = mean(deviation.end[leg == "post.2"])) %>%
  gather(leg,deviation.end,baseline,tDCS,post.1,post.2) %>% # gather new columns to use as factor
  mutate(leg = factor(leg, levels = c("baseline", "tDCS", "post.1", "post.2"))) # reorder factor levels
```

Make a separate data frame where the mean deviation is subtracted from subsequent blocks

``` r
accMeanBaseline <- accMeanLeg %>%
  group_by(subject,stimulation,direction,type) %>% # for each condition, subtract baseline scores and make new columns
  summarise(baseline = deviation.end[leg == "baseline"] - deviation.end[leg == "baseline"],
           tDCS = deviation.end[leg == "tDCS"] - deviation.end[leg == "baseline"],
           post.1 = deviation.end[leg == "post.1"] - deviation.end[leg == "baseline"],
           post.2 = deviation.end[leg == "post.2"] - deviation.end[leg == "baseline"]) %>%
  gather(leg, deviation.end.baseline, baseline, tDCS, post.1, post.2)  %>% # gather new columns to use as factor 
  mutate(leg = factor(leg, levels = c("baseline", "tDCS", "post.1", "post.2"))) # reorder factor levels
```

### Figure 6

``` r
baselineLabelAcc <- accMeanLeg %>%
  group_by(stimulation,direction,type,leg) %>%
  filter(leg == "baseline") %>%
  summarise(deviation.end.baseline = round(mean(deviation.end), digits = 2))
```

Plot mean *lateral saccade* deviation per participant

``` r
plot_sub_acc_lateral <- ggplot(filter(accMeanLeg, type == "lateral"), aes(leg, deviation.end)) +
  facet_grid(stimulation ~ direction) +
  geom_line(aes(group = subject, color = subject), alpha = 0.5) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), linetype = direction), geom = "line") +
  stat_summary(fun.y = mean, geom = "point", aes(shape = stimulation)) +
  guides(colour = FALSE) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous(expression("mean deviation " ( degree)), breaks = seq(0.5,2,0.5), limits = c(0.3,2.2)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank(), legend.position = "none") +
  ggtitle("Lateral saccades")
```

Plot mean *center saccade* deviation per participant

``` r
plot_sub_acc_center <- ggplot(filter(accMeanLeg, type == "center"), aes(leg, deviation.end)) +
  facet_grid(stimulation ~ direction) +
  geom_line(aes(group = subject, color = subject), alpha = 0.5) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), linetype = direction), geom = "line") +
  stat_summary(fun.y = mean, geom = "point", aes(shape = stimulation)) +
  guides(colour = FALSE) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous(expression("mean deviation " ( degree)), breaks = seq(0.5,2,0.5), limits = c(0.3,2.2)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank(), legend.position = "none") +
  ggtitle("Center saccades")
```

Plot group mean *lateral saccade* deviation as change from baseline

``` r
plot_group_acc_lateral <- ggplot(filter(accMeanBaseline, type == "lateral"), aes(leg, deviation.end.baseline)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = base_line_size) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "line", position = position_dodge(width = 0.5), size = .75) + 
  stat_summary(fun.data = mean_cl_normal, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "linerange", position = position_dodge(width = 0.5), show.legend = FALSE) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, shape = stimulation), geom = "point", position = position_dodge(width = 0.5), size = 2) +
  geom_text(data = subset(baselineLabelAcc, type=="lateral"), aes(y = c(.05,.025,.1,.075), label=deviation.end.baseline, color = stimulation), position = position_dodge(width = 0.5), size = mm_to_pt) +
  scale_colour_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous(expression("deviation change from baseline " ( degree)), breaks = seq(-.2,.2,.1)) +
  coord_cartesian(ylim = c(-.22,.22)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank())
```

Plot group mean *center saccade* deviation as change from baseline

``` r
plot_group_acc_center <- ggplot(filter(accMeanBaseline, type == "center"), aes(leg, deviation.end.baseline)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = base_line_size) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "line", position = position_dodge(width = 0.5), size = .75) + 
  stat_summary(fun.data = mean_cl_normal, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "linerange", position = position_dodge(width = 0.5), show.legend = FALSE) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, shape = stimulation), geom = "point", position = position_dodge(width = 0.5), size = 2) +
  geom_text(data = subset(baselineLabelAcc, type=="center"), aes(y = c(.05,.025,.1,.075), label=deviation.end.baseline, color = stimulation), position = position_dodge(width = 0.5), size = mm_to_pt) +
  scale_colour_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous(expression("deviation change from baseline " ( degree)), breaks = seq(-.2,.2,.1)) +
  coord_cartesian(ylim = c(-.22,.22)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank())
```

Combine all the plots:

``` r
legend_fig6 <- get_legend(plot_group_acc_lateral)
figure_6 <- plot_grid(
  plot_sub_acc_lateral, plot_sub_acc_center,
  plot_group_acc_lateral + theme(legend.position = "none"),
  plot_group_acc_center + theme(legend.position = "none"), rel_heights = c(1,.75))
figure_6 <- plot_grid(figure_6,legend_fig6,rel_widths = c(1,1/10))
figure_6
```

<img src="sacc-tDCS_outline_files/figure-markdown_github/Figure 6-1.png" width="75%" style="display: block; margin: auto;" />

Save the plot:

``` r
ggsave("fig/figure_6.pdf", plot = figure_6, width = 180, height = 112.5, units = "mm")
ggsave("fig/figure_6.png", plot = figure_6, width = 180, height = 112.5, units = "mm")
```

### Statistics

#### Baseline tests

``` r
accMeanLeg %>%
  filter(leg == "baseline") %>%
  group_by(direction,type) %>% 
  nest() %>% 
  mutate(stats = map(data, ~t.test(formula = deviation.end~stimulation, paired = TRUE, data =.))) %>% # run t-test on the data frames
  mutate(tidy_model = map(stats, tidy)) %>%
  unnest(tidy_model, .drop = TRUE) %>% 
  kable(.)
```

| direction | type    |    estimate|   statistic|    p.value|  parameter|    conf.low|   conf.high| method        | alternative |
|:----------|:--------|-----------:|-----------:|----------:|----------:|-----------:|-----------:|:--------------|:------------|
| left      | lateral |  -0.0434955|  -0.9392156|  0.3566066|         25|  -0.1388738|   0.0518827| Paired t-test | two.sided   |
| left      | center  |  -0.1052918|  -2.3796609|  0.0252737|         25|  -0.1964194|  -0.0141643| Paired t-test | two.sided   |
| right     | lateral |  -0.0141028|  -0.4198800|  0.6781598|         25|  -0.0832780|   0.0550724| Paired t-test | two.sided   |
| right     | center  |  -0.0576332|  -1.9185260|  0.0665325|         25|  -0.1195025|   0.0042361| Paired t-test | two.sided   |

#### Differences after baseline

Because there might be a significant baseline difference for center saccades, also print the numerical values for the blocks after, to see whether they continue to differ

``` r
accMeanLeg %>%
  filter(leg != "baseline", type == "center") %>%
  group_by(stimulation,direction,leg) %>% 
  summarise(mean = mean(deviation.end)) %>%
  kable(.)            
```

| stimulation | direction | leg    |       mean|
|:------------|:----------|:-------|----------:|
| anodal      | left      | tDCS   |  0.9149277|
| anodal      | left      | post.1 |  0.8752206|
| anodal      | left      | post.2 |  0.8978057|
| anodal      | right     | tDCS   |  0.7192357|
| anodal      | right     | post.1 |  0.7110237|
| anodal      | right     | post.2 |  0.7465372|
| cathodal    | left      | tDCS   |  0.9129682|
| cathodal    | left      | post.1 |  0.8395323|
| cathodal    | left      | post.2 |  0.8357657|
| cathodal    | right     | tDCS   |  0.7227287|
| cathodal    | right     | post.1 |  0.7062460|
| cathodal    | right     | post.2 |  0.7483820|

#### Lateral saccades

##### Classical ANOVAs

Repeated measures ANOVA matching [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045).

**Data**: Baseline subtracted, lateral saccades

**Dependent measure**: saccade endpoint deviation

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

Prepare data:

``` r
meanAccLateralStats <- accMeanBaseline %>%
  ungroup() %>%
  filter(leg != "baseline", type == "lateral") %>%
  select(-type) %>%
  mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) %>%
  mutate(subject = factor(subject))
```

``` r
aovMeanAccLateral <- ezANOVA(data = data.frame(meanAccLateralStats), dv = deviation.end.baseline, wid = subject, within = .(stimulation,leg,direction), type = 3)

kable(aovMeanAccLateral$ANOVA)
```

|     | Effect                    |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  2.0266072|  0.1669269|          |  0.0178886|
| 3   | leg                       |    2|   50|  0.8985792|  0.4136202|          |  0.0021712|
| 4   | direction                 |    1|   25|  1.5739544|  0.2212365|          |  0.0059023|
| 5   | stimulation:leg           |    2|   50|  0.5863132|  0.5601540|          |  0.0015026|
| 6   | stimulation:direction     |    1|   25|  0.1283222|  0.7231851|          |  0.0005433|
| 7   | leg:direction             |    2|   50|  0.1400779|  0.8696305|          |  0.0002306|
| 8   | stimulation:leg:direction |    2|   50|  0.2773853|  0.7589209|          |  0.0002720|

``` r
kable(aovMeanAccLateral$`Mauchly's Test for Sphericity`)
```

|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.6913883|  0.0119307| \*       |
| 5   | stimulation:leg           |  0.9959286|  0.9522226|          |
| 7   | leg:direction             |  0.9738561|  0.7276748|          |
| 8   | stimulation:leg:direction |  0.7556773|  0.0346766| \*       |

``` r
kable(aovMeanAccLateral$`Sphericity Corrections`)
```

|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.7641686|  0.3907386|                |  0.8038793|  0.3951392|                |
| 5   | stimulation:leg           |  0.9959451|  0.5595027|                |  1.0819912|  0.5601540|                |
| 7   | leg:direction             |  0.9745222|  0.8646202|                |  1.0558164|  0.8696305|                |
| 8   | stimulation:leg:direction |  0.8036501|  0.7108179|                |  0.8504748|  0.7234385|                |

##### Bayesian ANOVA

Same design as the classical ANOVA (without order effect)

Bayes factors agains the null model:

``` r
bfMeanAccLateral = anovaBF(deviation.end.baseline~stimulation*leg*direction+subject, data = data.frame(meanAccLateralStats), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) 
bfMeanAccLateral = sort(bfMeanAccLateral, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfMeanAccLateral), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |         bf|
|-------------------------------------------------------------------------------------------------------------------------------|----------:|
| stimulation + subject                                                                                                         |  6.4592216|
| stimulation + direction + subject                                                                                             |  2.9022767|
| stimulation + direction + stimulation:direction + subject                                                                     |  0.5603978|
| direction + subject                                                                                                           |  0.4464701|
| stimulation + leg + subject                                                                                                   |  0.3425893|
| stimulation + direction + leg + subject                                                                                       |  0.1551713|
| leg + subject                                                                                                                 |  0.0535084|
| stimulation + direction + stimulation:direction + leg + subject                                                               |  0.0303224|
| stimulation + leg + stimulation:leg + subject                                                                                 |  0.0294101|
| direction + leg + subject                                                                                                     |  0.0237864|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |  0.0131999|
| stimulation + direction + leg + direction:leg + subject                                                                       |  0.0107798|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |  0.0025799|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |  0.0021183|
| direction + leg + direction:leg + subject                                                                                     |  0.0015757|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |  0.0008764|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |  0.0001743|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |  0.0000190|

Inclusion Bayes factor across matched models:

``` r
kable(inclusionBF(bfMeanAccLateral, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |     6.4707119|
| direction                 |     0.4490838|
| stimulation:direction     |     0.1932297|
| leg                       |     0.0532511|
| stimulation:leg           |     0.0854756|
| direction:leg             |     0.0689793|
| stimulation:direction:leg |     0.1089252|

#### Center saccades

##### Classical ANOVA

Repeated measures ANOVA matching [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045).

**Data**: Baseline subtracted, center saccades

**Dependent measure**: saccade endpoint deviation

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

Prepare data:

``` r
meanAccCenterStats <- accMeanBaseline %>%
  ungroup() %>%
  filter(leg != "baseline", type == "center") %>%
  select(-type) %>%
  mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) %>%
  mutate(subject = factor(subject))
```

``` r
aovMeanAccCenter <- ezANOVA(data = data.frame(meanAccCenterStats), dv = deviation.end.baseline, wid = subject, within = .(stimulation,leg,direction), type = 3)

kable(aovMeanAccCenter$ANOVA)
```

|     | Effect                    |  DFn|  DFd|           F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|-----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  10.3422560|  0.0035735| \*       |  0.0698321|
| 3   | leg                       |    2|   50|   2.4132849|  0.0998788|          |  0.0064777|
| 4   | direction                 |    1|   25|   0.4250384|  0.5203830|          |  0.0041044|
| 5   | stimulation:leg           |    2|   50|   0.6101220|  0.5472793|          |  0.0012937|
| 6   | stimulation:direction     |    1|   25|   2.7995910|  0.1067576|          |  0.0126867|
| 7   | leg:direction             |    2|   50|   3.8854661|  0.0270094| \*       |  0.0071117|
| 8   | stimulation:leg:direction |    2|   50|   0.5883554|  0.5590374|          |  0.0011173|

``` r
kable(aovMeanAccCenter$`Mauchly's Test for Sphericity`)
```

|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.8436804|  0.1300575|          |
| 5   | stimulation:leg           |  0.9603861|  0.6156733|          |
| 7   | leg:direction             |  0.6656519|  0.0075677| \*       |
| 8   | stimulation:leg:direction |  0.8210183|  0.0938067|          |

``` r
kable(aovMeanAccCenter$`Sphericity Corrections`)
```

|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.8648128|  0.1084443|                |  0.9232826|  0.1046730|                |
| 5   | stimulation:leg           |  0.9618956|  0.5412918|                |  1.0404345|  0.5472793|                |
| 7   | leg:direction             |  0.7494296|  0.0403283| \*             |  0.7865648|  0.0380011| \*             |
| 8   | stimulation:leg:direction |  0.8481896|  0.5328387|                |  0.9034188|  0.5428502|                |

##### Bayesian ANOVA

Same design as the classical ANOVA

Bayes factors agains the null model:

``` r
bfMeanAccCenter = anovaBF(deviation.end.baseline~stimulation*leg*direction+subject, data = data.frame(meanAccCenterStats), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) 
bfMeanAccCenter = sort(bfMeanAccCenter, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfMeanAccCenter), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |            bf|
|-------------------------------------------------------------------------------------------------------------------------------|-------------:|
| stimulation + subject                                                                                                         |  4.533765e+04|
| stimulation + direction + stimulation:direction + subject                                                                     |  1.604652e+04|
| stimulation + direction + subject                                                                                             |  1.059494e+04|
| stimulation + leg + subject                                                                                                   |  4.256991e+03|
| stimulation + direction + stimulation:direction + leg + subject                                                               |  1.679144e+03|
| stimulation + direction + leg + subject                                                                                       |  1.079278e+03|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |  3.443323e+02|
| stimulation + leg + stimulation:leg + subject                                                                                 |  3.366993e+02|
| stimulation + direction + leg + direction:leg + subject                                                                       |  2.149325e+02|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |  1.326190e+02|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |  8.356877e+01|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |  2.639999e+01|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |  1.668678e+01|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |  3.503442e+00|
| direction + subject                                                                                                           |  2.378285e-01|
| leg + subject                                                                                                                 |  9.354690e-02|
| direction + leg + subject                                                                                                     |  2.346230e-02|
| direction + leg + direction:leg + subject                                                                                     |  4.146300e-03|

Inclusion Bayes factor across matched models:

``` r
kable(inclusionBF(bfMeanAccCenter, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |  4.524247e+04|
| direction                 |  2.354792e-01|
| stimulation:direction     |  1.520427e+00|
| leg                       |  9.746450e-02|
| stimulation:leg           |  7.867980e-02|
| direction:leg             |  2.024975e-01|
| stimulation:direction:leg |  1.327062e-01|

##### Follow-up tests

###### Main effect of stimulation

Follow-up One-sample t-test:

``` r
meanAccCenterStats %>%
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(deviation = mean(deviation.end.baseline)) %>% # average over all other variables (df is now still grouped per stimulation)
  summarise_if(is.numeric, funs(list(tidy(t.test(.))))) %>%  # run one-sample t-test for each stimulation condition, return tidy data frames
  unnest() %>% # unpack the list-column with data frame for each test
  kable(.)
```

| stimulation |    estimate|   statistic|    p.value|  parameter|    conf.low|   conf.high| method            | alternative |
|:------------|-----------:|-----------:|----------:|----------:|-----------:|-----------:|:------------------|:------------|
| anodal      |   0.0223864|   0.9712741|  0.3407160|         25|  -0.0250828|   0.0698555| One Sample t-test | two.sided   |
| cathodal    |  -0.0755975|  -3.1803958|  0.0038987|         25|  -0.1245524|  -0.0266426| One Sample t-test | two.sided   |

Follow-up Bayesian one-sample t-test:

``` r
meanAccCenterStats %>%
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(deviation = mean(deviation.end.baseline)) %>% # average over all other variables
  spread(stimulation, deviation) %>% # make separate columns with test data
  summarise_if(is.numeric, funs(extractBF(ttestBF(.), onlybf = TRUE))) %>% # run Bayesian t-test on each column, keeping only the BF
  gather(stimulation,BF,anodal,cathodal) %>% # make row for each stimulation condition
  kable(.)
```

| stimulation |          BF|
|:------------|-----------:|
| anodal      |   0.3173718|
| cathodal    |  10.5264920|

Saccade endpoint variability
----------------------------

*See the `accuracy.nb.html` notebook for more explanation and exploration of the saccade endpoint variability data*

Calculate endpoint variability, as standard deviation of the horizontal component (x-coordinate) of saccades

``` r
varMeanLeg <- preproc %>%
  group_by(subject,stimulation,leg,direction,type) %>% 
  summarise(std.deviation.x = sd(deviation.end.x)) %>% # standard deviation
  ungroup() %>%
  mutate(leg = as.character(leg), # edit leg factor to match other data frames
         leg = replace(leg, leg == "pre", "baseline"),
         leg = factor(leg, levels = c("baseline", "tDCS", "post.1", "post.2"))
         )
```

Make a separate data frame where the mean variability is subtracted from subsequent blocks

``` r
varMeanBaseline <- varMeanLeg %>%
  group_by(subject,stimulation,direction,type) %>% # for each condition, subtract baseline scores and make new columns
  summarise(baseline = std.deviation.x[leg == "baseline"] - std.deviation.x[leg == "baseline"],
           tDCS = std.deviation.x[leg == "tDCS"] - std.deviation.x[leg == "baseline"],
           post.1 = std.deviation.x[leg == "post.1"] - std.deviation.x[leg == "baseline"],
           post.2 = std.deviation.x[leg == "post.2"] - std.deviation.x[leg == "baseline"]) %>%
  gather(leg, std.deviation.x.baseline, baseline, tDCS, post.1, post.2)  %>% # gather new columns to use as factor 
  mutate(leg = factor(leg, levels = c("baseline", "tDCS", "post.1", "post.2"))) # reorder factor levels
```

### Figure 7

``` r
baselineLabelVar <- varMeanLeg %>%
  group_by(stimulation,direction,type,leg) %>%
  filter(leg == "baseline") %>%
  summarise(std.deviation.x.baseline = round(mean(std.deviation.x), digits = 2))
```

Plot mean *lateral saccade* variability per participant

``` r
plot_sub_var_lateral <- ggplot(filter(varMeanLeg, type == "lateral"), aes(leg, std.deviation.x)) +
  facet_grid(stimulation ~ direction) +
  geom_line(aes(group = subject, color = subject), alpha = 0.5) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), linetype = direction), geom = "line") +
  stat_summary(fun.y = mean, geom = "point", aes(shape = stimulation)) +
  guides(colour = FALSE) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous(expression("standard deviation " ( degree)), breaks = seq(0.5,2,0.5), limits = c(0,2)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank(), legend.position = "none") +
  ggtitle("Lateral saccades")
```

Plot mean *center saccade* variability per participant

``` r
plot_sub_var_center <- ggplot(filter(varMeanLeg, type == "center"), aes(leg, std.deviation.x)) +
  facet_grid(stimulation ~ direction) +
  geom_line(aes(group = subject, color = subject), alpha = 0.5) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), linetype = direction), geom = "line") +
  stat_summary(fun.y = mean, geom = "point", aes(shape = stimulation)) +
  guides(colour = FALSE) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous(expression("standard deviation " ( degree)), breaks = seq(0.5,2,0.5), limits = c(0,2)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank(), legend.position = "none") +
  ggtitle("Center saccades")
```

Plot group mean *lateral saccade* variability as change from baseline

``` r
plot_group_var_lateral <- ggplot(filter(varMeanBaseline, type == "lateral"), aes(leg, std.deviation.x.baseline)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = base_line_size) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "line", position = position_dodge(width = 0.5), size = .75) + 
  stat_summary(fun.data = mean_cl_normal, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "linerange", position = position_dodge(width = 0.5), show.legend = FALSE) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, shape = stimulation), geom = "point", position = position_dodge(width = 0.5), size = 2) +
  geom_text(data = subset(baselineLabelVar, type=="lateral"), aes(y = c(.05,.025,.1,.075), label=(std.deviation.x.baseline), color = stimulation), position = position_dodge(width = 0.5), size = mm_to_pt) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_colour_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_y_continuous(expression("SD change from baseline " ( degree)), breaks = seq(-.2,.2,.1)) +
  coord_cartesian(ylim = c(-.2,.2)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank())
```

Plot group mean *center saccade* variability as change from baseline

``` r
plot_group_var_center <- ggplot(filter(varMeanBaseline, type == "center"), aes(leg, std.deviation.x.baseline)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = base_line_size) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "line", position = position_dodge(width = 0.5), size = .75) + 
  stat_summary(fun.data = mean_cl_normal, aes(group = interaction(stimulation, direction), color = stimulation, linetype = direction), geom = "linerange", position = position_dodge(width = 0.5), show.legend = FALSE) +
  stat_summary(fun.y = mean, aes(group = interaction(stimulation, direction), color = stimulation, shape = stimulation), geom = "point", position = position_dodge(width = 0.5), size = 2) +
  geom_text(data = subset(baselineLabelVar, type=="center"), aes(y = c(.05,.025,.1,.075), label=(std.deviation.x.baseline), color = stimulation), position = position_dodge(width = 0.5), size = mm_to_pt) +
  scale_colour_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_discrete("time", labels = c("baseline", "tDCS", "post-1", "post-2")) +
  scale_y_continuous(expression("SD change from baseline " ( degree)), breaks = seq(-.2,.2,.1)) +
  coord_cartesian(ylim = c(-.2,.2)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5), axis.title.x = element_blank())
```

Combine all the plots:

``` r
legend_fig7 <- get_legend(plot_group_acc_lateral)
figure_7 <- plot_grid(
  plot_sub_var_lateral, plot_sub_var_center,
  plot_group_var_lateral + theme(legend.position = "none"),
  plot_group_var_center + theme(legend.position = "none"), rel_heights = c(1,.75))
figure_7 <- plot_grid(figure_7,legend_fig7,rel_widths = c(1,1/10))
figure_7
```

<img src="sacc-tDCS_outline_files/figure-markdown_github/Figure 7-1.png" width="75%" style="display: block; margin: auto;" />

Save the plot:

``` r
ggsave("fig/figure_7.pdf", plot = figure_7, width = 180, height = 112.5, units = "mm")
ggsave("fig/figure_7.png", plot = figure_7, width = 180, height = 112.5, units = "mm")
```

### Statistics

#### Baseline tests

``` r
varMeanLeg %>%
  filter(leg == "baseline") %>%
  group_by(direction,type) %>% 
  nest() %>% 
  mutate(stats = map(data, ~t.test(formula = std.deviation.x~stimulation, paired = TRUE, data =.))) %>% # run t-test on the data frames
  mutate(tidy_model = map(stats, tidy)) %>%
  unnest(tidy_model, .drop = TRUE) %>% 
  kable(.)
```

| direction | type    |    estimate|   statistic|    p.value|  parameter|    conf.low|  conf.high| method        | alternative |
|:----------|:--------|-----------:|-----------:|----------:|----------:|-----------:|----------:|:--------------|:------------|
| left      | lateral |  -0.0046841|  -0.0964036|  0.9239687|         25|  -0.1047541|  0.0953859| Paired t-test | two.sided   |
| left      | center  |  -0.0436919|  -0.9898113|  0.3317503|         25|  -0.1346033|  0.0472195| Paired t-test | two.sided   |
| right     | lateral |  -0.0130790|  -0.2683130|  0.7906600|         25|  -0.1134716|  0.0873137| Paired t-test | two.sided   |
| right     | center  |  -0.0597899|  -1.3048234|  0.2038368|         25|  -0.1541625|  0.0345827| Paired t-test | two.sided   |

#### Lateral saccades

##### Classical ANOVA

Repeated measures ANOVA matching [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045).

**Data**: Baseline subtracted, lateral saccades

**Dependent measure**: saccade endpoint variability (x-coordinate standard deviation)

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

Prepare data:

``` r
meanVarLateralStats <- varMeanBaseline %>%
  ungroup() %>%
  filter(leg != "baseline", type == "lateral") %>%
  select(-type) %>%
  mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) %>%
  mutate(subject = factor(subject))
```

``` r
aovMeanVarLateral <- ezANOVA(data = data.frame(meanVarLateralStats), dv = std.deviation.x.baseline, wid = subject, within = .(stimulation,leg,direction), type = 3)

kable(aovMeanVarLateral$ANOVA)
```

|     | Effect                    |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  1.2189602|  0.2800786|          |  0.0138324|
| 3   | leg                       |    2|   50|  0.4518157|  0.6390443|          |  0.0011604|
| 4   | direction                 |    1|   25|  0.2471961|  0.6234009|          |  0.0010690|
| 5   | stimulation:leg           |    2|   50|  1.1183082|  0.3348684|          |  0.0030580|
| 6   | stimulation:direction     |    1|   25|  0.1204196|  0.7314841|          |  0.0004919|
| 7   | leg:direction             |    2|   50|  0.9403214|  0.3972980|          |  0.0023341|
| 8   | stimulation:leg:direction |    2|   50|  0.1873223|  0.8297557|          |  0.0003276|

``` r
kable(aovMeanVarLateral$`Mauchly's Test for Sphericity`)
```

|     | Effect                    |         W|          p| p&lt;.05 |
|-----|:--------------------------|---------:|----------:|:---------|
| 3   | leg                       |  0.994425|  0.9351138|          |
| 5   | stimulation:leg           |  0.994966|  0.9412372|          |
| 7   | leg:direction             |  0.908637|  0.3167270|          |
| 8   | stimulation:leg:direction |  0.903796|  0.2970605|          |

``` r
kable(aovMeanVarLateral$`Sphericity Corrections`)
```

|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.9944559|  0.6379765|                |  1.0801686|  0.6390443|                |
| 5   | stimulation:leg           |  0.9949912|  0.3346792|                |  1.0808237|  0.3348684|                |
| 7   | leg:direction             |  0.9162854|  0.3907925|                |  0.9851512|  0.3961957|                |
| 8   | stimulation:leg:direction |  0.9122389|  0.8102650|                |  0.9802676|  0.8255939|                |

##### Bayesian ANOVA

Same design as the classical ANOVA

Bayes factors agains the null model:

``` r
bfMeanVarLateral = anovaBF(std.deviation.x.baseline~stimulation*leg*direction+subject, data = data.frame(meanVarLateralStats), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) 
bfMeanVarLateral = sort(bfMeanVarLateral, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfMeanVarLateral), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |         bf|
|-------------------------------------------------------------------------------------------------------------------------------|----------:|
| stimulation + subject                                                                                                         |  1.6065491|
| stimulation + direction + subject                                                                                             |  0.2426688|
| direction + subject                                                                                                           |  0.1495108|
| stimulation + leg + subject                                                                                                   |  0.0688596|
| stimulation + direction + stimulation:direction + subject                                                                     |  0.0448681|
| leg + subject                                                                                                                 |  0.0427400|
| stimulation + direction + leg + subject                                                                                       |  0.0101775|
| stimulation + leg + stimulation:leg + subject                                                                                 |  0.0072123|
| direction + leg + subject                                                                                                     |  0.0063480|
| stimulation + direction + stimulation:direction + leg + subject                                                               |  0.0020752|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |  0.0011029|
| stimulation + direction + leg + direction:leg + subject                                                                       |  0.0009751|
| direction + leg + direction:leg + subject                                                                                     |  0.0005826|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |  0.0002069|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |  0.0001855|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |  0.0001048|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |  0.0000190|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |  0.0000031|

Inclusion Bayes factor across matched models:

``` r
kable(inclusionBF(bfMeanVarLateral, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |     1.6087893|
| direction                 |     0.1503684|
| stimulation:direction     |     0.1856838|
| leg                       |     0.0427784|
| stimulation:leg           |     0.1050887|
| direction:leg             |     0.0937721|
| stimulation:direction:leg |     0.1642401|

#### Center saccades

##### Classical ANOVA

Repeated measures ANOVA matching [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045).

**Data**: Baseline subtracted, center saccades

**Dependent measure**: saccade endpoint variability (x-coordinate standard deviation)

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

Prepare data:

``` r
meanVarCenterStats <- varMeanBaseline %>%
  ungroup() %>%
  filter(leg != "baseline", type == "center") %>%
  select(-type) %>%
  mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) %>%
  mutate(subject = factor(subject))
```

``` r
aovMeanVarCenter <- ezANOVA(data = data.frame(meanVarCenterStats), dv = std.deviation.x.baseline, wid = subject, within = .(stimulation,leg,direction), type = 3)

kable(aovMeanVarCenter$ANOVA)
```

|     | Effect                    |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  3.8851203|  0.0598820|          |  0.0399400|
| 3   | leg                       |    2|   50|  2.7073528|  0.0764935|          |  0.0075993|
| 4   | direction                 |    1|   25|  0.6021177|  0.4450504|          |  0.0033319|
| 5   | stimulation:leg           |    2|   50|  1.1805422|  0.3155254|          |  0.0036640|
| 6   | stimulation:direction     |    1|   25|  0.1745759|  0.6796436|          |  0.0010192|
| 7   | leg:direction             |    2|   50|  2.1601871|  0.1259451|          |  0.0035285|
| 8   | stimulation:leg:direction |    2|   50|  0.4735067|  0.6255786|          |  0.0007475|

``` r
kable(aovMeanVarCenter$`Mauchly's Test for Sphericity`)
```

|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.8960127|  0.2677749|          |
| 5   | stimulation:leg           |  0.8569616|  0.1568686|          |
| 7   | leg:direction             |  0.6679634|  0.0078892| \*       |
| 8   | stimulation:leg:direction |  0.8089994|  0.0785921|          |

``` r
kable(aovMeanVarCenter$`Sphericity Corrections`)
```

|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.9058075|  0.0824883|                |  0.9725125|  0.0781996|                |
| 5   | stimulation:leg           |  0.8748612|  0.3116036|                |  0.9353175|  0.3136546|                |
| 7   | leg:direction             |  0.7507301|  0.1405427|                |  0.7880908|  0.1383445|                |
| 8   | stimulation:leg:direction |  0.8396302|  0.5929099|                |  0.8932129|  0.6044505|                |

##### Bayesian ANOVA

Same design as the classical ANOVA

Bayes factors agains the null model:

``` r
bfMeanVarCenter = anovaBF(std.deviation.x.baseline~stimulation*leg*direction+subject, data = data.frame(meanVarCenterStats), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) 
bfMeanVarCenter = sort(bfMeanVarCenter, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfMeanVarCenter), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |           bf|
|-------------------------------------------------------------------------------------------------------------------------------|------------:|
| stimulation + subject                                                                                                         |  141.5596032|
| stimulation + direction + subject                                                                                             |   31.0170327|
| stimulation + leg + subject                                                                                                   |   17.3145217|
| stimulation + direction + stimulation:direction + subject                                                                     |    6.2783439|
| stimulation + direction + leg + subject                                                                                       |    3.8046477|
| stimulation + leg + stimulation:leg + subject                                                                                 |    1.9607544|
| stimulation + direction + stimulation:direction + leg + subject                                                               |    0.7921655|
| stimulation + direction + leg + direction:leg + subject                                                                       |    0.4464083|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |    0.4353829|
| direction + subject                                                                                                           |    0.2145015|
| leg + subject                                                                                                                 |    0.1151435|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |    0.0932772|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |    0.0841059|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |    0.0474962|
| direction + leg + subject                                                                                                     |    0.0246581|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |    0.0097800|
| direction + leg + direction:leg + subject                                                                                     |    0.0032830|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |    0.0011420|

Inclusion Bayes factor across matched models:

``` r
kable(inclusionBF(bfMeanVarCenter, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |   143.0054521|
| direction                 |     0.2191801|
| stimulation:direction     |     0.2030063|
| leg                       |     0.1224590|
| stimulation:leg           |     0.1134795|
| direction:leg             |     0.1147686|
| stimulation:direction:leg |     0.1167715|

##### Follow-up tests

###### Main effect of stimulation

Follow-up One-sample t-test:

``` r
meanVarCenterStats %>%
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(std.deviation.x. = mean(std.deviation.x.baseline)) %>% # average over all other variables (df is now still grouped per stimulation)
  summarise_if(is.numeric, funs(list(tidy(t.test(.))))) %>%  # run one-sample t-test for each stimulation condition, return tidy data frames
  unnest() %>% # unpack the list-column with data frame for each test
  kable(.)
```

| stimulation |    estimate|  statistic|    p.value|  parameter|    conf.low|  conf.high| method            | alternative |
|:------------|-----------:|----------:|----------:|----------:|-----------:|----------:|:------------------|:------------|
| anodal      |   0.0477794|   1.683151|  0.1047947|         25|  -0.0106845|  0.1062434| One Sample t-test | two.sided   |
| cathodal    |  -0.0340719|  -1.236236|  0.2278603|         25|  -0.0908350|  0.0226911| One Sample t-test | two.sided   |

Follow-up Bayesian one-sample t-test:

``` r
meanVarCenterStats %>%
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(std.deviation.x = mean(std.deviation.x.baseline)) %>% # average over all other variables
  spread(stimulation, std.deviation.x) %>% # make separate columns with test data
  summarise_if(is.numeric, funs(extractBF(ttestBF(.), onlybf = TRUE))) %>% # run Bayesian t-test on each column, keeping only the BF
  gather(stimulation,BF,anodal,cathodal) %>% # make row for each stimulation condition
  kable(.)
```

| stimulation |         BF|
|:------------|----------:|
| anodal      |  0.7138830|
| cathodal    |  0.4103342|

Quantile analysis
=================

*See the `RT_quantiles.nb.html` notebook for more explanation and exploration of the saccade latency distribution data*

Estimate quantiles and compute shift function:

``` r
qData <- groupData %>%
  filter(!is.na(latency)) %>% # discard missing saccades
  group_by(subject,leg,type,direction) %>% # for each condition
  nest(stimulation,latency) %>% # make a list_column "data" out of the stimulation and latency columns. Now, each group has its own data frame consisting of the stimulation and latency columns for that group
  mutate(shift = purrr::map(data, ~ shiftdhd_pbci(., formula = latency ~ stimulation, nboot = 2000))) # for each group, estimate quantiles and compute shift function, and store as a list column in "shift"

qData <- unnest(qData,shift) # unpack the list column to get the model results for each group in the original data frame

# Alternatively, instead of "nest ...", "mutate ..." and "unnest ...", you can call shiftdhd through dplyr::do (but this is basically depreated in favor of purrr:map)
# do(shiftdhd(.[,c("stimulation","latency")], formula = latency ~ stimulation, nboot = 100)) # estimate quantiles and compute shift function
```

Note that this takes quite a while to compute with a large number of bootstrap samples, so we will load the result from disk in the next code chunk. By default, the chunk above will not run becuase of the `eval=FALSE` statement; remove this to execute it and compute the result from scratch.

``` r
qData <- read_csv(here("data", "sacc-tDCS_quantiles.csv")) %>%
  filter(!(subject %in% subs2exclude)) %>% # exclude subjects
    mutate(leg = replace(leg, leg == "pre", "baseline"), # rename and reoder levels
         leg = replace(leg, leg == "post.1", "post-1"),
         leg = replace(leg, leg == "post.2", "post-2"),
         leg = factor(leg, levels = c("baseline", "tDCS", "post-1", "post-2")))
```

``` r
qStats <- qData %>%
  group_by(subject,leg,type,direction) %>%
  mutate(deco = c(seq(1,5),seq(4,1))) %>% # add code of deciles to data frame
  group_by(leg,type,direction,q) %>% 
  mutate(anodal = mean(anodal)) %>% # mean of "anodal" quantiles OVER subjects
  group_by(leg,type,direction) %>% 
  mutate(anodal_median = median(anodal)) # median for plotting
```

For each quantile, count which subjects show significant effects, and in which direction

``` r
qSig <- qData %>%
  group_by(subject,leg,type,direction,q) %>%
  mutate(significance = NA,
         # If (anodal - cathodal) quantile difference is positive, latency for anodal > latency for cathodal
         significance = replace(significance, p_value < p_crit & difference > 0, "cathodal.faster"),
         # If (anodal - cathodal) quantile difference is negative, latency for anodal < latency for cathodal
         significance = replace(significance, p_value < p_crit & difference < 0, "anodal.faster")
  )
```

Figure 5
--------

Shift functions for *lateral saccades*

``` r
plot_shift_lateral <- ggplot(filter(qStats, type == "lateral"), aes(anodal, difference)) +
   facet_grid(leg ~ direction) +
   geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
   geom_vline(aes(xintercept = anodal_median), linetype = "dashed", alpha = 0.5) +
   stat_summary(fun.data = mean_cl_normal, geom = "linerange", colour = "black", size = 0.5) +
   stat_summary(fun.y = mean, geom = "line", size = 1, colour = "grey50", alpha = 0.5) +
   stat_summary(fun.y = mean, geom = "point", aes(fill = deco), size = 2, colour = "black", shape = 21) +
   scale_fill_gradient(low = "white", high = "grey30", guide = FALSE) +
   scale_x_continuous("anodal deciles (ms)", limits = c(80,200), breaks = seq(100,200,25)) +
   scale_y_continuous("anodal - cathodal deciles (ms)")
```

Number of subjects with signifcant difference per quantile for lateral saccades

``` r
plot_sig_lateral <- ggplot(filter(qSig, type == "lateral", !is.na(significance)), aes(q, fill = significance)) +
  facet_grid(leg ~ direction) +
  geom_bar(position = "stack") +
  stat_bin(binwidth = .1, geom = "text", size = mm_to_pt, aes(label = ..count..), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_continuous("decile", breaks = seq(0.1,0.9,0.1), labels = seq(1,9,1)) +
  scale_y_continuous("number of participants with significant difference", limits = c(0, length(unique(qSig$subject))), breaks = c(0,10,20,length(unique(qSig$subject))))
```

Shift functions for *center saccades*

``` r
plot_shift_center <- ggplot(filter(qStats, type == "center"), aes(anodal, difference)) +
   facet_grid(leg ~ direction) +
   geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
   geom_vline(aes(xintercept = anodal_median), linetype = "dashed", alpha = 0.5) +
   stat_summary(fun.data = mean_cl_normal, geom = "linerange", colour = "black", size = 0.5) +
   stat_summary(fun.y = mean, geom = "line", size = 1, colour = "grey50", alpha = 0.5) +
   stat_summary(fun.y = mean, geom = "point", aes(fill = deco), size = 2, colour = "black", shape = 21) +
   scale_fill_gradient(low = "white", high = "grey30", guide = FALSE) +
   scale_x_continuous("anodal deciles (ms)", limits = c(80,200), breaks = seq(100,200,25)) +
   scale_y_continuous("anodal - cathodal deciles (ms)")
```

Number of subjects with signifcant difference per quantile for lateral saccades

``` r
plot_sig_center <- ggplot(filter(qSig, type == "center", !is.na(significance)), aes(q, fill = significance)) +
  facet_grid(leg ~ direction) +
  geom_bar(position = "stack") +
  stat_bin(binwidth = .1, geom = "text", size = mm_to_pt, aes(label = ..count..), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_continuous("decile", breaks = seq(0.1,0.9,0.1), labels = seq(1,9,1)) +
  scale_y_continuous("number of participants with significant difference", limits = c(0, length(unique(qSig$subject))), breaks = c(0,10,20,length(unique(qSig$subject))))
```

``` r
figure_5_lateral <- plot_grid(plot_shift_lateral, plot_sig_lateral + theme(legend.position = "none"), rel_widths = c(1.33,1)) # top row
figure_5_center <- plot_grid(plot_shift_center,plot_sig_center + theme(legend.position = "none"), rel_widths = c(1.33,1)) # bottom row

#title objects
title_lateral <- ggdraw() + draw_label("Lateral saccades", fontface = 'bold', size = base_font_size)
title_center <- ggdraw() + draw_label("Center saccades", fontface = 'bold', size = base_font_size)

#combine with titles
figure_5_lateral <- plot_grid(title_lateral, figure_5_lateral, ncol = 1, rel_heights = c(0.05, 1))
figure_5_center <- plot_grid(title_center, figure_5_center, ncol = 1, rel_heights = c(0.05, 1))
figure_5 <- plot_grid(figure_5_lateral, figure_5_center, nrow = 2)
figure_5
```

<img src="sacc-tDCS_outline_files/figure-markdown_github/Figure 5-1.png" width="75%" style="display: block; margin: auto;" />

Save the plot:

``` r
ggsave("fig/figure_5.pdf", plot = figure_5, width = 180, height = 180, units = "mm")
ggsave("fig/figure_5.png", plot = figure_5, width = 180, height = 180 , units = "mm")
```

Supplementary results
=====================

tDCS adverse effect questionnaire
---------------------------------

*See the `questionnaires.nb.html` notebook for more explanation and exploration of the questionnaire data*

``` r
# Load the data frame
dataFile <- here("data", "tdcs_sensations.csv")
sensData <- read_csv2(dataFile, col_types = cols(
  session = col_factor(c("first","second")),
  stimulation = col_factor(c("anodal","cathodal"))
))
kable(head(sensData)) # show data frame
```

| subject | session | stimulation |  itching|  tingling|  burning|  pain|  headache|  fatigue|  dizziness|  nausea|  conf.itching|  conf.tingling|  conf.burning|  conf.pain|  conf.headache|  conf.fatigue|  conf.dizziness|  conf.nausea| felt.more |
|:--------|:--------|:------------|--------:|---------:|--------:|-----:|---------:|--------:|----------:|-------:|-------------:|--------------:|-------------:|----------:|--------------:|-------------:|---------------:|------------:|:----------|
| S01     | first   | cathodal    |        2|         4|        4|     2|         0|        0|          0|       0|             3|              4|             4|          4|              0|             0|               0|            0| cathode   |
| S01     | second  | anodal      |        1|         1|        1|     0|         0|        1|          0|       0|             4|              4|             4|          0|              0|             1|               0|            0| anode     |
| S02     | first   | anodal      |        1|         2|        2|     0|         0|        1|          0|       0|             4|              4|             4|          0|              0|             1|               0|            0| anode     |
| S02     | second  | cathodal    |        1|         2|        2|     0|         0|        0|          0|       0|             4|              4|             4|          0|              0|             0|               0|            0| equal     |
| S03     | first   | cathodal    |        0|         0|        0|     0|         0|        0|          0|       0|             0|              0|             0|          0|              0|             0|               0|            0| equal     |
| S03     | second  | anodal      |        0|         0|        0|     0|         0|        0|          0|       0|             0|              0|             0|          0|              0|             0|               0|            0| equal     |

Participants were asked to which degree the following sensations were present during stimulation: *tingling*, *itching sensation*, *burning sensation*, *pain*, *headache*, *fatigue*, *dizziness* and *nausea*. Each was rated on a scale from 0-4:

1.  none
2.  a little
3.  moderate
4.  strong
5.  very strong

They also rated their confidence *that the sensations were caused by the stimulation* on a scale from 0-4 (columns starting with `conf.`):

1.  n/a (meaning they rated the sensation a 0 on the previous scale)
2.  unlikely
3.  possibly
4.  likely
5.  very likely

**Factors**:

-   *subject*: subject ID (`S01`, `S02`, etc)
-   *session*: Whether data are from the `first` or `second` session
-   *stimulation*: Whether data are from the `anodal` or `cathodal` session

Calculate how many anodal and cathodal sessions were rated:

``` r
idxComplete <- rowSums(is.na(sensData)) != ncol(sensData) - 3 # rows that do not have all NAs (except the 3 factor columns)
# calculate number of questionnaires completed per stimulation type
nAnodal <- sum(sensData$stimulation == "anodal" & idxComplete) 
nCathodal <- sum(sensData$stimulation == "cathodal" & idxComplete)
```

### Figure S1

Prepare the sensation intensity and confidence ratings.

``` r
# Make long form data frame of sensation intensity
ratings <- sensData %>%
  select(everything(), -contains("conf"), -felt.more) %>% # drop other columns
  gather(sensation, rating, itching:nausea) # make long form

# Make long form data frame of sensation confidence
confidence <- sensData %>%
  select(contains("conf"), subject, session, stimulation) %>% 
  gather(sensation, confidence, conf.itching:conf.nausea) %>%
  mutate(sensation = str_replace(sensation, "conf.", "")) # get rid of "conf." prefix so it matches the sensation intensity tibble

# Join the two data frames
sensDataLong <- dplyr::full_join(ratings,confidence)
```

Make plot of sensation distribution:

``` r
plot_sensation <- ggplot(sensDataLong, aes(rating, fill = stimulation)) +
    facet_wrap(~sensation, nrow = 1) +
    geom_bar(position = "stack") +
    stat_bin(binwidth = 1, geom = "text", size = mm_to_pt, aes(label = ..count..), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#F25F5C", "#4B93B1")) +
    xlim(0.5,4.5) + # exclude "0" ratings that were not present
    ylim(0,35) + # bound plot at max number of ratings
    ylab("number of sessions") +
  ggtitle("Sensation intensity")
```

Make plot of confidence distribution:

``` r
plot_confidence <- ggplot(sensDataLong, aes(confidence, fill = stimulation)) +
    facet_wrap(~sensation, nrow = 1) +
    geom_bar(position = "stack") +
    stat_bin(binwidth = 1, geom = "text", size = mm_to_pt, aes(label = ..count..), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#F25F5C", "#4B93B1")) +
    xlim(0.5,4.5) +
    ylim(0,35) +
    xlab("ratings") +
    ylab("number of sessions") +
    ggtitle("Confidence")
```

Put together and show plot:

``` r
legend_S1 <- get_legend(plot_sensation)
figure_S1 <- plot_grid(
  plot_sensation + theme(legend.position = "none"),
  plot_confidence + theme(legend.position = "none"),
  nrow = 2)
figure_S1 <- plot_grid(figure_S1,legend_S1, rel_widths = c(1,1/10))
figure_S1
```

<img src="sacc-tDCS_outline_files/figure-markdown_github/tDCS AE plot-1.png" width="75%" style="display: block; margin: auto;" />

Save the plot:

``` r
ggsave("fig/figure_S1.pdf", plot = figure_S1, width = 180, height = 112.5, units = "mm")
ggsave("fig/figure_S1.png", plot = figure_S1, width = 180, height = 112.5, units = "mm")
```

### Statistics

Mann-Whitney U tests for difference in sensation intensity ratings between anodal and cathodal sessions.

``` r
sensationList <- c("itching", "tingling", "burning", "pain", "headache", "fatigue", "dizziness", "nausea")
senseTests <- data.frame(sensation = sensationList, p.value = NA) # initialize results data frame
for (item in sensationList) {
  testData <- sensData[[item]] # extract column with test dat
  tmp <- wilcox.test(testData[sensData$stimulation == "anodal"], testData[sensData$stimulation == "cathodal"])
  senseTests$p.value[senseTests$sensation %in% item] <- tmp$p.value # put p-value in row of results data frame
}
kable(senseTests)
```

| sensation |    p.value|
|:----------|----------:|
| itching   |  0.5787360|
| tingling  |  0.3652154|
| burning   |  0.1913439|
| pain      |  0.9287263|
| headache  |  0.9648458|
| fatigue   |  0.4330014|
| dizziness |  0.7203021|
| nausea    |        NaN|

Note that the p-value for "nausea" is undefined because it was never reported (i.e. all ranks are tied at 0).

Frontal eye field coordinates
-----------------------------

*See the `frontal_eye_field.nb.html` notebook for more explanation and exploration of the frontal eye field coordinate data*

### Native space

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

### MNI space (Table S1)

Load in the coordinates that were transformed to MNI space.

#### Table S1

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

#### Descriptives

Calculate descriptive statistics over subjects:

``` r
mniCoords %>%
  gather(dimension, coord, MNI_X:MNI_Z) %>%
  group_by(dimension) %>%
  summarise_at(vars(coord), funs(mean, min, max, sd)) %>% 
  kable(.)
```

| dimension |       mean|       min|       max|        sd|
|:----------|----------:|---------:|---------:|---------:|
| MNI\_X    |  31.470735|   25.7061|  41.21270|  4.036810|
| MNI\_Y    |  -1.806643|  -10.0611|   8.39057|  3.931623|
| MNI\_Z    |  51.645608|   43.4371|  58.96710|  3.922225|
