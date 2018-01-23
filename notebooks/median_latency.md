sacc-tDCS: Median saccade latency
================
Leon Reteig

-   [Load data](#load-data)
    -   [Eye data](#eye-data)
    -   [Subject metadata](#subject-metadata)
-   [Inspect median data](#inspect-median-data)
    -   [Remove outliers](#remove-outliers)
    -   [Data per block](#data-per-block)
    -   [15-minute intervals (collapse 3 blocks)](#minute-intervals-collapse-3-blocks)
-   [Median reaction time plots and analyses](#median-reaction-time-plots-and-analyses)
    -   [Line plot per leg over all subjects](#line-plot-per-leg-over-all-subjects)
    -   [Subtract baseline](#subtract-baseline)
-   [Statistics](#statistics)
    -   [Frequentist](#frequentist)
    -   [Bayesian](#bayesian)

R notebook for analyses of median saccade latency in the `sacc-tDCS` dataset. Previous processing:

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
library(ez) # ANOVA
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
library(broom) # transform model output into a data frame
library(knitr) # R markdown output (html, pdf, etc.)
# set default output and figure options
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 7, fig.asp = 0.618, out.width = "75%", fig.align = "center")


source(here("src", "lib", "InclusionBF.R"))

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
    ##  [1] knitr_1.15.1         broom_0.4.2          BayesFactor_0.9.12-2
    ##  [4] Matrix_1.2-9         coda_0.19-1          ez_4.4-0            
    ##  [7] dplyr_0.5.0          purrr_0.2.2          readr_1.1.0         
    ## [10] tidyr_0.6.1          tibble_1.3.0         ggplot2_2.2.1       
    ## [13] tidyverse_1.1.1      here_0.1            
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtools_3.5.0       pbapply_1.3-2      reshape2_1.4.2    
    ##  [4] splines_3.4.0      haven_1.0.0        lattice_0.20-35   
    ##  [7] colorspace_1.3-2   htmltools_0.3.6    yaml_2.1.14       
    ## [10] mgcv_1.8-17        nloptr_1.0.4       foreign_0.8-67    
    ## [13] DBI_0.6-1          modelr_0.1.0       readxl_1.0.0      
    ## [16] plyr_1.8.4         stringr_1.2.0      MatrixModels_0.4-1
    ## [19] munsell_0.4.3      gtable_0.2.0       cellranger_1.1.0  
    ## [22] rvest_0.3.2        mvtnorm_1.0-6      psych_1.7.3.21    
    ## [25] evaluate_0.10      forcats_0.2.0      SparseM_1.77      
    ## [28] quantreg_5.33      pbkrtest_0.4-7     parallel_3.4.0    
    ## [31] Rcpp_0.12.10       scales_0.4.1       backports_1.0.5   
    ## [34] jsonlite_1.4       lme4_1.1-13        mnormt_1.5-5      
    ## [37] hms_0.3            digest_0.6.12      stringi_1.1.5     
    ## [40] grid_3.4.0         rprojroot_1.2      tools_3.4.0       
    ## [43] magrittr_1.5       lazyeval_0.2.0     car_2.1-4         
    ## [46] MASS_7.3-47        xml2_1.1.1         lubridate_1.6.0   
    ## [49] assertthat_0.2.0   minqa_1.2.4        rmarkdown_1.5     
    ## [52] httr_1.2.1         R6_2.2.0           nnet_7.3-12       
    ## [55] nlme_3.1-131       compiler_3.4.0

Load data
=========

Eye data
--------

The .csv file with the eye tracking data was created in MATLAB.

``` r
# Load the data frame
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

``` r
# Load eye tracking data into data frame
dataFile <- here("data", "subject_info.csv")
subjectData <- read_csv2(dataFile, col_names = TRUE, progress = FALSE, col_types = cols(
  session.order = col_factor(c("first.anodal", "first.cathodal"))
))
```

Subject metadata
----------------

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

The main use is to see if the nuisance factor *session.order* covaries with the factors of interest in the design. This could indicate the presence of carryover effects between the stimulation, or a difference in subgroups within the sample (see <http://www.jerrydallal.com/lhsp/crossovr.htm> for an introduction to these kinds of analyses.).

Inspect median data
===================

Remove outliers
---------------

``` r
tooFast <- 50
tooSlow <- 400
badFix <- 1.8
badSacc <- 8
subs2exclude <- c("S28","S16","S22","S21","S25")
```

-   S21 and S25 were tested &lt; 48h apart
-   S16, S22 and S28 had fewer than 50 saccades per condition after trial rejection

Criteria for outlier saccades:

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
# Remove outliers and subjects
groupData <- filter(groupData,
                    # outliers
                    latency >= tooFast,
                    latency <= tooSlow,
                    deviation.start <= badFix,
                    deviation.end.x <= badSacc,
                    # subjects
                    !(subject %in% subs2exclude),
                    # missing values
                    complete.cases(groupData)
)
```

Data per block
--------------

``` r
# Compute median in each condition
latencyMedian <- groupData %>%
  group_by(subject,stimulation,leg,block,type,direction) %>%
  summarise(latency = median(latency))
```

### Full factorial plot

``` r
# Plot out all the data
fullPlot <- ggplot(latencyMedian, aes(interaction(block,leg), latency, color = stimulation, shape = stimulation)) +
  facet_grid(type ~ direction) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5))
fullPlot
```

<img src="median_latency_files/figure-markdown_github/Full factorial plot-1.png" width="75%" style="display: block; margin: auto;" />

For some reason cathodal is always a little faster than anodal. There is probably a few ms difference between them on average, also in the baseline block. This is clearly just random variation, but it could be a problem, since the effects we expect are not much bigger...

The above plot also shows a lot of variability, so it might be best to average over each 3 consecutive blocks so the data come in 15-minute intervals.

15-minute intervals (collapse 3 blocks)
---------------------------------------

``` r
# Compute median per leg
latencyMedianLeg <- groupData %>%
  group_by(subject,stimulation,direction,type) %>% 
  summarise(baseline = median(latency[leg == "pre"]), # take average of 3 blocks, make new column
            tDCS = median(latency[leg == "tDCS"]),
            post.1 = median(latency[leg == "post" & block <= 3]),
            post.2 = median(latency[leg == "post" & block >= 4])) %>%
 gather(leg,latency,baseline,tDCS,post.1,post.2) %>% # gather new columns to use as factor
 mutate(leg = factor(leg, levels = c("baseline", "tDCS", "post.1", "post.2"))) # reorder factor levels
```

### Line plot per leg, individual subjects

Now make the same plot, but for the data collapsed over 3 blocks. Also draw the plot for each individual subject, so we can see which subjects drive the baseline difference, and in which direction the stimulation effect goes for each subject (if there is any).

``` r
kanaiSubsPlot <- ggplot(latencyMedianLeg, aes(leg, latency, color = stimulation, shape = type)) +         
  facet_wrap(~subject, ncol = 5) +
  stat_summary(fun.y = mean, geom = "point", size = 2) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation)) +
  theme(axis.text.x = element_text(angle = 22, vjust = .5))
kanaiSubsPlot
```

<img src="median_latency_files/figure-markdown_github/Line plot for each subject-1.png" width="75%" style="display: block; margin: auto;" />

There are a couple of subjects that show large differences in the baseline already, especially S01. The stimulation effects seem to not be very apparent, but they could also drown out at this scale.

### Compute magnitude of baseline difference

Let's look at the size of the baseline difference per subject.

``` r
baselineDiff <- latencyMedianLeg %>% 
  filter(leg == "baseline") %>% # keep only baseline data
  spread(stimulation, latency) %>% # make separate columns for anodal and cathodal
  mutate(latency.diff = anodal - cathodal) %>% # subtract the difference
  group_by(subject) %>% 
  summarise(latency.diff = mean(latency.diff))# keep the average difference per subject

kable(baselineDiff, caption = 'Difference between baseline saccade latencies in anodal and cathodal session')
```

| subject |  latency.diff|
|:--------|-------------:|
| S01     |        39.125|
| S02     |        18.625|
| S03     |         2.375|
| S04     |        -6.125|
| S05     |       -17.500|
| S06     |        -3.875|
| S07     |         5.625|
| S08     |        11.125|
| S09     |        23.875|
| S10     |        14.750|
| S11     |         1.000|
| S12     |         8.375|
| S13     |         4.000|
| S14     |       -10.375|
| S15     |         1.875|
| S17     |        -5.500|
| S18     |        -2.875|
| S19     |         1.000|
| S20     |         4.750|
| S24     |         3.625|
| S26     |         0.000|
| S27     |         0.500|
| S29     |       -13.625|
| S30     |         7.125|
| S32     |        10.875|
| S33     |        -1.625|

There are a few subjects with substantial latency differences between sessions, so it's probably a good idea to subtract the baseline for the statistical analyses. The mean difference is 3.7 ms. So on average, the cathodal session is a little faster, so the counterbalancing is not perfect, but the difference is not very high.

Median reaction time plots and analyses
=======================================

Here we simply extract median RTs for each condition and use a repeated measures ANOVA for statistical analysis, following [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045).

Line plot per leg over all subjects
-----------------------------------

Let's look at the group average plot for the first time.

``` r
kanaiPlot <- ggplot(latencyMedianLeg, aes(leg, latency, color = stimulation, shape = stimulation)) +         
  facet_grid(type ~ direction) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3)
kanaiPlot
```

<img src="median_latency_files/figure-markdown_github/Line plot per leg-1.png" width="75%" style="display: block; margin: auto;" />

All differences between anodal & cathodal seem to be &lt; 5 ms. Any differences that are there also 1) do not seem to differ much between the time blocks, or 2) be in the opposite direction for anodal and cathodal, with the possible exception of the first post-block in the center-right condition.

### Individual subjects

#### Anodal session

Let's look at the same plot for individual subjects, to see whether there are any consistent patterns or huge outliers there.

``` r
kanaiPlotSubsAnodal <- ggplot(latencyMedianLeg[latencyMedianLeg$stimulation == "anodal", ], aes(leg, latency)) +
  facet_grid(type ~ direction) +
  geom_line(aes(group=subject,color=subject)) +
  stat_summary(fun.y = mean, aes(group = stimulation), geom = "line") +
  stat_summary(fun.y = mean, geom = "point") +
  ggtitle("Anodal session")
kanaiPlotSubsAnodal
```

<img src="median_latency_files/figure-markdown_github/Line plot per subject - anodal-1.png" width="75%" style="display: block; margin: auto;" />

There do appear to be some "outliers", but they are fairly well balanced (below or above the average), and &gt;80% or so of subjects seem to cluster together. I don't see much evidence for systematic individual differences (e.g. slow subjects get faster; fast subjects get slower.)

#### Cathodal session

``` r
kanaiPlotSubsCathodal <- ggplot(latencyMedianLeg[latencyMedianLeg$stimulation == "cathodal", ], aes(leg, latency)) +
  facet_grid(type ~ direction) +
  geom_line(aes(group=subject,color=subject)) +
  stat_summary(fun.y = mean, aes(group = stimulation), geom = "line") +
  stat_summary(fun.y = mean, geom = "point") +
  ggtitle("Cathodal session")
kanaiPlotSubsCathodal
```

<img src="median_latency_files/figure-markdown_github/Line plot per subject - cathodal-1.png" width="75%" style="display: block; margin: auto;" />

There seem to be less clear individual outliers, but the spread also seems a bit bigger (e.g. right-lateral)

Subtract baseline
-----------------

Baseline differences (while small here) are not informative and could obscure real changes from baseline within subjects. Also it's hard to see and compare the magnitude of the effects, also for instance because center saccades are much faster than lateral saccades. Therefore, subtract the baseline from each subsequent measurement.

``` r
latencyMedianBaseline <- latencyMedianLeg %>%
  group_by(subject,stimulation,direction,type) %>% # for each condition, subtract baseline scores and make new columns
  summarise(tDCS = latency[leg == "tDCS"] - latency[leg == "baseline"], 
           post.1 = latency[leg == "post.1"] - latency[leg == "baseline"],
           post.2 = latency[leg == "post.2"] - latency[leg == "baseline"]) %>%
  gather(leg, latency, tDCS, post.1, post.2)  %>% # gather new columns to use as factor 
  mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) # reorder factor levels
```

### Baseline data

Before we look at the change-from-baseline data, it would be good to check the baseline data itself.

#### Baseline reliability

Analyzing change-from-baseline only makes sense if the baseline data are not very noisy. Subtracting a noisy measure from the data will increase the noise in the data, and thereby [reduce power](http://datacolada.org/39). One way to assess this is to look at the correlation between the baseline blocks in both sessions. If there is a high correlation, the pro-saccade task produces comparable latencies each time it is performed, and thus the baseline data should on average not be very noisy.

Let's compute the correlation between all baselines (i.e. for the factors DIRECTION and TYPE).

``` r
baselineCorr <- latencyMedianLeg %>%
  filter(leg == "baseline") %>% # keep only baseline data
  group_by(direction,type) %>% # for each of these 4 pairs
  spread(stimulation,latency) %>% # make 2 columns with baseline latencies of anodal and cathodal session
  nest() %>% # create a separate data frame of all remaining columns. These data frames are stored in a list-column named "data"
  mutate(stats = map(data, ~cor.test(formula = ~ anodal + cathodal, data =.))) %>% # run correlation test on baselines from each condition
  mutate(tidy_model = map(stats, tidy)) %>% # force the four test outputs to also be data frames in a list column ("tidy_model")
  unnest(tidy_model, .drop = TRUE)  # unpack the list-column with results, drop the list-column we used to organize the data
```

Now a scatter plot of the data, with the results of the test:

``` r
latencyMedianLeg %>%
  filter(leg == "baseline") %>%
  spread(stimulation,latency) %>%
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column on session order from other data frame
  ggplot(aes(anodal,cathodal)) +
    facet_grid(type ~ direction) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_smooth(method = "lm") +
    geom_point(aes(color=session.order)) +
    xlim(80,220) + ylim(80,220) +
    geom_text(data = baselineCorr, x = 100, y = 200, aes(label = paste("italic(r) == ", round(estimate,2))), parse = TRUE) +
    labs(title = "Baseline in anodal and cathodal sessions", subtitle = "scatterplot of baseline latencies")
```

<img src="median_latency_files/figure-markdown_github/baseline scatterplots-1.png" width="75%" style="display: block; margin: auto;" />

The correlations are quite high, at least high enough for the baseline subtraction to help power rather than hurt (r's &gt; 0.5).

There are a few outliers though: some subjects deviate quite a bit from the idealized line (dashed).

When the data is split for session order, a clear pattern becomes apparent: latencies are generally faster in the 2nd session (for those who got cathodal first, are above the diagonal, meaning anodal &lt; cathodal. The reverse pattern is visible for those who received anodal first). So in that sense there is a difference between two groups of subjects.

#### Baseline differences

If the prosaccade task is indeed reliable, there shouldn't be any significant differences at the group-level between the baseline blocks of both sessions. If there are, this may be problematic, because all the change-scores are expressed relative to the baseline. A significant difference between two conditions in later blocks could then be driven by the differences that already occur in the baseline.

Let's plot the baseline scores in each session for each subject:

``` r
latencyMedianLeg %>%
  filter(leg == "baseline") %>%
  ggplot(aes(stimulation, latency)) +
    facet_grid(type ~ direction) +
    stat_summary(fun.y = mean, geom = "line", aes(group = 1), size = 1.5) +
    geom_line(aes(group = subject, colour = subject), alpha = 0.5, position = position_dodge(width = 0.1)) +
    labs(title = "Baseline in anodal and cathodal sessions", subtitle = "linked stripchart of baseline latencies")
```

<img src="median_latency_files/figure-markdown_github/baseline linked stripchart-1.png" width="75%" style="display: block; margin: auto;" />

Most lines appear fairly flat. However, there is quite some variability around the mean, and some subjects show quite a steep difference. There appears to be some sort of regression to the mean also: extreme scores in either session show a steeper slope than those scores that are less extreme.

Let's look at the pairwise differences in more detail:

``` r
latencyMedianLeg %>%
  filter(leg == "baseline") %>%
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column on session order from other data frame
  group_by(subject,direction,type,session.order) %>%
  summarise(latency.diff = latency[stimulation == "anodal"] - latency[stimulation == "cathodal"]) %>%
  ggplot(aes(factor(0), latency.diff)) +
    facet_grid(type ~ direction) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_summary(fun.data = mean_cl_normal) +
    stat_summary(fun.y = mean, aes(label=round(..y.., digits=2), x = 1.3), geom = "label", alpha = 0.5) +
    geom_point(shape = 21, aes(colour = session.order), position = position_jitter(width=.1)) +
    labs(title = "Baseline in anodal and cathodal sessions", subtitle = "anodal - cathodal")
```

<img src="median_latency_files/figure-markdown_github/baseline difference stripchart-1.png" width="75%" style="display: block; margin: auto;" />

The sequence effect we observed earlier is very clear when we split the data by session order again, particularly in the lateral condition (i.e. latencies are always faster in the 2nd session).

Somewhat worryingly, the baseline differences are in the range of the tDCS effect size we expect, particularly in the center condition. That said, the differences do cluster around 0, there are really only a few subjects that fall far outside the boundary (one is really extreme).

The CIs appear to overlap with 0 though, so the differences shouldn't be significant. Let's test that:

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

Indeed, they are not, although the center conditions come close.

### Line plots per leg from baseline

``` r
kanaiPlotBase <- ggplot(latencyMedianBaseline, aes(leg, latency, color = stimulation, shape = stimulation)) +         
  facet_grid(type ~ direction) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3) +
  coord_cartesian(ylim = c(-10,10))
kanaiPlotBase
```

<img src="median_latency_files/figure-markdown_github/Line plot from baseline-1.png" width="75%" style="display: block; margin: auto;" />

Generally the differences between the stimulation conditions are small (&lt; 4 ms). Except for the center-left condition, but this condition also had the largest baseline difference...

### Individual subjects

#### Anodal session

``` r
kanaiPlotBaseSubsAnodal <- ggplot(latencyMedianBaseline[latencyMedianBaseline$stimulation == "anodal", ], aes(leg, latency)) +
  facet_grid(type ~ direction) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(group=subject,color=subject)) +
  stat_summary(fun.y = mean, aes(group = stimulation), geom = "line") +
  stat_summary(fun.y = mean, geom = "point") + 
  ggtitle("Anodal difference from baseline")
kanaiPlotBaseSubsAnodal
```

<img src="median_latency_files/figure-markdown_github/Line plot from baseline per subject - anodal-1.png" width="75%" style="display: block; margin: auto;" />

#### Cathodal session

``` r
kanaiPlotBaseSubsCathodal <- ggplot(latencyMedianBaseline[latencyMedianBaseline$stimulation == "cathodal", ], aes(leg, latency)) +
  facet_grid(type ~ direction) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(group=subject,color=subject)) +
  stat_summary(fun.y = mean, aes(group = stimulation), geom = "line") +
  stat_summary(fun.y = mean, geom = "point") + 
  ggtitle("Cathodal difference from baseline")
kanaiPlotBaseSubsCathodal
```

<img src="median_latency_files/figure-markdown_github/Line plot from baseline per subject - cathodal-1.png" width="75%" style="display: block; margin: auto;" />

All of this is split pretty much 50-50, hence the average difference hovering around 0.

Statistics
==========

``` r
# Make "subject" a factor, so we can model the repeated measures
latencyMedianBaseline <- latencyMedianBaseline %>%
  ungroup() %>% # remove grouping info, because we need to refactor
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column from on session order from other data frame
  mutate(subject = factor(subject)) # refactor

latencyMedianLeg <- latencyMedianLeg %>%
  ungroup() %>%
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column from on session order from other data frame
  mutate(subject = factor(subject))
```

Frequentist
-----------

### Omnibus anova - saccade latency

**Data**: Outlier trials removed, collapsed into 15-minute intervals.

**Dependent measure**: saccadic latency

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (baseline, tDCS, post.1, post.2)
-   TYPE (lateral vs. center)
-   DIRECTION (left vs. right)

``` r
modelOmni <- ezANOVA(data = data.frame(latencyMedianLeg), # Repeated over subjects; type 3 sums of squares (cf. SPSS)
                        dv = .(latency), wid = .(subject), within = .(stimulation, leg, type, direction), type = 3)
kable(modelOmni$ANOVA)
```

|     | Effect                         |  DFn|  DFd|           F|          p| p&lt;.05 |        ges|
|-----|:-------------------------------|----:|----:|-----------:|----------:|:---------|----------:|
| 2   | stimulation                    |    1|   25|   1.8399816|  0.1870776|          |  0.0046328|
| 3   | leg                            |    3|   75|   0.4395422|  0.7253793|          |  0.0004232|
| 4   | type                           |    1|   25|  28.5410095|  0.0000154| \*       |  0.1295560|
| 5   | direction                      |    1|   25|   0.5725232|  0.4563300|          |  0.0015576|
| 6   | stimulation:leg                |    3|   75|   1.4343849|  0.2394786|          |  0.0007140|
| 7   | stimulation:type               |    1|   25|   0.4929953|  0.4890804|          |  0.0003222|
| 8   | leg:type                       |    3|   75|   4.7880531|  0.0041682| \*       |  0.0015673|
| 9   | stimulation:direction          |    1|   25|   0.6423151|  0.4304257|          |  0.0001620|
| 10  | leg:direction                  |    3|   75|   1.6943654|  0.1754445|          |  0.0003035|
| 11  | type:direction                 |    1|   25|   0.9424122|  0.3409576|          |  0.0008141|
| 12  | stimulation:leg:type           |    3|   75|   0.6925407|  0.5594502|          |  0.0000998|
| 13  | stimulation:leg:direction      |    3|   75|   2.2130816|  0.0935459|          |  0.0003838|
| 14  | stimulation:type:direction     |    1|   25|   0.3842442|  0.5409498|          |  0.0000109|
| 15  | leg:type:direction             |    3|   75|   6.8295412|  0.0003940| \*       |  0.0004416|
| 16  | stimulation:leg:type:direction |    3|   75|   1.0933932|  0.3572648|          |  0.0000761|

``` r
kable(modelOmni$`Mauchly's Test for Sphericity`)
```

|     | Effect                         |          W|          p| p&lt;.05 |
|-----|:-------------------------------|----------:|----------:|:---------|
| 3   | leg                            |  0.2579115|  0.0000058| \*       |
| 6   | stimulation:leg                |  0.3762730|  0.0003186| \*       |
| 8   | leg:type                       |  0.1267739|  0.0000000| \*       |
| 10  | leg:direction                  |  0.3612309|  0.0002082| \*       |
| 12  | stimulation:leg:type           |  0.5303766|  0.0102889| \*       |
| 13  | stimulation:leg:direction      |  0.5729125|  0.0216415| \*       |
| 15  | leg:type:direction             |  0.7933796|  0.3595870|          |
| 16  | stimulation:leg:type:direction |  0.8238579|  0.4676623|          |

``` r
kable(modelOmni$`Sphericity Corrections`)
```

|     | Effect                         |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:-------------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                            |  0.5904749|  0.6231869|                |  0.6322248|  0.6364737|                |
| 6   | stimulation:leg                |  0.5986297|  0.2486500|                |  0.6420287|  0.2482301|                |
| 8   | leg:type                       |  0.4789969|  0.0234582| \*             |  0.5002439|  0.0218426| \*             |
| 10  | leg:direction                  |  0.6068634|  0.1971371|                |  0.6519485|  0.1948609|                |
| 12  | stimulation:leg:type           |  0.6922238|  0.5100514|                |  0.7560487|  0.5218275|                |
| 13  | stimulation:leg:direction      |  0.7236356|  0.1150365|                |  0.7949444|  0.1090952|                |
| 15  | leg:type:direction             |  0.8657955|  0.0008106| \*             |  0.9750655|  0.0004504| \*             |
| 16  | stimulation:leg:type:direction |  0.9008895|  0.3541166|                |  1.0205912|  0.3572648|                |

#### Main effect: type

``` r
latencyMedianLeg %>%
  group_by(subject,type) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(type, latency)) +
  stat_summary(fun.data = mean_cl_normal, size = 1) +
  geom_jitter(width = 0.25)
```

<img src="median_latency_files/figure-markdown_github/Omni Main effect of type-1.png" width="75%" style="display: block; margin: auto;" />

This simply reflects that center saccades are faster than lateral saccades, because the location of the target is known.

#### Interaction: leg by type

``` r
latencyMedianLeg %>%
  group_by(subject,leg,type) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(leg, latency, shape = type)) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line", aes(group = type, linetype = type))
```

<img src="median_latency_files/figure-markdown_github/Omni Interaction leg by type-1.png" width="75%" style="display: block; margin: auto;" />

The effect of Type changes over time: lateral saccades steadily become slower; center saccades vary more randomly.

#### Interaction: leg by type by direction

``` r
latencyMedianLeg %>%
  group_by(subject,leg,type,direction) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(leg, latency, shape = type)) +
  facet_wrap(~direction) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line", aes(group = type, linetype = type))
```

<img src="median_latency_files/figure-markdown_github/Omni Interaction leg by type by direction-1.png" width="75%" style="display: block; margin: auto;" />

The effect over time looks non-linear for left-center saccades. The "tDCS" block is the deviant here, but there is no interaction with stimulation.

Also, for right-lateral saccades, there is a big change that only emerges in the 2nd post-block.

### ANOVA matching Kanai et al. (2012) - lateral saccades

#### Without session order

Differing from the previous omnibus analysis, [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045) analysed shifts from baseline and only had lateral saccades.

**Data**:

-   Outliers removed
-   Collapsed into 15-minute intervals
-   Subtract the baseline from each subsequent block
-   Discard center, keep only lateral saccades

**Dependent measure**: saccadic latency

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

``` r
modelKanai <- ezANOVA(data = data.frame(filter(latencyMedianBaseline, type == "lateral")),
                        dv = .(latency), wid = .(subject), within = .(stimulation,leg,direction), type = 3)

# OR, without the EZ package:
# modelKanai=aov(latency~stimulation*leg*direction + Error(subject/(stimulation*leg*direction)),data=latencyMedianBaselineLateral)
# summary(modelKanai)

kable(modelKanai$ANOVA)
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
kable(modelKanai$`Mauchly's Test for Sphericity`)
```

|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.5213492|  0.0004032| \*       |
| 5   | stimulation:leg           |  0.6360361|  0.0043831| \*       |
| 7   | leg:direction             |  0.7573179|  0.0355909| \*       |
| 8   | stimulation:leg:direction |  0.9123864|  0.3327710|          |

``` r
kable(modelKanai$`Sphericity Corrections`)
```

|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.6762922|  0.0597258|                |  0.7012858|  0.0578252|                |
| 5   | stimulation:leg           |  0.7331572|  0.6087290|                |  0.7674993|  0.6179798|                |
| 7   | leg:direction             |  0.8047111|  0.0434681| \*             |  0.8517313|  0.0406111| \*             |
| 8   | stimulation:leg:direction |  0.9194442|  0.0903598|                |  0.9889658|  0.0859306|                |

##### Main effect of leg

``` r
latencyMedianBaseline %>%
  filter(type == "lateral") %>%
  group_by(subject,leg) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(leg, latency)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.data = mean_cl_normal, size = 1) +
  stat_summary(fun.y = mean, geom = "line", aes(group = 1)) +
  geom_jitter(width = 0.25)
```

<img src="median_latency_files/figure-markdown_github/Kanai Main effect of leg-1.png" width="75%" style="display: block; margin: auto;" />

Saccades become slower over time with respect to the baseline. Note that this effect becomes just non-significant when correcting for sphericity.

##### Interaction: Leg by direction

``` r
latencyMedianBaseline %>%
  filter(type == "lateral") %>%
  group_by(subject,leg,direction) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(leg, latency, shape = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line", aes(group = direction, linetype = direction))
```

<img src="median_latency_files/figure-markdown_github/Kanai Interaction leg by direction-1.png" width="75%" style="display: block; margin: auto;" />

This main effect of leg is stronger for right saccades, but does not occur until the 2nd post-block.

##### Interaction: stimulation by leg by direction

This interaction is not significant, but it is of primary interest, so let's still look at it in more detail.

``` r
latencyMedianBaseline %>%
  filter(type == "lateral") %>%
  group_by(subject,stimulation,leg,direction) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(leg, latency, shape = stimulation, color = stimulation)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~direction) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation))
```

<img src="median_latency_files/figure-markdown_github/Kanai Interaction stimulation by leg by direction-1.png" width="75%" style="display: block; margin: auto;" />

Interestingly, subjects get slower over time in all conditions, but this trend is reversed for left-saccades in the anodal session, as predicted. However, the effect barely reaches 1 ms...

#### With session order

Add an additional factor SESSION ORDER, which creates two groups: those subjects who received anodal tDCS in the first session vs. those who received cathodal tDCS in the first session. Note that these groups are not exactly balanced, which might affect (correcting for) violations of sphericity:

``` r
latencyMedianBaseline %>%
  group_by(session.order) %>%
  summarize(count = n_distinct(subject)) %>%
  kable(.)
```

| session.order  |  count|
|:---------------|------:|
| first.anodal   |     14|
| first.cathodal |     12|

**Data**:

-   Outliers removed
-   Collapsed into 15-minute intervals
-   Subtract the baseline from each subsequent block
-   Discard center, keep only lateral saccades

**Dependent measure**: saccadic latency

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)
-   SESSION ORDER (first anodal vs. first cathodal)

``` r
modelKanaiOrder <- ezANOVA(data = data.frame(filter(latencyMedianBaseline, type == "lateral")), dv = .(latency), 
          wid = .(subject), within = .(stimulation,leg,direction),  between = session.order, type = 3)
kable(modelKanaiOrder$ANOVA)
```

|     | Effect                                  |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:----------------------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | session.order                           |    1|   24|  0.1433284|  0.7083227|          |  0.0020996|
| 3   | stimulation                             |    1|   24|  0.6340723|  0.4336719|          |  0.0068970|
| 5   | leg                                     |    2|   48|  3.1888328|  0.0500837|          |  0.0158928|
| 7   | direction                               |    1|   24|  0.1014256|  0.7528804|          |  0.0004663|
| 4   | session.order:stimulation               |    1|   24|  2.9615300|  0.0981361|          |  0.0314179|
| 6   | session.order:leg                       |    2|   48|  0.9521914|  0.3930621|          |  0.0047991|
| 8   | session.order:direction                 |    1|   24|  0.1559521|  0.6963989|          |  0.0007169|
| 9   | stimulation:leg                         |    2|   48|  0.2675333|  0.7663984|          |  0.0004736|
| 11  | stimulation:direction                   |    1|   24|  0.5243058|  0.4760104|          |  0.0014194|
| 13  | leg:direction                           |    2|   48|  3.4028814|  0.0414916| \*       |  0.0023321|
| 10  | session.order:stimulation:leg           |    2|   48|  3.6820297|  0.0325327| \*       |  0.0064786|
| 12  | session.order:stimulation:direction     |    1|   24|  0.0512017|  0.8229013|          |  0.0001388|
| 14  | session.order:leg:direction             |    2|   48|  0.5664097|  0.5713066|          |  0.0003889|
| 15  | stimulation:leg:direction               |    2|   48|  2.5819394|  0.0860974|          |  0.0030905|
| 16  | session.order:stimulation:leg:direction |    2|   48|  0.7368106|  0.4839727|          |  0.0008839|

``` r
kable(modelKanaiOrder$`Mauchly's Test for Sphericity`)
```

|     | Effect                                  |          W|          p| p&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------|
| 5   | leg                                     |  0.5041662|  0.0003798| \*       |
| 6   | session.order:leg                       |  0.5041662|  0.0003798| \*       |
| 9   | stimulation:leg                         |  0.7028372|  0.0173312| \*       |
| 10  | session.order:stimulation:leg           |  0.7028372|  0.0173312| \*       |
| 13  | leg:direction                           |  0.7680908|  0.0481110| \*       |
| 14  | session.order:leg:direction             |  0.7680908|  0.0481110| \*       |
| 15  | stimulation:leg:direction               |  0.8975070|  0.2883601|          |
| 16  | session.order:stimulation:leg:direction |  0.8975070|  0.2883601|          |

``` r
kable(modelKanaiOrder$`Sphericity Corrections`)
```

|     | Effect                                  |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 5   | leg                                     |  0.6685235|  0.0724403|                |  0.6933380|  0.0704904|                |
| 6   | session.order:leg                       |  0.6685235|  0.3628233|                |  0.6933380|  0.3657030|                |
| 9   | stimulation:leg                         |  0.7709133|  0.7088451|                |  0.8136384|  0.7209771|                |
| 10  | session.order:stimulation:leg           |  0.7709133|  0.0453943| \*             |  0.8136384|  0.0426597| \*             |
| 13  | leg:direction                           |  0.8117481|  0.0526303|                |  0.8622305|  0.0493806| \*             |
| 14  | session.order:leg:direction             |  0.8117481|  0.5370589|                |  0.8622305|  0.5469276|                |
| 15  | stimulation:leg:direction               |  0.9070352|  0.0920383|                |  0.9770103|  0.0875336|                |
| 16  | session.order:stimulation:leg:direction |  0.9070352|  0.4720869|                |  0.9770103|  0.4811399|                |

##### Interaction: session order by stimulation by direction

``` r
latencyMedianBaseline %>%
  filter(type == "lateral") %>%
  group_by(subject,session.order,stimulation,direction) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(direction, latency, shape = stimulation, color = stimulation)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~session.order) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation))
```

<img src="median_latency_files/figure-markdown_github/Kanai Interaction session order by stimulation by direction-1.png" width="75%" style="display: block; margin: auto;" />

The session order by stimulation direction interaction can be construed as a main effect of session. It appears that in the 2nd session, people slow down more than in the 1st. The difference between left and right also seems more pronounced, especially when the 2nd session is anodal.

This latter three-way interaction is very difficult to interpret though. In addition, the stimulation by direction interaction is not significant without session order, so we do not have to worry that this effect is present but confounded by session order.

### ANOVA matching Kanai et al. (2012) - center saccades

#### Without session order

Repeat the same ANOVA, but now for center saccades (which Kanai did not have).

``` r
modelKanaiCenter <- ezANOVA(data = data.frame(filter(latencyMedianBaseline, type == "center")),
                        dv = .(latency), wid = .(subject), within = .(stimulation,leg,direction), type = 3)

kable(modelKanaiCenter$ANOVA)
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
kable(modelKanaiCenter$`Mauchly's Test for Sphericity`)
```

|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.5173954|  0.0003680| \*       |
| 5   | stimulation:leg           |  0.8690648|  0.1856204|          |
| 7   | leg:direction             |  0.9230394|  0.3825099|          |
| 8   | stimulation:leg:direction |  0.7823535|  0.0525819|          |

``` r
kable(modelKanaiCenter$`Sphericity Corrections`)
```

|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.6744887|  0.1241893|                |  0.6991962|  0.1225297|                |
| 5   | stimulation:leg           |  0.8842239|  0.8774896|                |  0.9465498|  0.8902341|                |
| 7   | leg:direction             |  0.9285390|  0.0579604|                |  0.9999608|  0.0537210|                |
| 8   | stimulation:leg:direction |  0.8212564|  0.1600417|                |  0.8713552|  0.1575428|                |

##### Main effect of direction

``` r
latencyMedianBaseline %>%
  filter(type == "center") %>%
  group_by(subject,direction) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(direction, latency)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "line", aes(group = 1), size = 2) +
  geom_line(aes(colour = subject, group = subject))
```

<img src="median_latency_files/figure-markdown_github/Kanai-Center Main effect of direction-1.png" width="75%" style="display: block; margin: auto;" />

This seems to be a really tiny and variable effect, but apparently left saccades get somewhat slower and right saccades somewhat faster compared to baseline.

#### With session order

``` r
modelKanaiCenterOrder <- ezANOVA(data = data.frame(filter(latencyMedianBaseline, type == "center")), dv = .(latency), 
          wid = .(subject), within = .(stimulation,leg,direction),  between = session.order, type = 3)
kable(modelKanaiCenterOrder$ANOVA)
```

|     | Effect                                  |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:----------------------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | session.order                           |    1|   24|  0.1166399|  0.7356800|          |  0.0025860|
| 3   | stimulation                             |    1|   24|  3.3377548|  0.0801702|          |  0.0258116|
| 5   | leg                                     |    2|   48|  2.3943806|  0.1020456|          |  0.0066197|
| 7   | direction                               |    1|   24|  6.5039074|  0.0175653| \*       |  0.0189442|
| 4   | session.order:stimulation               |    1|   24|  0.9506165|  0.3392880|          |  0.0074896|
| 6   | session.order:leg                       |    2|   48|  0.5737494|  0.5672253|          |  0.0015943|
| 8   | session.order:direction                 |    1|   24|  5.5454061|  0.0270448| \*       |  0.0161976|
| 9   | stimulation:leg                         |    2|   48|  0.0965518|  0.9081388|          |  0.0001683|
| 11  | stimulation:direction                   |    1|   24|  2.4540824|  0.1303112|          |  0.0068754|
| 13  | leg:direction                           |    2|   48|  2.8467499|  0.0678686|          |  0.0016019|
| 10  | session.order:stimulation:leg           |    2|   48|  0.5117829|  0.6026597|          |  0.0008915|
| 12  | session.order:stimulation:direction     |    1|   24|  3.2609785|  0.0835059|          |  0.0091154|
| 14  | session.order:leg:direction             |    2|   48|  1.1056981|  0.3392615|          |  0.0006228|
| 15  | stimulation:leg:direction               |    2|   48|  2.2553582|  0.1158348|          |  0.0013961|
| 16  | session.order:stimulation:leg:direction |    2|   48|  1.7398007|  0.1864427|          |  0.0010773|

``` r
kable(modelKanaiCenterOrder$`Mauchly's Test for Sphericity`)
```

|     | Effect                                  |          W|          p| p&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------|
| 5   | leg                                     |  0.5132456|  0.0004664| \*       |
| 6   | session.order:leg                       |  0.5132456|  0.0004664| \*       |
| 9   | stimulation:leg                         |  0.8520126|  0.1585365|          |
| 10  | session.order:stimulation:leg           |  0.8520126|  0.1585365|          |
| 13  | leg:direction                           |  0.8959726|  0.2827415|          |
| 14  | session.order:leg:direction             |  0.8959726|  0.2827415|          |
| 15  | stimulation:leg:direction               |  0.8087403|  0.0870572|          |
| 16  | session.order:stimulation:leg:direction |  0.8087403|  0.0870572|          |

``` r
kable(modelKanaiCenterOrder$`Sphericity Corrections`)
```

|     | Effect                                  |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 5   | leg                                     |  0.6726060|  0.1234676|                |  0.6980931|  0.1217371|                |
| 6   | session.order:leg                       |  0.6726060|  0.5031425|                |  0.6980931|  0.5090948|                |
| 9   | stimulation:leg                         |  0.8710897|  0.8835886|                |  0.9334806|  0.8962764|                |
| 10  | session.order:stimulation:leg           |  0.8710897|  0.5781127|                |  0.9334806|  0.5903977|                |
| 13  | leg:direction                           |  0.9057746|  0.0738118|                |  0.9754789|  0.0693693|                |
| 14  | session.order:leg:direction             |  0.9057746|  0.3351508|                |  0.9754789|  0.3382625|                |
| 15  | stimulation:leg:direction               |  0.8394475|  0.1256368|                |  0.8953942|  0.1221703|                |
| 16  | session.order:stimulation:leg:direction |  0.8394475|  0.1923135|                |  0.8953942|  0.1903485|                |

###### Interaction: session order by direction

There is a main effect of direction in the ANOVA without session order, so this effect might be confounded.

``` r
latencyMedianBaseline %>%
  filter(type == "center") %>%
  group_by(subject,session.order,direction) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(direction, latency, shape = session.order)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line", aes(group = session.order, linetype = session.order))
```

<img src="median_latency_files/figure-markdown_github/Kanai-Center-Session-order Interaction with direction-1.png" width="75%" style="display: block; margin: auto;" />

The main effect of direction is present for both groups, but it is much stronger in the group that first receives cathodal stimulation. This could be construed as a *carryover* effect, i.e. that the effect of cathodal stimulation was not yet washed-out by the second session and modulated the effect of anodal (and vice versa). Or, it could simply be that there is a difference between the two groups (anodal first vs. cathodal first) that has nothing to do with stimulation.

Bayesian
--------

### Linear mixed effects matching Kanai - saccade latency

#### Test against the null model

``` r
bfKanaiLateral = anovaBF(latency~stimulation*leg*direction+subject, data = data.frame(filter(latencyMedianBaseline, type == "lateral")), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiLateral = sort(bfKanaiLateral, decreasing = TRUE) # sort such that winning model is at the top
```

First we compare all models to the most simple (null) model, which is the intercept only + random effect model: `latency ~ subject`. This does not test for effects of SUBJECT but models it as a nuisance factor (`whichRandom = "subject"`). In addition, to decrease the model space, we do not consider models that have an interaction without the corresponding main effects (`whichModels = "withmain"`).

``` r
kable(select(extractBF(bfKanaiLateral), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |         bf|
|-------------------------------------------------------------------------------------------------------------------------------|----------:|
| leg + subject                                                                                                                 |  0.8369410|
| stimulation + subject                                                                                                         |  0.7471023|
| stimulation + leg + subject                                                                                                   |  0.6548069|
| direction + subject                                                                                                           |  0.1329814|
| direction + leg + subject                                                                                                     |  0.1156958|
| stimulation + direction + subject                                                                                             |  0.0994426|
| stimulation + direction + leg + subject                                                                                       |  0.0871115|
| stimulation + leg + stimulation:leg + subject                                                                                 |  0.0470574|
| stimulation + direction + stimulation:direction + subject                                                                     |  0.0221274|
| stimulation + direction + stimulation:direction + leg + subject                                                               |  0.0198278|
| direction + leg + direction:leg + subject                                                                                     |  0.0108875|
| stimulation + direction + leg + direction:leg + subject                                                                       |  0.0087366|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |  0.0063617|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |  0.0018489|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |  0.0013926|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |  0.0006423|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |  0.0001351|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |  0.0000248|

The winning model is the one with just a main effect of LEG, with a Bayes factor of 0.8. However, so even for this model there is less evidence than for the null model, which by definition has a Bayes Factor of 1.

The conventional interpretation for Bayes Factors is the following:

-   BF<sub>10</sub> &lt; 0.1: strong evidence for H0
-   0.1 &lt; BF<sub>10</sub> &lt; .33: moderate evidence for H0
-   0.333 &lt; BF<sub>10</sub> &lt; 1: anecdotal evidence for H0

-   BF<sub>10</sub> = 1: equivalent evidence for H0 and H1

-   1 BF<sub>10</sub> &lt; 3: anecdotal evidence for H1
-   3 &lt; BF<sub>10</sub> &lt; 10: moderate evidence for H1
-   BF<sub>10</sub> &gt; 10: strong evidence H1

So to switch between expressing the Bayes Factor as evidence for H1 (BF<sub>10</sub>) or H0 (BF<sub>01</sub>), you simply invert it (divide by 1).

We can compute the evidence for a particular effect by comparing this winning model with the best-fitting model that does *not* contain the effect. We can compute the evidence for *absence* of a particular effect by comparing the winning model with the best-fitting model that *does* contain the effect.

The evidence for the effect of LEG can be quantified by comparing the Bayes factors of the first and the 3rd model (because that is the first one that does *not* contain LEG as a factor): 1.3. This constitues only anocdotal evidence for the presence of a leg effect, even though it was significant in the classical ANOVA. But, given that the leg model itself is such a poor fit, in the Bayesian analysis this weak evidence is no surprise.

The evidence for the *absence* of the STIMULATION by DIRECTION interaction can be quantified by comparing the Bayes factors of the first and the 9th model (because that is the first one that *does* contain the interaction): 37.8. This constitues strong evidence for the absence of the interaction.

#### Test against the full model

Another option for quantifying evidence for a particular effect is to compare the full model to a model where that effect is omitted (`whichModels = top")`. The full model is `stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject`.

``` r
bfKanaiLateralFull = anovaBF(latency~stimulation*leg*direction+subject, data = data.frame(filter(latencyMedianBaseline, type == "lateral")), whichModels = "top", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiLateralFull
```

    ## Bayes factor top-down analysis
    ## --------------
    ## When effect is omitted from stimulation + direction + leg + stimulation:direction + stimulation:leg + direction:leg + stimulation:direction:leg + subject , BF is...
    ## [1] Omit direction:leg:stimulation : 4.955423  <U+00B1>11.9%
    ## [2] Omit direction:leg             : 8.543017  <U+00B1>11.68%
    ## [3] Omit leg:stimulation           : 11.8334   <U+00B1>11.89%
    ## [4] Omit direction:stimulation     : 3.805628  <U+00B1>11.64%
    ## [5] Omit leg                       : 0.9760963 <U+00B1>12.17%
    ## [6] Omit direction                 : 6.304331  <U+00B1>11.75%
    ## [7] Omit stimulation               : 1.285336  <U+00B1>21.89%
    ## 
    ## Against denominator:
    ##   latency ~ stimulation + direction + leg + stimulation:direction + stimulation:leg + direction:leg + stimulation:direction:leg +     subject 
    ## ---
    ## Bayes factor type: BFlinearModel, JZS

Removing the LEG effect from the model yields a higher Bayes Factor, so removing this effect actually improved the model, although only by a little bit. The evidence thus goes in the direction of the null; when expressed in favor of the alternative, the Bayes Factir becomes less than one: 1  0.976 `=` 1

Similarly, removing the DIRECTION by STIMULATION effect improves the model, and this time a bit more: in the range of moderate evidence for the null.

#### Bayesian model averaging

The problem with the 1st option (comparing single models, for example to the winning model) is that for designs with many factors (and therefore models), it becomes risky and a bit subjective to base conclusions on just two models. In this case, the winning model is even a bad fit, so it doesn't seem appropriate to use it as a benchmark.

The problem with the 2nd option (comparing against the full model) is actually very apparent in this dataset: the full model is a terrible fit, as it comes in last. There is even a lot of evidence against it when compared to the null model!

One solution is to do Bayesian Model Averaging: compare multiple models and aggregate the Bayes Factors.

In the [JASP stats package](http://jasp-stats.org), this analysis is also called the "inclusion Bayes factor". Briefly, it compares all models that include the effect of interest vs. all models that do not. For examples and a conceptual explanation, see example 5 in: [Bayesian inference for psychology. Part II: Example applications with JASP](https://doi.org/10.3758/s13423-017-1323-7).

There's also a different way of calculating inclusion Bayes factors ([conceptualized by Sebastiaan Mathot](https://www.cogsci.nl/blog/interpreting-bayesian-repeated-measures-in-jasp)) that is currently being implemented in JASP. This is called the "inclusion Bayes factor across matched models", and is more selective in the set of models that's being compared than the standard inclusion BF. Briefly, the inclusion BF across matched models compares:

-   all models that include the effect of interest, but NO interactions with the effect of interest, VERSUS
-   the models that result from stripping the effect of interest from this set of models.

Let's compare the default inclusion Bayes Factors for all effects:

``` r
# Default inclusion Bayes Factors
kable(inclusionBF(bfKanaiLateral))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |     0.2890214|
| direction                 |     0.0551289|
| stimulation:direction     |     0.0262216|
| leg                       |     0.3196411|
| stimulation:leg           |     0.0322399|
| direction:leg             |     0.0127990|
| stimulation:direction:leg |     0.0001179|

Doing the analysis this way, we also find strong evidence against most of these effects; particularly the interactions.

My preferred approach is the inclusion Bayes factor across matched models. The evidence is qualitatively similar, but less strong, probably because we aren't including as many poorly fitting models in the model comparison (which have a very large BF<sub>10</sub> and therefore a lot of effect on the composite Bayes Factor):

``` r
kable(inclusionBF(bfKanaiLateral, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |     0.7618390|
| direction                 |     0.1343900|
| stimulation:direction     |     0.2240872|
| leg                       |     0.8564833|
| stimulation:leg           |     0.0719756|
| direction:leg             |     0.0965773|
| stimulation:direction:leg |     0.1838927|

### Linear mixed effects matching Kanai - center saccades

``` r
bfKanaiCenter = anovaBF(latency~stimulation*leg*direction+subject, data = data.frame(filter(latencyMedianBaseline, type == "center")), whichModels="withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiCenter <- sort(bfKanaiCenter, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfKanaiCenter), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |           bf|
|-------------------------------------------------------------------------------------------------------------------------------|------------:|
| stimulation + direction + subject                                                                                             |  536.0919637|
| stimulation + direction + stimulation:direction + subject                                                                     |  401.9954582|
| stimulation + direction + leg + subject                                                                                       |   96.5673070|
| stimulation + direction + stimulation:direction + leg + subject                                                               |   71.4297259|
| stimulation + subject                                                                                                         |   56.3933677|
| stimulation + leg + subject                                                                                                   |    9.4872949|
| stimulation + direction + leg + direction:leg + subject                                                                       |    9.0981205|
| direction + subject                                                                                                           |    8.0433468|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |    6.7501094|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |    6.4583019|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |    4.7668502|
| direction + leg + subject                                                                                                     |    1.3128487|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |    0.7327462|
| stimulation + leg + stimulation:leg + subject                                                                                 |    0.6264923|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |    0.4730413|
| leg + subject                                                                                                                 |    0.1575411|
| direction + leg + direction:leg + subject                                                                                     |    0.1387820|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |    0.0654546|

Interestingly, there is a lot more evidence across the board, especially for the models that feature Stimulation and Direction.

``` r
kable(inclusionBF(bfKanaiCenter, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |    66.4291778|
| direction                 |     9.5836353|
| stimulation:direction     |     0.7480027|
| leg                       |     0.1783263|
| stimulation:leg           |     0.0675387|
| direction:leg             |     0.0952325|
| stimulation:direction:leg |     0.1383697|

There seems to be a mismatch with the classical ANOVA: An effect of stimulation receives strong support, whereas it was non-significant. The effect of direction was significant, but it is less strongly supported by the inclusion Bayes Factor (although the evidence still goes in the right direction).

#### Main effect of stimulation

``` r
latencyMedianBaseline %>%
  filter(type == "center") %>%
  group_by(subject,stimulation) %>%
  summarise(latency = mean(latency)) %>%
  ggplot(aes(stimulation, latency)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "line", aes(group = 1), size = 2) +
  geom_line(aes(colour = subject, group = subject))
```

<img src="median_latency_files/figure-markdown_github/BF Kanai-Center Main effect of stimulation-1.png" width="75%" style="display: block; margin: auto;" />

It is easy to see why this effect is non-significant: the average difference is tiny and there is a lot of variability. The plot shows one major outlier though in terms of the effect size (S01): let's see what happens to the Bayes Factor if we take their data out.

``` r
latencyNoS01 <- latencyMedianBaseline %>%
  filter(subject != "S01") %>%
  mutate(subject = factor(subject))

bfKanaiCenterNoS01 = anovaBF(latency~stimulation*leg*direction+subject, data = data.frame(filter(latencyNoS01, type == "center")), whichModels="withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiCenterNoS01 <- sort(bfKanaiCenterNoS01, decreasing = TRUE) # sort such that winning model is at the top
kable(select(extractBF(bfKanaiCenterNoS01), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |            bf|
|-------------------------------------------------------------------------------------------------------------------------------|-------------:|
| stimulation + direction + subject                                                                                             |  1414.9173286|
| direction + subject                                                                                                           |   593.9260628|
| stimulation + direction + leg + subject                                                                                       |   500.8782055|
| stimulation + direction + stimulation:direction + subject                                                                     |   355.6182797|
| direction + leg + subject                                                                                                     |   207.5273825|
| stimulation + direction + stimulation:direction + leg + subject                                                               |   129.1617435|
| stimulation + direction + leg + direction:leg + subject                                                                       |    66.2674131|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |    34.6189068|
| direction + leg + direction:leg + subject                                                                                     |    24.8226424|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |    15.0123673|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |     8.7836953|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |     4.5161916|
| stimulation + subject                                                                                                         |     2.0143682|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |     1.0355578|
| stimulation + leg + subject                                                                                                   |     0.6306025|
| leg + subject                                                                                                                 |     0.2998792|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |     0.1608650|
| stimulation + leg + stimulation:leg + subject                                                                                 |     0.0414637|

``` r
kable(inclusionBF(bfKanaiCenterNoS01, models = "matched"))
```

| effect                    |  Bayes.factor|
|:--------------------------|-------------:|
| stimulation               |     2.3982184|
| direction                 |   690.3290091|
| stimulation:direction     |     0.2521335|
| leg                       |     0.3541737|
| stimulation:leg           |     0.0688191|
| direction:leg             |     0.1267400|
| stimulation:direction:leg |     0.1553414|

This completely abolishes the strong support for Stimulation, and greatly enhances the support for Direction, bringing the Bayesian and the classical ANOVAs more in line. At present it is unclear why the classical and Bayesian analyses differ in this regard. When simulating this case for normally distributed data, the Bayes Factors and p-values track each other nicely (see [discussion on JASP/BayesFactor forum](http://forum.cogsci.nl/index.php?p=/discussion/3596/large-bayes-factor-changes-with-exclusion-of-single-subject-bayesian-anova)). This suggests there must be some assumption that is not met in this particular dataset, which is causing the divergence between the analyses.

Still, let's do some follow-up tests with all of the data (i.e. including this subject) to see whether the anodal or cathodal change scores are significantly different from 0 on their own.

Bayesian one-sample t-tests:

``` r
latencyMedianBaseline %>%
  filter(type == "center") %>% # keep only center saccades
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(deviation.end = mean(latency)) %>% # average over all other variables
  spread(stimulation,deviation.end) %>% # make separate columns with test data
  summarise_if(is.numeric, funs(extractBF(ttestBF(.), onlybf = TRUE))) %>% # run Bayesian t-test on each column, keeping only the BF
  gather(stimulation,BF,anodal,cathodal) %>% # make row for each stimulation condition
  kable(.)
```

| stimulation |         BF|
|:------------|----------:|
| anodal      |  0.3218926|
| cathodal    |  0.2840670|

Frequentist one-sample t-tests:

``` r
latencyMedianBaseline %>%
  filter(type == "center") %>% # keep only center saccades
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(deviation.end = mean(latency)) %>% # average over all other variables (df is now still grouped per stimulation)
  summarise_if(is.numeric, funs(list(tidy(t.test(.))))) %>%  # run one-sample t-test for each stimulation condition, return tidy data frames
  unnest() %>% # unpack the list-column with data frame for each test
  kable(.)
```

| stimulation |   estimate|   statistic|    p.value|  parameter|   conf.low|  conf.high| method            | alternative |
|:------------|----------:|-----------:|----------:|----------:|----------:|----------:|:------------------|:------------|
| anodal      |  -1.608974|  -0.9875463|  0.3328370|         25|  -4.964508|   1.746559| One Sample t-test | two.sided   |
| cathodal    |   1.096154|   0.8335773|  0.4124131|         25|  -1.612139|   3.804446| One Sample t-test | two.sided   |

So neither are actually significant or have a BF with evidence for the alternative.
