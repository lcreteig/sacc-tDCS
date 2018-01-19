sacc-tDCS: Saccade accuracy
================
Leon Reteig

-   [Load data](#load-data)
    -   [Load eye data](#load-eye-data)
    -   [Subject metadata](#subject-metadata)
-   [Preprocess data](#preprocess-data)
    -   [Outliers](#outliers)
    -   [Cut into 15-minute sections](#cut-into-15-minute-sections)
-   [Saccade end point deviation](#saccade-end-point-deviation)
    -   [Prepare data frame for plotting & statistics](#prepare-data-frame-for-plotting-statistics)
    -   [Plot](#plot)
        -   [With baseline block](#with-baseline-block)
        -   [Baseline subtracted](#baseline-subtracted)
    -   [Statistics](#statistics)
        -   [Frequentist](#frequentist)
        -   [Bayesian](#bayesian)
-   [Saccade end point variability](#saccade-end-point-variability)
    -   [Calculate endpoint variability](#calculate-endpoint-variability)
    -   [Prepare data frame for plotting & statistics](#prepare-data-frame-for-plotting-statistics-1)
    -   [Plot](#plot-1)
        -   [With baseline block](#with-baseline-block-1)
        -   [Baseline subtracted](#baseline-subtracted-1)
    -   [Statistics](#statistics-1)
        -   [Frequentist](#frequentist-1)
        -   [Bayesian](#bayesian-1)

R notebook for inspection of data and analyses of saccade end point deviation and variability in the `sacc-tDCS` dataset. Previous processing:

-   Raw data were parsed into events (saccades, fixations, etc.) by the EyeLink data were collected on.
-   Events were extracted and saccade measures were computed with a MATLAB script.

``` r
# Load some libraries
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
library(forcats) # manipulatin factors
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
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 6, fig.asp = 0.618, out.width = "70%", fig.align = "center")

source("src/lib/InclusionBF.R")

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
    ##  [7] forcats_0.2.0        dplyr_0.5.0          purrr_0.2.2         
    ## [10] readr_1.1.0          tidyr_0.6.1          tibble_1.3.0        
    ## [13] ggplot2_2.2.1        tidyverse_1.1.1     
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
    ## [25] evaluate_0.10      SparseM_1.77       quantreg_5.33     
    ## [28] pbkrtest_0.4-7     parallel_3.4.0     Rcpp_0.12.10      
    ## [31] scales_0.4.1       backports_1.0.5    jsonlite_1.4      
    ## [34] lme4_1.1-13        mnormt_1.5-5       hms_0.3           
    ## [37] digest_0.6.12      stringi_1.1.5      grid_3.4.0        
    ## [40] rprojroot_1.2      tools_3.4.0        magrittr_1.5      
    ## [43] lazyeval_0.2.0     car_2.1-4          MASS_7.3-47       
    ## [46] xml2_1.1.1         lubridate_1.6.0    assertthat_0.2.0  
    ## [49] minqa_1.2.4        rmarkdown_1.5      httr_1.2.1        
    ## [52] R6_2.2.0           nnet_7.3-12        nlme_3.1-131      
    ## [55] compiler_3.4.0

Load data
=========

Load eye data
-------------

The .csv file with the eye tracking data was created in MATLAB.

``` r
# Load the data frame
dataFile <- file.path("data", "sacc-tDCS_data.csv")
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

Subject metadata
----------------

``` r
# Load eye tracking data into data frame
dataFile <- file.path("data", "subject_info.csv")
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

The main use is to see if the nuisance factor *session.order* covaries with the factors of interest in the design. This could indicate the presence of carryover effects between the stimulation, or a difference in subgroups within the sample (see <http://www.jerrydallal.com/lhsp/crossovr.htm> for an introduction to these kinds of analyses.).

Preprocess data
===============

Outliers
--------

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

Cut into 15-minute sections
---------------------------

Cut the post-block into two so we have four 15-minute intervals: one before, one during, and two after stimulation.

``` r
# Split the "post" leg into two
groupData <- mutate(groupData,
                    leg = as.character(leg), # cannot edit leg if it's still a factor
                    leg = replace(leg, leg == "post" & block <= 3, "post.1"),
                    leg = replace(leg, block > 3, "post.2"),
                    leg = factor(leg, levels = c("pre", "tDCS", "post.1", "post.2")) # refactor and order levels
                    )
```

Saccade end point deviation
===========================

One estimate of the accuracy of saccades is the mean landing position with respect to the target location. [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045) also examined this, but found no effects of tDCS.

The simplest measure (which Kanai et al. (2012) also used) is the Euclidian distance (shortest straight line) between the saccade end point and the center of the target stimulus. We already have the deviations in the x- and y- directions in degrees of visual angle. Now we just need to calculate the length of the vector.

``` r
# Calculate end point deviation
devData <- mutate(groupData, deviation.end = sqrt(deviation.end.x^2 + deviation.end.y^2))
```

Prepare data frame for plotting & statistics
--------------------------------------------

Average over three blocks:

``` r
devData <- devData %>%
  group_by(subject,stimulation,direction,type) %>% 
  summarise(baseline = mean(deviation.end[leg == "pre"]), # take average of 3 blocks, make new column
            tDCS = mean(deviation.end[leg == "tDCS"]),
            post.1 = mean(deviation.end[leg == "post.1"]),
            post.2 = mean(deviation.end[leg == "post.2"])) %>%
gather(leg, deviation.end, baseline, tDCS, post.1, post.2)  %>% # gather new columns to use as factor 
mutate(leg = factor(leg, levels = c("baseline", "tDCS", "post.1", "post.2"))) # reorder factor levels
```

Subtract the baseline from each average:

``` r
# Subtract baseline
devDataBase <- devData %>%
  group_by(subject,stimulation,direction,type) %>% 
  # subtract baseline block from others, make new column
  summarise(tDCS = deviation.end[leg == "tDCS"] - deviation.end[leg == "baseline"], 
            post.1 = deviation.end[leg == "post.1"] - deviation.end[leg == "baseline"],
            post.2 = deviation.end[leg == "post.2"] - deviation.end[leg == "baseline"]) %>%
gather(leg, deviation.end, tDCS, post.1, post.2) %>% # gather new columns to use as factor 
mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) # reorder factor levels
```

Plot
----

### With baseline block

``` r
kanaiPlotDev <- ggplot(devData, aes(leg, deviation.end, color = stimulation, shape = stimulation)) +         
  facet_grid(type ~ direction) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3)
kanaiPlotDev
```

<img src="accuracy_files/figure-markdown_github/Line plot per leg - deviation-1.png" width="70%" style="display: block; margin: auto;" />

At first glance there don't seem to be many differences that are larger than the baseline differences and/or relate clearly to the polarity or timing of stimulation.

Let's look at the individual subject data:

``` r
kanaiPlotSubsAnodal <- ggplot(devData[devData$stimulation == "anodal", ], aes(leg, deviation.end)) +
  facet_grid(type ~ direction) +
  geom_line(aes(group = subject,color = subject)) +
  stat_summary(fun.y = mean, aes(group = stimulation), geom = "line") +
  stat_summary(fun.y = mean, geom = "point") +
  ggtitle("Anodal session")
kanaiPlotSubsAnodal
```

<img src="accuracy_files/figure-markdown_github/Line plot per subject - deviation anodal-1.png" width="70%" style="display: block; margin: auto;" />

``` r
kanaiPlotSubsCathodal <- ggplot(devData[devData$stimulation == "cathodal", ], aes(leg, deviation.end)) +
  facet_grid(type ~ direction) +
  geom_line(aes(group = subject,color = subject)) +
  stat_summary(fun.y = mean, aes(group = stimulation), geom = "line") +
  stat_summary(fun.y = mean, geom = "point") +
  ggtitle("Cathodal session")
kanaiPlotSubsCathodal
```

<img src="accuracy_files/figure-markdown_github/Line plot per subject - deviation cathodal-1.png" width="70%" style="display: block; margin: auto;" />

There are definitely some outliers, but mostly in terms of overall offset / baseline differences.

#### Baseline reliability

Scatterplot and correlation of baseline data in the two sessions:

``` r
baselineCorrDev <- devData %>%
  filter(leg == "baseline") %>% 
  group_by(direction,type) %>% 
  spread(stimulation,deviation.end) %>% 
  nest() %>% 
  mutate(stats = map(data, ~cor.test(formula = ~ anodal + cathodal, data =.))) %>% # run correlation test on baselines from each condition
  mutate(tidy_model = map(stats, tidy)) %>% 
  unnest(tidy_model, .drop = TRUE)
```

``` r
devData %>%
  filter(leg == "baseline") %>%
  spread(stimulation,deviation.end) %>%
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column on session order from other data frame
  ggplot(aes(anodal,cathodal)) +
    facet_grid(type ~ direction) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_smooth(method = "lm") +
    geom_point(aes(color=session.order)) +
    xlim(0,2) + ylim(0,2) +
    geom_text(data = baselineCorrDev, x = 0.2, y = 1.5, aes(label = paste("italic(r) == ", round(estimate,2))), parse = TRUE) +
    labs(title = "Baseline in anodal and cathodal sessions", subtitle = "scatterplot of baseline endpoint deviation")
```

<img src="accuracy_files/figure-markdown_github/baseline scatterplots - deviation-1.png" width="70%" style="display: block; margin: auto;" />

The correlations are not as high as for the latency data, but still reasonable. The sequence effect we observed in the latency data is not so prominent, so apparently there's less of a practice effect in saccade enpdoint deviation.

#### Baseline differences

``` r
devData %>%
  filter(leg == "baseline") %>%
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column on session order from other data frame
  group_by(subject,direction,type,session.order) %>%
  summarise(deviation.end.diff = deviation.end[stimulation == "anodal"] - deviation.end[stimulation == "cathodal"]) %>%
  ggplot(aes(factor(0), deviation.end.diff)) +
    facet_grid(type ~ direction) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_summary(fun.data = mean_cl_normal) +
    stat_summary(fun.y = mean, aes(label=round(..y.., digits=2), x = 1.3), geom = "label", alpha = 0.5) +
    geom_point(shape = 21, aes(colour = session.order), position = position_jitter(width=.1)) +
    labs(title = "Baseline in anodal and cathodal sessions", subtitle = "anodal - cathodal")
```

<img src="accuracy_files/figure-markdown_github/baseline difference stripchart - deviation-1.png" width="70%" style="display: block; margin: auto;" />

The baseline differences are not so extreme, except for the center (left) condition: that seems quite large and consistent over subjects.

``` r
devData %>%
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

Indeed, in the center-left condition the baseline difference is significant, and it's at trend in the center-right condition.

### Baseline subtracted

``` r
kanaiPlotDevBase <- ggplot(devDataBase, aes(leg, deviation.end, color = stimulation, shape = stimulation)) +      
  facet_grid(type ~ direction) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3)
kanaiPlotDevBase
```

<img src="accuracy_files/figure-markdown_github/Line plot from baseline - deviation-1.png" width="70%" style="display: block; margin: auto;" />

This clearly shows that all the changes are quite tiny (less than 0.15 degrees of visual angle). There appears to be a clear difference in between the anodal and cathodal change scores for center saccades (and maybe for left-lateral saccades).

However, we know that the difference between anodal and cathodal is actually maximal in the baseline. Thus it remains unclear whether the baseline differences are spurious and the effect is real, or whether the "effect" is driven by the baseline difference (i.e. something akin to regression to the mean).

Statistics
----------

``` r
# Make "subject" a factor, so we can model the repeated measures
devDataBase <- devDataBase %>%
  ungroup() %>% # remove any grouping info, because we need to refactor
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column on session order from other data frame
  mutate(subject = factor(subject)) # refactor
```

### Frequentist

#### ANOVA matching Kanai et al. (2012) - lateral saccades

##### Without session order

**Data**:

-   Outliers removed
-   Collapsed into 15-minute intervals
-   Subtract the baseline from each subsequent block
-   Discard center, keep only lateral saccades

**Dependent measure**: saccade end point deviation

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

``` r
modelKanai <- ezANOVA(data = data.frame(filter(devDataBase, type == "lateral")),
                        dv = .(deviation.end), wid = .(subject), within = .(stimulation,leg,direction), type = 3)

kable(modelKanai)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect                    |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  2.0266072|  0.1669269|          |  0.0178886|
| 3   | leg                       |    2|   50|  0.8985792|  0.4136202|          |  0.0021712|
| 4   | direction                 |    1|   25|  1.5739544|  0.2212365|          |  0.0059023|
| 5   | stimulation:leg           |    2|   50|  0.5863132|  0.5601540|          |  0.0015026|
| 6   | stimulation:direction     |    1|   25|  0.1283222|  0.7231851|          |  0.0005433|
| 7   | leg:direction             |    2|   50|  0.1400779|  0.8696305|          |  0.0002306|
| 8   | stimulation:leg:direction |    2|   50|  0.2773853|  0.7589209|          |  0.0002720|

</td>
<td>
|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.6913883|  0.0119307| \*       |
| 5   | stimulation:leg           |  0.9959286|  0.9522226|          |
| 7   | leg:direction             |  0.9738561|  0.7276748|          |
| 8   | stimulation:leg:direction |  0.7556773|  0.0346766| \*       |

</td>
<td>
|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.7641686|  0.3907386|                |  0.8038793|  0.3951392|                |
| 5   | stimulation:leg           |  0.9959451|  0.5595027|                |  1.0819912|  0.5601540|                |
| 7   | leg:direction             |  0.9745222|  0.8646202|                |  1.0558164|  0.8696305|                |
| 8   | stimulation:leg:direction |  0.8036501|  0.7108179|                |  0.8504748|  0.7234385|                |

</td>
</tr>
</tbody>
</table>
##### With session order

Add an additional factor SESSION ORDER, which creates two groups: those subjects who received anodal tDCS in the first session vs. those who received cathodal tDCS in the first session. Note that these groups are not exactly balanced, which might affect (correcting for) violations of sphericity:

``` r
devDataBase %>%
  group_by(session.order) %>%
  summarize(count = n_distinct(subject)) %>%
  kable(.)
```

| session.order  |  count|
|:---------------|------:|
| first.anodal   |     14|
| first.cathodal |     12|

``` r
modelKanaiOrder <- ezANOVA(data = data.frame(filter(devDataBase, type == "lateral")), dv = .(deviation.end), 
          wid = .(subject), within = .(stimulation,leg,direction),  between = session.order, type = 3)
kable(modelKanaiOrder)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect                                  |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:----------------------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | session.order                           |    1|   24|  0.2049021|  0.6548582|          |  0.0033510|
| 3   | stimulation                             |    1|   24|  2.2664083|  0.1452548|          |  0.0204847|
| 5   | leg                                     |    2|   48|  0.8688109|  0.4259409|          |  0.0022474|
| 7   | direction                               |    1|   24|  1.5405279|  0.2265344|          |  0.0062015|
| 4   | session.order:stimulation               |    1|   24|  1.1398125|  0.2963138|          |  0.0104081|
| 6   | session.order:leg                       |    2|   48|  0.1072400|  0.8985246|          |  0.0002780|
| 8   | session.order:direction                 |    1|   24|  0.0363616|  0.8503745|          |  0.0001473|
| 9   | stimulation:leg                         |    2|   48|  0.6189985|  0.5427260|          |  0.0016790|
| 11  | stimulation:direction                   |    1|   24|  0.2640686|  0.6120378|          |  0.0010552|
| 13  | leg:direction                           |    2|   48|  0.0746248|  0.9281991|          |  0.0001222|
| 10  | session.order:stimulation:leg           |    2|   48|  0.3828206|  0.6839985|          |  0.0010391|
| 12  | session.order:stimulation:direction     |    1|   24|  3.3325860|  0.0803899|          |  0.0131556|
| 14  | session.order:leg:direction             |    2|   48|  1.9582493|  0.1522160|          |  0.0031959|
| 15  | stimulation:leg:direction               |    2|   48|  0.2604897|  0.7717565|          |  0.0002689|
| 16  | session.order:stimulation:leg:direction |    2|   48|  0.5209071|  0.5973007|          |  0.0005376|

</td>
<td>
|     | Effect                                  |          W|          p| p&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------|
| 5   | leg                                     |  0.6837007|  0.0126171| \*       |
| 6   | session.order:leg                       |  0.6837007|  0.0126171| \*       |
| 9   | stimulation:leg                         |  0.9937065|  0.9299690|          |
| 10  | session.order:stimulation:leg           |  0.9937065|  0.9299690|          |
| 13  | leg:direction                           |  0.9622853|  0.6426792|          |
| 14  | session.order:leg:direction             |  0.9622853|  0.6426792|          |
| 15  | stimulation:leg:direction               |  0.7507047|  0.0369739| \*       |
| 16  | session.order:stimulation:leg:direction |  0.7507047|  0.0369739| \*       |

</td>
<td>
|     | Effect                                  |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 5   | leg                                     |  0.7597056|  0.4008100|                |  0.8003634|  0.4056562|                |
| 6   | session.order:leg                       |  0.7597056|  0.8449885|                |  0.8003634|  0.8559635|                |
| 9   | stimulation:leg                         |  0.9937459|  0.5417737|                |  1.0831863|  0.5427260|                |
| 10  | session.order:stimulation:leg           |  0.9937459|  0.6827125|                |  1.0831863|  0.6839985|                |
| 13  | leg:direction                           |  0.9636560|  0.9226131|                |  1.0461526|  0.9281991|                |
| 14  | session.order:leg:direction             |  0.9636560|  0.1540577|                |  1.0461526|  0.1522160|                |
| 15  | stimulation:leg:direction               |  0.8004512|  0.7225967|                |  0.8487521|  0.7357368|                |
| 16  | session.order:stimulation:leg:direction |  0.8004512|  0.5582804|                |  0.8487521|  0.5685106|                |

</td>
</tr>
</tbody>
</table>
#### ANOVA matching Kanai et al. (2012) - center saccades

##### Without session order

Repeat the same ANOVA, but now discard the lateral and keep only center saccades (which Kanai did not have).

``` r
modelKanaiCenter <- ezANOVA(data = data.frame(filter(devDataBase, type == "center")),
                        dv = .(deviation.end), wid = .(subject), within = .(stimulation,leg,direction), type = 3)

kable(modelKanaiCenter)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect                    |  DFn|  DFd|           F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|-----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  10.3422560|  0.0035735| \*       |  0.0698321|
| 3   | leg                       |    2|   50|   2.4132849|  0.0998788|          |  0.0064777|
| 4   | direction                 |    1|   25|   0.4250384|  0.5203830|          |  0.0041044|
| 5   | stimulation:leg           |    2|   50|   0.6101220|  0.5472793|          |  0.0012937|
| 6   | stimulation:direction     |    1|   25|   2.7995910|  0.1067576|          |  0.0126867|
| 7   | leg:direction             |    2|   50|   3.8854661|  0.0270094| \*       |  0.0071117|
| 8   | stimulation:leg:direction |    2|   50|   0.5883554|  0.5590374|          |  0.0011173|

</td>
<td>
|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.8436804|  0.1300575|          |
| 5   | stimulation:leg           |  0.9603861|  0.6156733|          |
| 7   | leg:direction             |  0.6656519|  0.0075677| \*       |
| 8   | stimulation:leg:direction |  0.8210183|  0.0938067|          |

</td>
<td>
|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.8648128|  0.1084443|                |  0.9232826|  0.1046730|                |
| 5   | stimulation:leg           |  0.9618956|  0.5412918|                |  1.0404345|  0.5472793|                |
| 7   | leg:direction             |  0.7494296|  0.0403283| \*             |  0.7865648|  0.0380011| \*             |
| 8   | stimulation:leg:direction |  0.8481896|  0.5328387|                |  0.9034188|  0.5428502|                |

</td>
</tr>
</tbody>
</table>
###### Main effect of stimulation

``` r
devDataBase %>%
  filter(type == "center") %>%
  group_by(subject,stimulation) %>%
  summarise(deviation.end = mean(deviation.end)) %>%
  ggplot(aes(stimulation, deviation.end)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.data = mean_cl_normal, size = 1) +
  geom_jitter(width = 0.25)
```

<img src="accuracy_files/figure-markdown_github/Kanai-Center Main effect of stimulation-1.png" width="70%" style="display: block; margin: auto;" />

The accuracy in the cathodal session improves from baseline for most subjects; anodal stays the same or slightly worsens.

Let's do some follow-up tests to see whether the anodal or cathodal change scores are significantly different from 0 on their own.

Frequentist one-sample t-tests:

``` r
devDataBase %>%
  filter(type == "center") %>% # keep only center saccades
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(deviation.end = mean(deviation.end)) %>% # average over all other variables (df is now still grouped per stimulation)
  summarise_if(is.numeric, funs(list(tidy(t.test(.))))) %>%  # run one-sample t-test for each stimulation condition, return tidy data frames
  unnest() %>% # unpack the list-column with data frame for each test
  kable(.)
```

| stimulation |    estimate|   statistic|    p.value|  parameter|    conf.low|   conf.high| method            | alternative |
|:------------|-----------:|-----------:|----------:|----------:|-----------:|-----------:|:------------------|:------------|
| anodal      |   0.0223864|   0.9712741|  0.3407160|         25|  -0.0250828|   0.0698555| One Sample t-test | two.sided   |
| cathodal    |  -0.0755975|  -3.1803958|  0.0038987|         25|  -0.1245524|  -0.0266426| One Sample t-test | two.sided   |

The cathodal effect is highly signifcant, but the anodal is not.

###### Interaction: Leg by direction

``` r
devDataBase %>%
  filter(type == "center") %>%
  group_by(subject,leg,direction) %>%
  summarise(deviation.end = mean(deviation.end)) %>%
  ggplot(aes(leg, deviation.end, shape = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line", aes(group = direction, linetype = direction))
```

<img src="accuracy_files/figure-markdown_github/Kanai Interaction leg by direction-1.png" width="70%" style="display: block; margin: auto;" />

The accuracy of leftward saccades is improved for left saccades in the final half hour of task performance, more so than for right saccades, for which the accuracy goes back to baseline eventually.

##### With session order

``` r
modelKanaiCenterOrder <- ezANOVA(data = data.frame(filter(devDataBase, type == "center")), dv = .(deviation.end), 
          wid = .(subject), within = .(stimulation,leg,direction),  between = session.order, type = 3)
kable(modelKanaiCenterOrder)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect                                  |  DFn|  DFd|           F|          p| p&lt;.05 |        ges|
|-----|:----------------------------------------|----:|----:|-----------:|----------:|:---------|----------:|
| 2   | session.order                           |    1|   24|   0.9091324|  0.3498501|          |  0.0093825|
| 3   | stimulation                             |    1|   24|  11.2906179|  0.0025998| \*       |  0.0769854|
| 5   | leg                                     |    2|   48|   2.6985313|  0.0775127|          |  0.0074970|
| 7   | direction                               |    1|   24|   0.3758142|  0.5456181|          |  0.0039524|
| 4   | session.order:stimulation               |    1|   24|   1.7960865|  0.1927370|          |  0.0130944|
| 6   | session.order:leg                       |    2|   48|   1.3375601|  0.2720914|          |  0.0037301|
| 8   | session.order:direction                 |    1|   24|   0.1079995|  0.7452836|          |  0.0011390|
| 9   | stimulation:leg                         |    2|   48|   0.6721289|  0.5153601|          |  0.0014931|
| 11  | stimulation:direction                   |    1|   24|   3.6456189|  0.0682423|          |  0.0155529|
| 13  | leg:direction                           |    2|   48|   3.6359248|  0.0338605| \*       |  0.0072533|
| 10  | session.order:stimulation:leg           |    2|   48|   1.0528744|  0.3568522|          |  0.0023369|
| 12  | session.order:stimulation:direction     |    1|   24|   3.8042220|  0.0628957|          |  0.0162185|
| 14  | session.order:leg:direction             |    2|   48|   0.0826374|  0.9208157|          |  0.0001660|
| 15  | stimulation:leg:direction               |    2|   48|   0.4549433|  0.6371915|          |  0.0008799|
| 16  | session.order:stimulation:leg:direction |    2|   48|   1.7824344|  0.1791825|          |  0.0034386|

</td>
<td>
|     | Effect                                  |          W|          p| p&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------|
| 5   | leg                                     |  0.8326608|  0.1217259|          |
| 6   | session.order:leg                       |  0.8326608|  0.1217259|          |
| 9   | stimulation:leg                         |  0.9647617|  0.6619581|          |
| 10  | session.order:stimulation:leg           |  0.9647617|  0.6619581|          |
| 13  | leg:direction                           |  0.6671987|  0.0095265| \*       |
| 14  | session.order:leg:direction             |  0.6671987|  0.0095265| \*       |
| 15  | stimulation:leg:direction               |  0.8547089|  0.1644030|          |
| 16  | session.order:stimulation:leg:direction |  0.8547089|  0.1644030|          |

</td>
<td>
|     | Effect                                  |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 5   | leg                                     |  0.8566490|  0.0868332|                |  0.9160721|  0.0828527|                |
| 6   | session.order:leg                       |  0.8566490|  0.2711195|                |  0.9160721|  0.2716862|                |
| 9   | stimulation:leg                         |  0.9659611|  0.5105755|                |  1.0489826|  0.5153601|                |
| 10  | session.order:stimulation:leg           |  0.9659611|  0.3551055|                |  1.0489826|  0.3568522|                |
| 13  | leg:direction                           |  0.7502994|  0.0482991| \*             |  0.7892425|  0.0456985| \*             |
| 14  | session.order:leg:direction             |  0.7502994|  0.8692654|                |  0.7892425|  0.8793191|                |
| 15  | stimulation:leg:direction               |  0.8731405|  0.6112836|                |  0.9359565|  0.6245419|                |
| 16  | session.order:stimulation:leg:direction |  0.8731405|  0.1842757|                |  0.9359565|  0.1817943|                |

</td>
</tr>
</tbody>
</table>
There are some significant effects, but they do not interact with session order.

### Bayesian

*See the `median_latency.nb.html` notebook for more explanation of the Bayesian analyses*

#### Linear mixed effects matching Kanai - lateral saccades

Bayesian analogue of the frequentist repeated measures ANOVA (without order effect), with the same factors.

``` r
bfKanaiLateral = anovaBF(deviation.end~stimulation*leg*direction+subject, data = data.frame(filter(devDataBase, type == "lateral")), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiLateral = sort(bfKanaiLateral, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfKanaiLateral), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |         bf|
|-------------------------------------------------------------------------------------------------------------------------------|----------:|
| stimulation + subject                                                                                                         |  6.3193125|
| stimulation + direction + subject                                                                                             |  2.8777068|
| stimulation + direction + stimulation:direction + subject                                                                     |  0.5665412|
| direction + subject                                                                                                           |  0.4347378|
| stimulation + leg + subject                                                                                                   |  0.3496719|
| stimulation + direction + leg + subject                                                                                       |  0.1593865|
| leg + subject                                                                                                                 |  0.0537824|
| stimulation + direction + stimulation:direction + leg + subject                                                               |  0.0299003|
| stimulation + leg + stimulation:leg + subject                                                                                 |  0.0294942|
| direction + leg + subject                                                                                                     |  0.0233883|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |  0.0135322|
| stimulation + direction + leg + direction:leg + subject                                                                       |  0.0103230|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |  0.0025180|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |  0.0020019|
| direction + leg + direction:leg + subject                                                                                     |  0.0015891|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |  0.0008962|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |  0.0001833|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |  0.0000202|

Two models fare better than the null model: (1) a main effect of stimulation, and (2) a main effect of both stimulation and direction.

``` r
inclusionBF(bfKanaiLateral, models = "matched")
```

    ##                      effect Bayes.factor
    ## 1               stimulation   6.41983261
    ## 2                 direction   0.45261010
    ## 3     stimulation:direction   0.19633414
    ## 4                       leg   0.05501992
    ## 5           stimulation:leg   0.08457339
    ## 6             direction:leg   0.06555207
    ## 7 stimulation:direction:leg   0.11009814

There is moderate evidence for inclusion of an effect of stimulation, even though the classical analysis does not reach significance.

#### Linear mixed effects matching Kanai - center saccades

``` r
bfKanaiCenter = anovaBF(deviation.end~stimulation*leg*direction+subject, data = data.frame(filter(devDataBase, type == "center")), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiCenter = sort(bfKanaiCenter, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfKanaiCenter), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |            bf|
|-------------------------------------------------------------------------------------------------------------------------------|-------------:|
| stimulation + subject                                                                                                         |  4.059019e+04|
| stimulation + direction + stimulation:direction + subject                                                                     |  1.674061e+04|
| stimulation + direction + subject                                                                                             |  1.051873e+04|
| stimulation + leg + subject                                                                                                   |  4.207337e+03|
| stimulation + direction + stimulation:direction + leg + subject                                                               |  1.712539e+03|
| stimulation + direction + leg + subject                                                                                       |  1.079519e+03|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |  3.602672e+02|
| stimulation + leg + stimulation:leg + subject                                                                                 |  3.254140e+02|
| stimulation + direction + leg + direction:leg + subject                                                                       |  2.156527e+02|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |  1.275241e+02|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |  8.579462e+01|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |  2.707515e+01|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |  1.735635e+01|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |  3.278777e+00|
| direction + subject                                                                                                           |  2.396816e-01|
| leg + subject                                                                                                                 |  9.409320e-02|
| direction + leg + subject                                                                                                     |  2.307420e-02|
| direction + leg + direction:leg + subject                                                                                     |  4.186700e-03|

All the models with a main effect of stimulation are strongly supported.

``` r
inclusionBF(bfKanaiCenter, models = "matched")
```

    ##                      effect Bayes.factor
    ## 1               stimulation 4.159438e+04
    ## 2                 direction 2.589376e-01
    ## 3     stimulation:direction 1.591669e+00
    ## 4                       leg 1.031604e-01
    ## 5           stimulation:leg 7.698218e-02
    ## 6             direction:leg 2.064136e-01
    ## 7 stimulation:direction:leg 1.210991e-01

Overwhelming evidence for inclusion of a main effect of stimulation, which is in accord with the highly significant p-value.

Bayesian one-sample t-tests:

``` r
devDataBase %>% 
  filter(type == "center") %>% # keep only center saccades
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(deviation.end = mean(deviation.end)) %>% # average over all other variables
  spread(stimulation,deviation.end) %>% # make separate columns with test data
  summarise_if(is.numeric, funs(extractBF(ttestBF(.), onlybf = TRUE))) %>% # run Bayesian t-test on each column, keeping only the BF
  gather(stimulation,BF,anodal,cathodal) %>% # make row for each stimulation condition
  kable(.)
```

| stimulation |          BF|
|:------------|-----------:|
| anodal      |   0.3173718|
| cathodal    |  10.5264920|

The cathodal effect on its ownhas a BF<sub>10</sub> &gt; 10, but the anodal effect does not.

Saccade end point variability
=============================

In the motor literature, people often look at the spread in movement endpoints, as it's often believed that this is what the motor system is trying to optimize (i.e. minimize). [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045) also examined this, but found no effects of tDCS.

Calculate endpoint variability
------------------------------

Kanai et al. (2012) operationalized variability as the standard deviation of the x-coordinate of the saccade end point.

``` r
stdData <- groupData %>%
  group_by(subject,stimulation,leg,direction,type) %>% 
  summarise(std.deviation.x = sd(deviation.end.x))
```

This is a summary measure across trials, so we have one estimate per subject per condition:

``` r
kable(head(stdData))
```

| subject | stimulation | leg  | direction | type    |  std.deviation.x|
|:--------|:------------|:-----|:----------|:--------|----------------:|
| S01     | anodal      | pre  | left      | lateral |        0.6287700|
| S01     | anodal      | pre  | left      | center  |        0.6782264|
| S01     | anodal      | pre  | right     | lateral |        0.5883124|
| S01     | anodal      | pre  | right     | center  |        0.4780784|
| S01     | anodal      | tDCS | left      | lateral |        0.5671711|
| S01     | anodal      | tDCS | left      | center  |        0.6352989|

Prepare data frame for plotting & statistics
--------------------------------------------

``` r
stdData$leg <- fct_recode(stdData$leg, baseline = "pre") # recode factor to match deviation data frame
```

Subtract the baseline from each average:

``` r
# Subtract baseline
stdDataBase <- stdData %>%
  group_by(subject,stimulation,direction,type) %>% 
  # subtract baseline block from others, make new column
  summarise(tDCS = std.deviation.x[leg == "tDCS"] - std.deviation.x[leg == "baseline"], 
            post.1 = std.deviation.x[leg == "post.1"] - std.deviation.x[leg == "baseline"],
            post.2 = std.deviation.x[leg == "post.2"] - std.deviation.x[leg == "baseline"]) %>%
gather(leg, std.deviation.x, tDCS, post.1, post.2) %>% # gather new columns to use as factor 
mutate(leg = factor(leg, levels = c("tDCS", "post.1", "post.2"))) # reorder factor levels
```

Plot
----

### With baseline block

``` r
kanaiPlotStd <- ggplot(stdData, aes(leg, std.deviation.x, color = stimulation, shape = stimulation)) +         
  facet_grid(type ~ direction) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3) +
  ggtitle("Horizontal standard deviation")
kanaiPlotStd
```

<img src="accuracy_files/figure-markdown_github/Line plot per leg - variability-1.png" width="70%" style="display: block; margin: auto;" />

The changes here seem less pronounced than in the endpoint deviation data.

Let's look at the individual subject data:

``` r
kanaiPlotSubsAnodal <- ggplot(stdData[stdData$stimulation == "anodal", ], aes(leg, std.deviation.x)) +
  facet_grid(type ~ direction) +
  geom_line(aes(group = subject,color = subject)) +
  stat_summary(fun.y = mean, aes(group = stimulation), geom = "line") +
  stat_summary(fun.y = mean, geom = "point") +
  ggtitle("Anodal session")
kanaiPlotSubsAnodal
```

<img src="accuracy_files/figure-markdown_github/Line plot per subject - variability anodal-1.png" width="70%" style="display: block; margin: auto;" />

``` r
kanaiPlotSubsCathodal <- ggplot(stdData[stdData$stimulation == "cathodal", ], aes(leg, std.deviation.x)) +
  facet_grid(type ~ direction) +
  geom_line(aes(group = subject,color = subject)) +
  stat_summary(fun.y = mean, aes(group = stimulation), geom = "line") +
  stat_summary(fun.y = mean, geom = "point") +
  ggtitle("Cathodal session")
kanaiPlotSubsCathodal
```

<img src="accuracy_files/figure-markdown_github/Line plot per subject - variability cathodal-1.png" width="70%" style="display: block; margin: auto;" />

This measure seems particularly variable across subjects and also subject to quite a few spikes that only show up in a few conditions.

#### Baseline reliability

Scatterplot and correlation of baseline data in the two sessions:

``` r
baselineCorrStd <- stdData %>%
  filter(leg == "baseline") %>% 
  group_by(direction,type) %>% 
  spread(stimulation,std.deviation.x) %>% 
  nest() %>% 
  mutate(stats = map(data, ~cor.test(formula = ~ anodal + cathodal, data =.))) %>% # run correlation test on baselines from each condition
  mutate(tidy_model = map(stats, tidy)) %>% 
  unnest(tidy_model, .drop = TRUE)
```

``` r
stdData %>%
  filter(leg == "baseline") %>%
  spread(stimulation,std.deviation.x) %>%
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column on session order from other data frame
  ggplot(aes(anodal,cathodal)) +
    facet_grid(type ~ direction) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_smooth(method = "lm") +
    geom_point(aes(color=session.order)) +
    xlim(0.2,1.6) + ylim(0.2,1.6) +
    geom_text(data = baselineCorrStd, x = 0.5, y = 1.5, aes(label = paste("italic(r) == ", round(estimate,2))), parse = TRUE) +
    labs(title = "Baseline in anodal and cathodal sessions", subtitle = "scatterplot of baseline endpoint variability")
```

<img src="accuracy_files/figure-markdown_github/baseline scatterplots - variability-1.png" width="70%" style="display: block; margin: auto;" />

The correlations for this measure are quite low, especially for the center conditions. Apparently the standard deviation of saccade endpoints is not so reliable.

#### Baseline differences

``` r
stdData %>%
  filter(leg == "baseline") %>%
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column on session order from other data frame
  group_by(subject,direction,type,session.order) %>%
  summarise(std.deviation.x.diff = std.deviation.x[stimulation == "anodal"] - std.deviation.x[stimulation == "cathodal"]) %>%
  ggplot(aes(factor(0), std.deviation.x.diff)) +
    facet_grid(type ~ direction) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_summary(fun.data = mean_cl_normal) +
    stat_summary(fun.y = mean, aes(label=round(..y.., digits=2), x = 1.3), geom = "label", alpha = 0.5) +
    geom_point(shape = 21, aes(colour = session.order), position = position_jitter(width=.1)) +
    labs(title = "Baseline in anodal and cathodal sessions", subtitle = "anodal - cathodal")
```

<img src="accuracy_files/figure-markdown_github/baseline difference stripchart - variability-1.png" width="70%" style="display: block; margin: auto;" />

The average baseline differences are small, but the spread is quite large.

``` r
stdData %>%
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

On average none of the baselines differ significantly from each other.

### Baseline subtracted

``` r
kanaiPlotStdBase <- ggplot(stdDataBase, aes(leg, std.deviation.x, color = stimulation, shape = stimulation)) +         
  facet_grid(type ~ direction) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3)
kanaiPlotStdBase
```

<img src="accuracy_files/figure-markdown_github/Line plot from baseline - standard deviation-1.png" width="70%" style="display: block; margin: auto;" />

Here the changes are even tinier than the endpoint deviation data (&lt;.1 degree). If anything, the differences seem to grow more pronounced after tDCS.

Statistics
----------

``` r
# Make "subject" a factor, so we can model the repeated measures
stdDataBase <- stdDataBase %>%
  ungroup() %>% # remove any grouping info, because we need to refactor
  inner_join(., subjectData[ ,c("subject","session.order")], by = c("subject")) %>% # add column on session order from other data frame
  mutate(subject = factor(subject)) # refactor
```

### Frequentist

#### ANOVA matching Kanai et al. (2012) - lateral saccades

##### Without session order

**Data**:

-   Outliers removed
-   Collapsed into 15-minute intervals
-   Subtract the baseline from each subsequent block
-   Discard center, keep only lateral saccades

**Dependent measure**: saccade end point variability (horizontal standard deviation)

**Factors**:

-   STIMULATION (anodal vs. cathodal)
-   LEG (tDCS, post.1, post.2)
-   DIRECTION (left vs. right)

``` r
modelKanaiStd <- ezANOVA(data = data.frame(filter(stdDataBase, type == "lateral")),
                        dv = .(std.deviation.x), wid = .(subject), within = .(stimulation,leg,direction), type = 3)

kable(modelKanaiStd)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect                    |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  1.2189602|  0.2800786|          |  0.0138324|
| 3   | leg                       |    2|   50|  0.4518157|  0.6390443|          |  0.0011604|
| 4   | direction                 |    1|   25|  0.2471961|  0.6234009|          |  0.0010690|
| 5   | stimulation:leg           |    2|   50|  1.1183082|  0.3348684|          |  0.0030580|
| 6   | stimulation:direction     |    1|   25|  0.1204196|  0.7314841|          |  0.0004919|
| 7   | leg:direction             |    2|   50|  0.9403214|  0.3972980|          |  0.0023341|
| 8   | stimulation:leg:direction |    2|   50|  0.1873223|  0.8297557|          |  0.0003276|

</td>
<td>
|     | Effect                    |         W|          p| p&lt;.05 |
|-----|:--------------------------|---------:|----------:|:---------|
| 3   | leg                       |  0.994425|  0.9351138|          |
| 5   | stimulation:leg           |  0.994966|  0.9412372|          |
| 7   | leg:direction             |  0.908637|  0.3167270|          |
| 8   | stimulation:leg:direction |  0.903796|  0.2970605|          |

</td>
<td>
|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.9944559|  0.6379765|                |  1.0801686|  0.6390443|                |
| 5   | stimulation:leg           |  0.9949912|  0.3346792|                |  1.0808237|  0.3348684|                |
| 7   | leg:direction             |  0.9162854|  0.3907925|                |  0.9851512|  0.3961957|                |
| 8   | stimulation:leg:direction |  0.9122389|  0.8102650|                |  0.9802676|  0.8255939|                |

</td>
</tr>
</tbody>
</table>
##### With session order

``` r
modelKanaiStdOrder <- ezANOVA(data = data.frame(filter(stdDataBase, type == "lateral")), dv = .(std.deviation.x), 
          wid = .(subject), within = .(stimulation,leg,direction),  between = session.order, type = 3)
kable(modelKanaiStdOrder)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect                                  |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:----------------------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | session.order                           |    1|   24|  0.1607364|  0.6920265|          |  0.0017906|
| 3   | stimulation                             |    1|   24|  1.1906769|  0.2860271|          |  0.0143992|
| 5   | leg                                     |    2|   48|  0.4813409|  0.6209047|          |  0.0013051|
| 7   | direction                               |    1|   24|  0.2831863|  0.5995127|          |  0.0012899|
| 4   | session.order:stimulation               |    1|   24|  0.0246391|  0.8765826|          |  0.0003022|
| 6   | session.order:leg                       |    2|   48|  0.2627471|  0.7700350|          |  0.0007128|
| 8   | session.order:direction                 |    1|   24|  0.3163115|  0.5790487|          |  0.0014406|
| 9   | stimulation:leg                         |    2|   48|  1.1029347|  0.3401590|          |  0.0032123|
| 11  | stimulation:direction                   |    1|   24|  0.2860170|  0.5977048|          |  0.0010464|
| 13  | leg:direction                           |    2|   48|  0.8050343|  0.4530163|          |  0.0020275|
| 10  | session.order:stimulation:leg           |    2|   48|  0.0483124|  0.9528823|          |  0.0001411|
| 12  | session.order:stimulation:direction     |    1|   24|  4.5903060|  0.0425038| \*       |  0.0165327|
| 14  | session.order:leg:direction             |    2|   48|  1.2562264|  0.2939184|          |  0.0031603|
| 15  | stimulation:leg:direction               |    2|   48|  0.2086854|  0.8123830|          |  0.0003866|
| 16  | session.order:stimulation:leg:direction |    2|   48|  0.1838631|  0.8326329|          |  0.0003406|

</td>
<td>
|     | Effect                                  |          W|          p| p&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------|
| 5   | leg                                     |  0.9929565|  0.9219294|          |
| 6   | session.order:leg                       |  0.9929565|  0.9219294|          |
| 9   | stimulation:leg                         |  0.9948738|  0.9426101|          |
| 10  | session.order:stimulation:leg           |  0.9948738|  0.9426101|          |
| 13  | leg:direction                           |  0.8816782|  0.2349990|          |
| 14  | session.order:leg:direction             |  0.8816782|  0.2349990|          |
| 15  | stimulation:leg:direction               |  0.9042201|  0.3141614|          |
| 16  | session.order:stimulation:leg:direction |  0.9042201|  0.3141614|          |

</td>
<td>
|     | Effect                                  |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 5   | leg                                     |  0.9930058|  0.6196034|                |  1.0822730|  0.6209047|                |
| 6   | session.order:leg                       |  0.9930058|  0.7685199|                |  1.0822730|  0.7700350|                |
| 9   | stimulation:leg                         |  0.9949000|  0.3399524|                |  1.0846107|  0.3401590|                |
| 10  | session.order:stimulation:leg           |  0.9949000|  0.9522849|                |  1.0846107|  0.9528823|                |
| 13  | leg:direction                           |  0.8941970|  0.4412480|                |  0.9614309|  0.4488867|                |
| 14  | session.order:leg:direction             |  0.8941970|  0.2920372|                |  0.9614309|  0.2933216|                |
| 15  | stimulation:leg:direction               |  0.9125920|  0.7926861|                |  0.9837647|  0.8089070|                |
| 16  | session.order:stimulation:leg:direction |  0.9125920|  0.8133016|                |  0.9837647|  0.8292325|                |

</td>
</tr>
</tbody>
</table>
The interaction with session order, stimulation, and direction is significant. However, the stimulation:direction interaction was not significant in the ANOVA without the session order factor, so we should interpret this with caution. In addition, an interaction of session order and stimulation could just as well reflect a main effect of session (1 vs. 2): there's no way to distinguish between these possibilities.

#### ANOVA matching Kanai et al. (2012) - center saccades

##### Without session order

``` r
modelKanaiStdCenter <- ezANOVA(data = data.frame(filter(stdDataBase, type == "center")),
                        dv = .(std.deviation.x), wid = .(subject), within = .(stimulation,leg,direction), type = 3)

kable(modelKanaiStdCenter)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect                    |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:--------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | stimulation               |    1|   25|  3.8851203|  0.0598820|          |  0.0399400|
| 3   | leg                       |    2|   50|  2.7073528|  0.0764935|          |  0.0075993|
| 4   | direction                 |    1|   25|  0.6021177|  0.4450504|          |  0.0033319|
| 5   | stimulation:leg           |    2|   50|  1.1805422|  0.3155254|          |  0.0036640|
| 6   | stimulation:direction     |    1|   25|  0.1745759|  0.6796436|          |  0.0010192|
| 7   | leg:direction             |    2|   50|  2.1601871|  0.1259451|          |  0.0035285|
| 8   | stimulation:leg:direction |    2|   50|  0.4735067|  0.6255786|          |  0.0007475|

</td>
<td>
|     | Effect                    |          W|          p| p&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------|
| 3   | leg                       |  0.8960127|  0.2677749|          |
| 5   | stimulation:leg           |  0.8569616|  0.1568686|          |
| 7   | leg:direction             |  0.6679634|  0.0078892| \*       |
| 8   | stimulation:leg:direction |  0.8089994|  0.0785921|          |

</td>
<td>
|     | Effect                    |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:--------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 3   | leg                       |  0.9058075|  0.0824883|                |  0.9725125|  0.0781996|                |
| 5   | stimulation:leg           |  0.8748612|  0.3116036|                |  0.9353175|  0.3136546|                |
| 7   | leg:direction             |  0.7507301|  0.1405427|                |  0.7880908|  0.1383445|                |
| 8   | stimulation:leg:direction |  0.8396302|  0.5929099|                |  0.8932129|  0.6044505|                |

</td>
</tr>
</tbody>
</table>
##### Main effect of stimulation

This effect is just non-significant, but let's inspect anyway:

``` r
stdDataBase %>%
  filter(type == "center") %>%
  group_by(subject,stimulation) %>%
  summarise(std.deviation.x = mean(std.deviation.x)) %>%
  ggplot(aes(stimulation, std.deviation.x)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.data = mean_cl_normal, size = 1) +
  geom_jitter(width = 0.25)
```

<img src="accuracy_files/figure-markdown_github/Kanai-Center-variability Main effect of stimulation-1.png" width="70%" style="display: block; margin: auto;" />

This resembles the difference found for the saccade endpoint deviation, except here the larger and more consistent effect seems to be in the anodal condition.

##### With session order

``` r
modelKanaiStdCenterOrder <- ezANOVA(data = data.frame(filter(stdDataBase, type == "center")), dv = .(std.deviation.x), 
          wid = .(subject), within = .(stimulation,leg,direction),  between = session.order, type = 3)
kable(modelKanaiStdCenterOrder)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect                                  |  DFn|  DFd|          F|          p| p&lt;.05 |        ges|
|-----|:----------------------------------------|----:|----:|----------:|----------:|:---------|----------:|
| 2   | session.order                           |    1|   24|  0.1624780|  0.6904536|          |  0.0015282|
| 3   | stimulation                             |    1|   24|  4.3268347|  0.0483642| \*       |  0.0451874|
| 5   | leg                                     |    2|   48|  2.8313514|  0.0688095|          |  0.0084055|
| 7   | direction                               |    1|   24|  0.8197631|  0.3742496|          |  0.0045103|
| 4   | session.order:stimulation               |    1|   24|  1.5157716|  0.2301895|          |  0.0163088|
| 6   | session.order:leg                       |    2|   48|  0.6228833|  0.5406747|          |  0.0018614|
| 8   | session.order:direction                 |    1|   24|  2.1829361|  0.1525517|          |  0.0119209|
| 9   | stimulation:leg                         |    2|   48|  1.1504478|  0.3250667|          |  0.0038706|
| 11  | stimulation:direction                   |    1|   24|  0.2458431|  0.6245275|          |  0.0014926|
| 13  | leg:direction                           |    2|   48|  2.1991247|  0.1219515|          |  0.0036576|
| 10  | session.order:stimulation:leg           |    2|   48|  0.0383245|  0.9624300|          |  0.0001294|
| 12  | session.order:stimulation:direction     |    1|   24|  1.0500393|  0.3157162|          |  0.0063444|
| 14  | session.order:leg:direction             |    2|   48|  1.5933477|  0.2138058|          |  0.0026527|
| 15  | stimulation:leg:direction               |    2|   48|  0.4598347|  0.6341403|          |  0.0007630|
| 16  | session.order:stimulation:leg:direction |    2|   48|  0.7980825|  0.4560741|          |  0.0013235|

</td>
<td>
|     | Effect                                  |          W|          p| p&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------|
| 5   | leg                                     |  0.8967279|  0.2854944|          |
| 6   | session.order:leg                       |  0.8967279|  0.2854944|          |
| 9   | stimulation:leg                         |  0.8553298|  0.1657816|          |
| 10  | session.order:stimulation:leg           |  0.8553298|  0.1657816|          |
| 13  | leg:direction                           |  0.6851706|  0.0129326| \*       |
| 14  | session.order:leg:direction             |  0.6851706|  0.0129326| \*       |
| 15  | stimulation:leg:direction               |  0.7903814|  0.0668532|          |
| 16  | session.order:stimulation:leg:direction |  0.7903814|  0.0668532|          |

</td>
<td>
|     | Effect                                  |        GGe|    p\[GG\]| p\[GG\]&lt;.05 |        HFe|    p\[HF\]| p\[HF\]&lt;.05 |
|-----|:----------------------------------------|----------:|----------:|:---------------|----------:|----------:|:---------------|
| 5   | leg                                     |  0.9063947|  0.0747228|                |  0.9762321|  0.0702670|                |
| 6   | session.order:leg                       |  0.9063947|  0.5258855|                |  0.9762321|  0.5370486|                |
| 9   | stimulation:leg                         |  0.8736141|  0.3204306|                |  0.9365284|  0.3229043|                |
| 10  | session.order:stimulation:leg           |  0.8736141|  0.9470365|                |  0.9365284|  0.9553929|                |
| 13  | leg:direction                           |  0.7605549|  0.1363083|                |  0.8013684|  0.1338319|                |
| 14  | session.order:leg:direction             |  0.7605549|  0.2189988|                |  0.8013684|  0.2183303|                |
| 15  | stimulation:leg:direction               |  0.8267068|  0.5980305|                |  0.8801197|  0.6098886|                |
| 16  | session.order:stimulation:leg:direction |  0.8267068|  0.4356829|                |  0.8801197|  0.4424174|                |

</td>
</tr>
</tbody>
</table>
Here the effect does just reach significance, but there's no interaction with session order.

### Bayesian

Bayesian analogues of the frequentist repeated measures ANOVAs (without order effect), with the same factors.

#### Linear mixed effects matching Kanai - lateral saccades

``` r
bfKanaiStd = anovaBF(std.deviation.x~stimulation*leg*direction+subject, data = data.frame(filter(stdDataBase, type == "lateral")), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiStd = sort(bfKanaiStd, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfKanaiStd), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |         bf|
|-------------------------------------------------------------------------------------------------------------------------------|----------:|
| stimulation + subject                                                                                                         |  1.6435119|
| stimulation + direction + subject                                                                                             |  0.2506199|
| direction + subject                                                                                                           |  0.1522912|
| stimulation + leg + subject                                                                                                   |  0.0688916|
| stimulation + direction + stimulation:direction + subject                                                                     |  0.0461304|
| leg + subject                                                                                                                 |  0.0423792|
| stimulation + direction + leg + subject                                                                                       |  0.0104902|
| stimulation + leg + stimulation:leg + subject                                                                                 |  0.0075101|
| direction + leg + subject                                                                                                     |  0.0063375|
| stimulation + direction + stimulation:direction + leg + subject                                                               |  0.0019310|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |  0.0010857|
| stimulation + direction + leg + direction:leg + subject                                                                       |  0.0010284|
| direction + leg + direction:leg + subject                                                                                     |  0.0006035|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |  0.0001965|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |  0.0001836|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |  0.0001012|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |  0.0000194|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |  0.0000021|

``` r
# Inclusion Bayes Factors
inclusionBF(bfKanaiStd, models = "matched")
```

    ##                      effect Bayes.factor
    ## 1               stimulation   1.64324522
    ## 2                 direction   0.15234609
    ## 3     stimulation:direction   0.18403395
    ## 4                       leg   0.04204598
    ## 5           stimulation:leg   0.10800287
    ## 6             direction:leg   0.09660458
    ## 7 stimulation:direction:leg   0.10831713

Across the board, there is only marginal support for an effect of stimulation. For the interaction between stimulation and direction, the BF approaches moderate evidence for the null.

#### Linear mixed effects matching Kanai - center saccades

``` r
bfKanaiStdCenter = anovaBF(std.deviation.x~stimulation*leg*direction+subject, data = data.frame(filter(stdDataBase, type == "center")), whichModels = "withmain", whichRandom = "subject", progress = FALSE, iterations = 100000) # compute Bayes Factors
bfKanaiStdCenter = sort(bfKanaiStdCenter, decreasing = TRUE) # sort such that winning model is at the top
```

``` r
kable(select(extractBF(bfKanaiStdCenter), bf)) # show only the Bayes factors in a table
```

|                                                                                                                               |           bf|
|-------------------------------------------------------------------------------------------------------------------------------|------------:|
| stimulation + subject                                                                                                         |  141.3218729|
| stimulation + direction + subject                                                                                             |   31.3722123|
| stimulation + leg + subject                                                                                                   |   17.4798589|
| stimulation + direction + stimulation:direction + subject                                                                     |    6.2920364|
| stimulation + direction + leg + subject                                                                                       |    3.8156410|
| stimulation + leg + stimulation:leg + subject                                                                                 |    2.0747730|
| stimulation + direction + stimulation:direction + leg + subject                                                               |    0.8217815|
| stimulation + direction + leg + stimulation:leg + subject                                                                     |    0.4303726|
| stimulation + direction + leg + direction:leg + subject                                                                       |    0.4166426|
| direction + subject                                                                                                           |    0.2142540|
| leg + subject                                                                                                                 |    0.1155812|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + subject                                             |    0.0874919|
| stimulation + direction + stimulation:direction + leg + direction:leg + subject                                               |    0.0843555|
| stimulation + direction + leg + stimulation:leg + direction:leg + subject                                                     |    0.0674698|
| direction + leg + subject                                                                                                     |    0.0247495|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + subject                             |    0.0100521|
| direction + leg + direction:leg + subject                                                                                     |    0.0026778|
| stimulation + direction + stimulation:direction + leg + stimulation:leg + direction:leg + stimulation:direction:leg + subject |    0.0011877|

Like for the saccade endpoint deviation data, models with stimulation as a factor receive some support, although to a less strong degree. In contrast to endpoint deviation though, here the classical analysis was (barely) non-significant, so there is a discrepancy between the Bayesian and Frequentist approaches.

``` r
# Inclusion Bayes Factors
inclusionBF(bfKanaiStdCenter, models = "matched")
```

    ##                      effect Bayes.factor
    ## 1               stimulation  143.2340660
    ## 2                 direction    0.2213517
    ## 3     stimulation:direction    0.2020843
    ## 4                       leg    0.1235159
    ## 5           stimulation:leg    0.1180532
    ## 6             direction:leg    0.1121996
    ## 7 stimulation:direction:leg    0.1181513

Again, especially considering the non-significant p-value, the support is quite strong.

Let's do some follow-up tests to see whether the anodal or cathodal change scores are significantly different from 0 on their own.

Bayesian one-sample t-tests:

``` r
stdDataBase %>%
  filter(type == "center") %>% # keep only center saccades
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(deviation.end = mean(std.deviation.x)) %>% # average over all other variables
  spread(stimulation,deviation.end) %>% # make separate columns with test data
  summarise_if(is.numeric, funs(extractBF(ttestBF(.), onlybf = TRUE))) %>% # run Bayesian t-test on each column, keeping only the BF
  gather(stimulation,BF,anodal,cathodal) %>% # make row for each stimulation condition
  kable(.)
```

| stimulation |         BF|
|:------------|----------:|
| anodal      |  0.7138830|
| cathodal    |  0.4103342|

Frequentist one-sample t-tests:

``` r
stdDataBase %>%
  filter(type == "center") %>% # keep only center saccades
  group_by(stimulation,subject) %>% # for each session and subject
  summarise(deviation.end = mean(std.deviation.x)) %>% # average over all other variables (df is now still grouped per stimulation)
  summarise_if(is.numeric, funs(list(tidy(t.test(.))))) %>%  # run one-sample t-test for each stimulation condition, return tidy data frames
  unnest() %>% # unpack the list-column with data frame for each test
  kable(.) 
```

| stimulation |    estimate|  statistic|    p.value|  parameter|    conf.low|  conf.high| method            | alternative |
|:------------|-----------:|----------:|----------:|----------:|-----------:|----------:|:------------------|:------------|
| anodal      |   0.0477794|   1.683151|  0.1047947|         25|  -0.0106845|  0.1062434| One Sample t-test | two.sided   |
| cathodal    |  -0.0340719|  -1.236236|  0.2278603|         25|  -0.0908350|  0.0226911| One Sample t-test | two.sided   |

So neither effect holds up on their own.
