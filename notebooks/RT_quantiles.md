sacc-tDCS: Saccade latency quantile analysis
================
Leon Reteig

-   [Setup](#setup)
-   [Shift function analysis](#shift-function-analysis)
    -   [Estimate quantiles and compute shift function](#estimate-quantiles-and-compute-shift-function)
    -   [Summary plots](#summary-plots)
    -   [Individual shift functions](#individual-shift-functions)

Setup
=====

R notebook for analysis of saccade latency quantiles in the `sacc-tDCS` dataset.

``` r
# Load some libraries
library(here) # file paths
```

    ## here() starts at /Volumes/research$/reteig/sacc-tDCS

``` r
library(tidyverse) # ggplot2, tibble, tidyr, readr, purrr, and dplyr
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
library(rogme) # robust graphical methods (quantile estimation, shift function)
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
    ##  [1] knitr_1.15.1     rogme_0.1.0.9000 dplyr_0.5.0      purrr_0.2.2     
    ##  [5] readr_1.1.0      tidyr_0.6.1      tibble_1.3.0     ggplot2_2.2.1   
    ##  [9] tidyverse_1.1.1  here_0.1        
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
# Load the data frame
groupData <- mutate(groupData,
                    leg = as.character(leg), # cannot edit leg if it's still a factor
                    leg = replace(leg, leg == "pre", "baseline"),
                    leg = replace(leg, leg == "post", "post-1"),
                    leg = replace(leg, block > 3, "post-2"),
                    leg = factor(leg, levels = c("baseline", "tDCS", "post-1", "post-2")) # refactor and order levels
                    )
```

``` r
kable(head(groupData))
```

| subject | stimulation | leg      |  block|  trial| type    | direction |  deviation.start|  deviation.end.x|  deviation.end.y|  amplitude|  latency|    drift.x|   drift.y|
|:--------|:------------|:---------|------:|------:|:--------|:----------|----------------:|----------------:|----------------:|----------:|--------:|----------:|---------:|
| S01     | anodal      | baseline |      1|      1| lateral | right     |         0.462897|         0.170455|       -0.0080638|    8.02463|      433|  0.0953736|  0.102814|
| S01     | anodal      | baseline |      1|      1| center  | left      |         0.459092|         1.032850|        0.0665268|    7.16262|      439|  0.0953736|  0.102814|
| S01     | anodal      | baseline |      1|      2| lateral | right     |         0.344561|        -0.344967|        0.2197400|    7.74873|      291|  0.0953736|  0.102814|
| S01     | anodal      | baseline |      1|      2| center  | left      |         0.550230|         0.361201|        0.3507760|    7.43233|      198|  0.0953736|  0.102814|
| S01     | anodal      | baseline |      1|      3| lateral | right     |         0.514736|        -0.588470|        0.1673250|    7.61080|      281|  0.0953736|  0.102814|
| S01     | anodal      | baseline |      1|      3| center  | left      |         0.620728|         1.576610|        0.4031910|    6.35043|      376|  0.0953736|  0.102814|

-   **subject**: subject ID
-   **stimulation**: Whether data are from the anodal or cathodal session
-   **leg**: Whether data are before (`pre`), during (`tDCS`), or after (`post.1`, `post.2`) tDCS
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

Shift function analysis
=======================

This analysis is based on Rand Wilcox's shift function (for dependent groups), as published in [Wilcox, R. R., & Erceg-Hurn, D. M. (2012). Comparing two dependent groups via quantiles. Journal of Applied Statistics, 39(12), 2655-2664](http://dx.doi.org/10.1080/02664763.2012.724665). It uses the implementation in the [`rogme` package](https://github.com/GRousselet/rogme). See the [accompanying paper and code](https://figshare.com/articles/Modern_graphical_methods_to_compare_two_groups_of_observations/4055970) for much more info.

Estimate quantiles and compute shift function
---------------------------------------------

Briefly, the code below performs the following steps:

-   Compute the deciles (10% quantiles) for each individual subject's distribution of saccade latencies, for each condition, using the [Harrell-Davis quantile estimator](https://garstats.wordpress.com/2016/06/09/the-harrell-davis-quantile-estimator/). This procedure uses a weighted average of order statistics (based on the data points) to estimate each quantile, so it tends to be more precise than simple methods, and it is said to perform well in a wide variety of cases.
-   Take two of these distributions (in this case, the *anodal* and *cathodal* sessions) and subtract the deciles to look at the difference. A plot of the decile differences vs. the deciles in one condition is called a [shift function](https://garstats.wordpress.com/2016/07/12/shift-function/).
-   Compute 95% confidence intervals with a percentile bootstrap of the decile differences. If the confidence interval of a particular decile difference does not include zero, the difference is significant. The confidence intervals are **not** corrected for multiple comparisons, but the p-values are (using Hochberg's method).

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

This takes quite a while to compute with a large number of bootstrap samples, so write the data file to disk. The current and next code chunk are not evaluated by default! (`eval=FALSE`)

``` r
write_csv(qData, here("data", "sacc-tDCS_quantiles.csv"))
```

``` r
subs2exclude <- subs2exclude <- c("S21","S25","S16","S22","S28") 
qData <- read_csv(here("data", "sacc-tDCS_quantiles.csv")) %>%
  filter(!(subject %in% subs2exclude)) %>% # exclude subjects
    mutate(leg = replace(leg, leg == "pre", "baseline"), # rename and reoder levels
         leg = replace(leg, leg == "post.1", "post-1"),
         leg = replace(leg, leg == "post.2", "post-2"),
         leg = factor(leg, levels = c("baseline", "tDCS", "post-1", "post-2")))
```

Summary plots
-------------

Add decile codes and anodal medians to data frame, for plotting:

``` r
qStats <- qData %>%
  group_by(subject,leg,type,direction) %>%
  mutate(deco = c(seq(1,5),seq(4,1))) %>% # add code of deciles to data frame
  group_by(leg,type,direction,q) %>% 
  mutate(anodal = mean(anodal)) %>% # mean of "anodal" quantiles OVER subjects
  group_by(leg,type,direction) %>% 
  mutate(anodal_median = median(anodal)) # median for plotting
```

For each quantile, count which subjects show significant effects, and in which direction:

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

### Lateral saccades

Shift functions:

``` r
ggplot(filter(qStats, type == "lateral"), aes(anodal, difference)) +
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

<img src="RT_quantiles_files/figure-markdown_github/Shift function lateral-1.png" width="75%" style="display: block; margin: auto;" />

It looks like the slowest saccades show the biggest difference for most conditions. Curiously, cathodal tDCS there leads to faster latencies. However, this pattern already seems to be present in the baseline.

Let's examine the number of subjects with a signifcant difference per quantile:

``` r
ggplot(filter(qSig, type == "lateral", !is.na(significance)), aes(q, fill = significance)) +
  facet_grid(leg ~ direction) +
  geom_bar(position = "stack") +
  stat_bin(binwidth = .1, geom = "text", size = 2, aes(label = ..count..), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_continuous("decile", breaks = seq(0.1,0.9,0.1), labels = seq(1,9,1)) +
  scale_y_continuous("number of participants with significant difference", limits = c(0, length(unique(qSig$subject))), breaks = c(0,10,20,length(unique(qSig$subject))))
```

<img src="RT_quantiles_files/figure-markdown_github/Significance shift function lateral-1.png" width="75%" style="display: block; margin: auto;" />

This effect in the right tail seems to be driven by a relatively small number of participants, and a sizable number of participants also shows an anodal effect there... Further, if anyhting, the most significant differences are in the baseline block...

### Center saccades

``` r
ggplot(filter(qStats, type == "center"), aes(anodal, difference)) +
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

<img src="RT_quantiles_files/figure-markdown_github/Shift function center-1.png" width="75%" style="display: block; margin: auto;" />

There the pattern appears opposite: the greatest difference is in the fastest saccades. It makes sense that these would be most impacted for center saccades. They could be impulsive errors that are increased with stimulation (though again, cathodal leads to faster latencies, which is unexpected). But again, this pattern is present across the board.

``` r
ggplot(filter(qSig, type == "center", !is.na(significance)), aes(q, fill = significance)) +
  facet_grid(leg ~ direction) +
  geom_bar(position = "stack") +
  stat_bin(binwidth = .1, geom = "text", size = 2, aes(label = ..count..), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#F25F5C", "#4B93B1")) +
  scale_x_continuous("decile", breaks = seq(0.1,0.9,0.1), labels = seq(1,9,1)) +
  scale_y_continuous("number of participants with significant difference", limits = c(0, length(unique(qSig$subject))), breaks = c(0,10,20,length(unique(qSig$subject))))
```

<img src="RT_quantiles_files/figure-markdown_github/Significance shift function center-1.png" width="75%" style="display: block; margin: auto;" />

The strongest effects are again significant in fewer subjects, and on the whole there is never an effect in one direction which is significant in more than half the sample.

Individual shift functions
--------------------------

Let's plot the shift functions for each subject for the left-lateral condition (of main interest as [Kanai et al. (2012)](http://dx.doi.org/10.3389/fpsyt.2012.00045) found effects of anodal tDCS there).

### Baseline

``` r
ggplot(filter(qStats, leg == "baseline", type  == "lateral", direction == "left"), aes(anodal, difference)) +
  facet_wrap(~subject, ncol = 5, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(aes(xintercept = anodal_median), linetype = "dashed", alpha = 0.5) +
  geom_linerange(aes(ymin=ci_lower, ymax=ci_upper), colour = "black", size = 0.5) +
  geom_line(size = 1, colour = "grey50", alpha = 0.5) +
  geom_point(aes(fill = deco), size = 2, colour = "black", shape = 21) +
  scale_fill_gradient(low = "white", high = "grey30", guide = FALSE) +
  xlab("anodal deciles") +
  ylab("anodal - cathodal deciles (ms)") +
  ggtitle("Left lateral saccades during baseline")
```

<img src="RT_quantiles_files/figure-markdown_github/Shift lateral left baseline-1.png" width="75%" style="display: block; margin: auto;" />

### During tDCS

``` r
ggplot(filter(qStats, leg == "tDCS", type  == "lateral", direction == "left"), aes(anodal, difference)) +
  facet_wrap(~subject, ncol = 5, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(aes(xintercept = anodal_median), linetype = "dashed", alpha = 0.5) +
  geom_linerange(aes(ymin=ci_lower, ymax=ci_upper), colour = "black", size = 0.5) +
  geom_line(size = 1, colour = "grey50", alpha = 0.5) +
  geom_point(aes(fill = deco), size = 2, colour = "black", shape = 21) +
  scale_fill_gradient(low = "white", high = "grey30", guide = FALSE) +
  xlab("anodal deciles") +
  ylab("anodal - cathodal deciles (ms)") +
  ggtitle("Left lateral saccades during tDCS")
```

<img src="RT_quantiles_files/figure-markdown_github/Shift lateral left tDCS-1.png" width="75%" style="display: block; margin: auto;" />

### After tDCS (0-15 min)

``` r
ggplot(filter(qStats, leg == "post-1", type  == "lateral", direction == "left"), aes(anodal, difference)) +
  facet_wrap(~subject, ncol = 5, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(aes(xintercept = anodal_median), linetype = "dashed", alpha = 0.5) +
  geom_linerange(aes(ymin=ci_lower, ymax=ci_upper), colour = "black", size = 0.5) +
  geom_line(size = 1, colour = "grey50", alpha = 0.5) +
  geom_point(aes(fill = deco), size = 2, colour = "black", shape = 21) +
  scale_fill_gradient(low = "white", high = "grey30", guide = FALSE) +
  xlab("anodal deciles") +
  ylab("anodal - cathodal deciles (ms)") +
  ggtitle("Left lateral saccades post-tDCS (0-15 min)")
```

<img src="RT_quantiles_files/figure-markdown_github/Shift lateral left post-1-1.png" width="75%" style="display: block; margin: auto;" />

### After tDCS (15-30 min)

``` r
ggplot(filter(qStats, leg == "post-2", type  == "lateral", direction == "left"), aes(anodal, difference)) +
  facet_wrap(~subject, ncol = 5, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(aes(xintercept = anodal_median), linetype = "dashed", alpha = 0.5) +
  geom_linerange(aes(ymin=ci_lower, ymax=ci_upper), colour = "black", size = 0.5) +
  geom_line(size = 1, colour = "grey50", alpha = 0.5) +
  geom_point(aes(fill = deco), size = 2, colour = "black", shape = 21) +
  scale_fill_gradient(low = "white", high = "grey30", guide = FALSE) +
  xlab("anodal deciles") +
  ylab("anodal - cathodal deciles (ms)") +
  ggtitle("Left lateral saccades post-tDCS (15-30 min)")
```

<img src="RT_quantiles_files/figure-markdown_github/Shift lateral left post-2-1.png" width="75%" style="display: block; margin: auto;" />
