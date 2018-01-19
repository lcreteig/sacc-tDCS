sacc-tDCS: Questionnaires
================
Leon Reteig

-   [PANAS](#panas)
    -   [Load data](#load-data)
    -   [Post-pre differences per item](#post-pre-differences-per-item)
    -   [Composite positive/negative scores](#composite-positivenegative-scores)
    -   [Statistics](#statistics)
-   [tDCS sensations](#tdcs-sensations)
    -   [Load data](#load-data-1)
    -   [Plot distributions](#plot-distributions)
    -   [Statistics](#statistics-1)

R notebook for analysis of questionnaires in the `sacc-tDCS` dataset.

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
library(ez) # ANOVA
library(knitr) # R markdown output (html, pdf, etc.)
# set default output and figure options
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.width = 6, fig.asp = 0.618, out.width = "70%", fig.align = "center")

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
    ## [1] knitr_1.15.1    ez_4.4-0        dplyr_0.5.0     purrr_0.2.2    
    ## [5] readr_1.1.0     tidyr_0.6.1     tibble_1.3.0    ggplot2_2.2.1  
    ## [9] tidyverse_1.1.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] reshape2_1.4.2     splines_3.4.0      haven_1.0.0       
    ##  [4] lattice_0.20-35    colorspace_1.3-2   htmltools_0.3.6   
    ##  [7] yaml_2.1.14        mgcv_1.8-17        nloptr_1.0.4      
    ## [10] foreign_0.8-67     DBI_0.6-1          modelr_0.1.0      
    ## [13] readxl_1.0.0       plyr_1.8.4         stringr_1.2.0     
    ## [16] MatrixModels_0.4-1 munsell_0.4.3      gtable_0.2.0      
    ## [19] cellranger_1.1.0   rvest_0.3.2        psych_1.7.3.21    
    ## [22] evaluate_0.10      forcats_0.2.0      SparseM_1.77      
    ## [25] quantreg_5.33      pbkrtest_0.4-7     parallel_3.4.0    
    ## [28] broom_0.4.2        Rcpp_0.12.10       scales_0.4.1      
    ## [31] backports_1.0.5    jsonlite_1.4       lme4_1.1-13       
    ## [34] mnormt_1.5-5       hms_0.3            digest_0.6.12     
    ## [37] stringi_1.1.5      grid_3.4.0         rprojroot_1.2     
    ## [40] tools_3.4.0        magrittr_1.5       lazyeval_0.2.0    
    ## [43] car_2.1-4          MASS_7.3-47        Matrix_1.2-9      
    ## [46] xml2_1.1.1         lubridate_1.6.0    assertthat_0.2.0  
    ## [49] minqa_1.2.4        rmarkdown_1.5      httr_1.2.1        
    ## [52] R6_2.2.0           nnet_7.3-12        nlme_3.1-131      
    ## [55] compiler_3.4.0

PANAS
=====

Load data
---------

In each session of the experiment, participants completed the [Positive and Negative Affect Scale (PANAS)](http://booksite.elsevier.com/9780123745170/Chapter%203/Chapter_3_Worksheet_3.1.pdf) twice:

1.  **Pre-measurement**: Before starting setup of the tDCS and eye tracker
2.  **Post-measurement**: After the electrodes were removed following task completion

As most participants were native Dutch speakers, we also used a [Dutch language version](http://www.ekgp.ugent.be/pages/nl/vragenlijsten/PANAS.pdf) of the PANAS, as reported by [Engelen et al. (2006)](http://link.springer.com/article/10.1007/BF03087979).

``` r
# Load the data frame
dataFile <- file.path("data", "PANAS.csv")
panasData <- read_csv2(dataFile, col_types = cols( # read in data; make columns into factors
  session = col_factor(c("first","second")),
  stimulation = col_factor(c("anodal","cathodal")),
  time = col_factor(c("pre","post"))
))
head(panasData) # show data frame
```

    ## # A tibble: 6 x 24
    ##   subject session stimulation   time pos.1.interested neg.2.distressed
    ##     <chr>  <fctr>      <fctr> <fctr>            <int>            <int>
    ## 1     S01   first    cathodal    pre                2                1
    ## 2     S01   first    cathodal   post                2                1
    ## 3     S01  second      anodal    pre                1                1
    ## 4     S01  second      anodal   post                1                1
    ## 5     S02   first      anodal    pre                4                1
    ## 6     S02   first      anodal   post                4                1
    ## # ... with 18 more variables: pos.3.excited <int>, neg.4.upset <int>,
    ## #   pos.5.strong <int>, neg.6.guilty <int>, neg.7.scared <int>,
    ## #   neg.8.hostile <int>, pos.9.enthusiastic <int>, pos.10.proud <int>,
    ## #   neg.11.irritable <int>, pos.12alert <int>, neg.13.ashamed <int>,
    ## #   pos.14.inspired <int>, neg.15.nervous <int>, pos.16.determined <int>,
    ## #   pos.17.attentive <int>, neg.18.jittery <int>, pos.19.active <int>,
    ## #   neg.20.afraid <int>

**Factors**:

-   *subject*: subject ID (`S01`, `S02`, etc)
-   *session*: Whether data are from the `first` or `second` session
-   *stimulation*: Whether data are from the `anodal` or `cathodal` session
-   *time*: Whether data are from the `pre` or `post` measurement

**DATA**:

The PANAS contains 20 items, split equally into positive and negative affect. Each item is rated on a Likert scale:

1.  very slightly or not at all
2.  a little
3.  moderately
4.  quite a bit
5.  extremely

To obtain a positive affect score, we sum items: 1 (interested), 3 (excited), 5 (strong), 9 (enthusiastic), 10 (proud), 12 (alert), 14 (inspired), 16 (determined), 17 (attentive) and 19 (active).

To obtain a negative affect score, we sum items: 2 (distressed), 4 (upset), 6 (guilty), 7 (scared), 8 (hostile), 11 (irritable), 13 (ashamed), 15 (nervous), 18 (jittery) and 20 (afraid).

Post-pre differences per item
-----------------------------

``` r
panasPrePost <- panasData %>%
  group_by(subject,session,stimulation) %>% # for each subject, session, stimulation combination
  select(-time) %>% # don't do subtraction for this column
  summarise_each(funs(diff(.))) # subtract pre and post
```

Let's plot the post-pre differences for each item, split for stimulation session.

``` r
panasPrePost %>%
  gather(item, score, pos.1.interested:neg.20.afraid) %>% # gather all PANAS columns so it can be used as a factor
  ggplot(aes(item, score)) +
      stat_summary(fun.y = mean, geom = "point", position = position_dodge(width = 0.6), aes(color = stimulation)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      coord_cartesian(ylim = c(-1,1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```

<img src="questionnaires_files/figure-markdown_github/Plot post-pre differences per item-1.png" width="70%" style="display: block; margin: auto;" />

Most scores become negative, meaning the `post` score was smaller than `pre`. So they generally feel "less" of everything. This seems particularly the case for positive affect (except for *determined*).

For negative items, people become especially less *jittery*, *nervous*, *irritable*. The biggest decreases in positive scores are in the *attentive*, *active*, and *enthusiastic* items.

The biggest differences between the anodal and cathodal sessions seem to be in the positive items (except for *irritable*).

For some positive items, anodal stimulation results in a more positive score, while cathodal results in a more negative score: (*determined* and *insprired*).

Composite positive/negative scores
----------------------------------

``` r
panasComposite <- panasData %>% 
  mutate(positive = rowSums(select(., contains("pos")))) %>% # sum all the positive scores together
  mutate(negative = rowSums(select(., contains("neg")))) %>% # sum all the negative scores together
  select(subject:time, positive, negative) %>% # drop the original PANAS columns
  gather(affect, score, positive, negative) # gather affect so it can be used as a factor
```

``` r
ggplot(panasComposite, aes(time,score)) +
facet_wrap(~affect) +
  geom_point(position = position_jitter(width = 0.2), aes(color = stimulation)) +
  stat_summary(fun.data = mean_cl_normal, position = position_dodge(width = .5), shape = 21, size = 1, aes(fill = stimulation), color = "black")
```

<img src="questionnaires_files/figure-markdown_github/Plot the composite scores-1.png" width="70%" style="display: block; margin: auto;" />

The negative scores are much lower than the positive overall.

The pre to post differences are larger for positive scores, particularly for the cathodal session.

Statistics
----------

### Positive scores

Repeated measures ANOVA with factors STIMULATION (anodal vs. cathodal) and TIME (pre vs. post)

``` r
# To run the ANOVA
panasCompositeStats <- panasComposite %>%
  group_by(subject) %>% # exclude all subjects with missing cases
  filter(sum(is.na(score)) == 0) %>%
  ungroup() %>%
  mutate(subject = factor(subject)) # "subject" must be factorized
```

``` r
modelPositive <- ezANOVA(data = data.frame(filter(panasCompositeStats, affect == "positive")), # Repeated over subjects; type 3 sums of squares (cf. SPSS)
                         dv = .(score), wid = .(subject), within = .(stimulation, time), type = 3)
kable(modelPositive)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect           |  DFn|  DFd|           F|          p| p&lt;.05 |        ges|
|-----|:-----------------|----:|----:|-----------:|----------:|:---------|----------:|
| 2   | stimulation      |    1|   29|   0.3362518|  0.5664791|          |  0.0016411|
| 3   | time             |    1|   29|  20.9326283|  0.0000825| \*       |  0.0662793|
| 4   | stimulation:time |    1|   29|   4.6731175|  0.0390292| \*       |  0.0070037|

</td>
</tr>
</tbody>
</table>
This confirms our observation that subjects' mood becomes less positive after stimulation, particularly in the cathodal session

### Negative scores

Repeated measures ANOVA with factors STIMULATION (anodal vs. cathodal) and TIME (pre vs. post)

``` r
modelNegative <- ezANOVA(data = data.frame(filter(panasCompositeStats, affect == "negative")), # Repeated over subjects; type 3 sums of squares (cf. SPSS)
                        dv = .(score), wid = .(subject), within = .(stimulation, time), type = 3)
kable(modelNegative)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>
|     | Effect           |  DFn|  DFd|           F|          p| p&lt;.05 |        ges|
|-----|:-----------------|----:|----:|-----------:|----------:|:---------|----------:|
| 2   | stimulation      |    1|   29|   0.0986395|  0.7557167|          |  0.0003892|
| 3   | time             |    1|   29|  22.4236431|  0.0000529| \*       |  0.0677905|
| 4   | stimulation:time |    1|   29|   1.3206320|  0.2598687|          |  0.0021153|

</td>
</tr>
</tbody>
</table>
Subjects' mood also become less negative, but this time there is no interaction with stimulation condition, so this is likely a pure effect of time.

tDCS sensations
===============

After each session, subjects also completed a custom questionnaire probing tDCS sensations.

Load data
---------

``` r
# Load the data frame
dataFile <- file.path("data", "tdcs_sensations.csv")
sensData <- read_csv2(dataFile, col_types = cols(
  session = col_factor(c("first","second")),
  stimulation = col_factor(c("anodal","cathodal"))
))
head(sensData) # show data frame
```

    ## # A tibble: 6 x 20
    ##   subject session stimulation itching tingling burning  pain headache
    ##     <chr>  <fctr>      <fctr>   <int>    <int>   <int> <int>    <int>
    ## 1     S01   first    cathodal       2        4       4     2        0
    ## 2     S01  second      anodal       1        1       1     0        0
    ## 3     S02   first      anodal       1        2       2     0        0
    ## 4     S02  second    cathodal       1        2       2     0        0
    ## 5     S03   first    cathodal       0        0       0     0        0
    ## 6     S03  second      anodal       0        0       0     0        0
    ## # ... with 12 more variables: fatigue <int>, dizziness <int>,
    ## #   nausea <int>, conf.itching <int>, conf.tingling <int>,
    ## #   conf.burning <int>, conf.pain <int>, conf.headache <int>,
    ## #   conf.fatigue <int>, conf.dizziness <int>, conf.nausea <int>,
    ## #   felt.more <chr>

They were asked to which degree the following sensations were present during stimulation: *tingling*, *itching sensation*, *burning sensation*, *pain*, *headache*, *fatigue*, *dizziness* and *nausea*. Each was rated on a scale from 0-4:

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

Finally, they filled in whether they felt one electrode more than the other (`felt.more`); and if so, which and for which sensations.

**Factors**:

-   *subject*: subject ID (`S01`, `S02`, etc)
-   *session*: Whether data are from the `first` or `second` session
-   *stimulation*: Whether data are from the `anodal` or `cathodal` session

``` r
idxComplete <- rowSums(is.na(sensData)) != ncol(sensData) - 3 # rows that do not have all NAs (except the 3 factor columns)
# calculate number of questionnaires completed per stimulation type
nAnodal <- sum(sensData$stimulation == "anodal" & idxComplete) 
nCathodal <- sum(sensData$stimulation == "cathodal" & idxComplete)
```

Plot distributions
------------------

### Sensation intensity

``` r
sensData %>%
  select(everything(), -contains("conf"), -felt.more) %>% # drop the confidence columns
  gather(sensation, rating, itching:nausea) %>% # gather sensations so they can be used as a factor
  ggplot(aes(rating, fill = stimulation)) +
    facet_wrap(~sensation, nrow = 2) +
    geom_histogram(position = "stack", binwidth = .5) +
    stat_bin(binwidth = 1, geom = "text", aes(label = ..count..), position = position_stack(vjust = 0.5)) +
    xlim(0.5,4.5) + ylim(0,30) + # exclude "0" ratings that were not present
    ggtitle(paste("Sensations over", nAnodal, "anodal sessions,", nCathodal, "cathodal sessions"))
```

<img src="questionnaires_files/figure-markdown_github/tDCS sensation distributions-1.png" width="70%" style="display: block; margin: auto;" />

*Nausea* was never experienced; *dizziness* is especially rare and mild, as are *headache* and *pain*.

*Burning*, *itching* and *tingling* are most frequently experienced and to a higher degree.

### Sensation confidence

``` r
sensData %>%
  select(contains("conf"), subject, session, stimulation) %>% # drop the rating columns
  gather(sensation, rating, conf.itching:conf.nausea) %>%
  ggplot(aes(rating, fill = stimulation)) +
    facet_wrap(~sensation, nrow = 2) +
    geom_histogram(position = "stack", binwidth = .5) +
    stat_bin(binwidth = 1, geom = "text", aes(label = ..count..), position = position_stack(vjust = 0.5)) +
    xlim(0.5,4.5) +  ylim(0,30) + # exclude "0" ratings that were not present
    ggtitle(paste("Confidence over", nAnodal, "anodal sessions,", nCathodal, "cathodal sessions"))
```

<img src="questionnaires_files/figure-markdown_github/tDCS confidence distributions-1.png" width="70%" style="display: block; margin: auto;" />

For "local" sensations like *burning*, *tingling* and *dizziness*, subjects have high confidence that these are due to tDCS. For more diffuse sensations, like *fatigue* and *headache*, ratings are very low, so their occurence might just be due to performing the task for an extended period of time.

Note that nausea was never reported, so the two "1" confidence ratings are techincally non-sensical.

### Sensation difference between anode and cathode

``` r
sensData %>%
  select(felt.more, subject, session, stimulation) %>% # keep only the relevant column
  mutate(felt.more = factor(felt.more, levels = c("anode", "equal", "cathode"))) %>% # keep colors consistent
  ggplot(aes(stimulation, fill = felt.more)) +
    geom_bar() +
    stat_count(geom = "text", aes(label = ..count..), position = position_stack(vjust = 0.5), color = "white")
```

<img src="questionnaires_files/figure-markdown_github/Sensations for anode and cathode-1.png" width="70%" style="display: block; margin: auto;" />

There's a pretty even split between which electrode people feel the most, so there don't appear to be worrisome biases.

During anodal stimulation, the anode is felt more than the cathode; during cathodal stimulation the cathode is felt more than the anode. In other words, the electrode over the FEF is always felt more than the forehead electrode. However, this pattern is more pronounced in the anodal session.

About a third of people in the cathodal session did not indicate they felt one more than the other; this is slightly less in the anodal session.

Statistics
----------

We will test for every sensation separately whether there was a difference between the anodal and cathodal sessions. Mann-Whitney-U tests are most appropriate here, as the data are ordinal (Likert) and do not look normally distributed.

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

None of the differences are significant. The p-value for nausea is undefined because it was rated 0 by everyone for both anodal and cathodal.
