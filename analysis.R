library("tidyr")
library("dplyr")
library("ggplot2")
library("ez")

## Load data
dataFile <- file.path("data", "saccLatency_nanmedian.txt")
groupDataWide <- read.table(dataFile, header = TRUE, sep = "\t")

## Restructure into tidy data frame
groupDataLong <- groupDataWide %>%
  select(-X) %>% # drop final empty column named "X"
  gather(ID,saccade.latency,-subject) %>% # make long form
  separate(ID, c("stimulation", "leg", "block", "type", "direction"), sep = "_") %>% # create one column for each factor
  mutate(stimulation = factor(stimulation), leg = factor(leg), block = factor(block), type = factor(type), direction = factor(direction)) %>% # mark columns as factors
  mutate(leg = factor(leg, levels = levels(leg)[c(2,3,1)])) %>% # reorder levels chronologically instead of alphabetically
  arrange(subject) # sort by subject

## Reduced data frame (following Kanai et al)
groupDataReduced <- groupDataLong %>%
  filter(type == "lateral") %>% ## keep only lateral saccades; Kanai did not have central
  group_by(subject,stimulation,direction) %>% # for each condition
  summarise(baseline = mean(saccade.latency[leg == "pre"]), # take average of 3 blocks, make new column
            tDCS = mean(saccade.latency[leg == "tDCS"]),
            post.1 = mean(saccade.latency[leg == "post" & block == c(1,2,3)]),
            post.2 = mean(saccade.latency[leg == "post" & block == c(4,5,6)])) %>%
  gather(leg,saccade.latency,baseline,tDCS,post.1,post.2) %>% # gather new columns to use as factor 
  mutate(leg = factor(leg), leg = factor(leg, levels = levels(leg)[c(1,4,2,3)])) %>% # make leg into factor and reorder levels
  arrange(subject,stimulation,leg) # sort by subject

## Minus baseline
fromBaseline <- groupDataReduced %>%
group_by(subject,stimulation,direction) %>% # for each condition, subtract baseline scores and make new columns
  summarise(tDCS = saccade.latency[leg == "tDCS"] - saccade.latency[leg == "baseline"], 
           post.1 = saccade.latency[leg == "post.1"] - saccade.latency[leg == "baseline"],
           post.2 = saccade.latency[leg == "post.2"] - saccade.latency[leg == "baseline"]) %>%
  gather(leg, saccade.latency, tDCS, post.1, post.2)  %>% # gather new columns to use as factor 
  mutate(leg = factor(leg), leg = factor(leg, levels = levels(leg)[c(3,1,2)])) # make leg into factor and reorder levels

## Plot full
fullPlot <- ggplot(groupDataLong, aes(interaction(block,leg), saccade.latency, color = stimulation, shape = stimulation)) +
  facet_grid(type ~ direction) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1)
  #stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3)
fullPlot

## Kanai plot
kanaiPlot <- ggplot(groupDataReduced, aes(leg, saccade.latency, color = stimulation, shape = stimulation)) +         
  facet_wrap(~direction) +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3)
kanaiPlot

## Individual subjects
kanaiSubsPlot <- ggplot(groupDataReduced, aes(leg, saccade.latency)) +         
  facet_wrap(stimulation ~ direction) +
  geom_line(aes(group=subject,color=subject)) +
  stat_summary(fun.y = mean, aes(group = direction), geom = "line") +
  stat_summary(fun.y = mean, geom = "point")
  #scale_y_continuous(limits = c(120,175))
  #coord_cartesian(ylim = c(120, 175)) 
kanaiSubsPlot

## Kanai plot from baseline
kanaiBasePlot <- ggplot(fromBaseline, aes(leg, saccade.latency, color = stimulation, shape = stimulation)) +
  facet_wrap(~direction) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun.y = mean, geom = "point", size = 3) +
  stat_summary(fun.y = mean, geom = "line", aes(group = stimulation), size = 1) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.3)
kanaiBasePlot

## Kanai plot from baseline for individual subjects
kanaiBaseSubsPlot <- ggplot(fromBaseline, aes(leg, saccade.latency)) +
  facet_wrap(stimulation ~ direction) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(group=subject,color=subject)) +
  stat_summary(fun.y = mean, aes(group = direction), geom = "line") +
  stat_summary(fun.y = mean, geom = "point")
kanaiBaseSubsPlot

## Stats
latencyModel <- ezANOVA(data = data.frame(groupDataReduced), 
                        dv = .(saccade.latency), wid = .(subject), within = .(stimulation, leg, direction), type = 3)

latencyBaseModel <- ezANOVA(data = data.frame(fromBaseline),
                           dv = .(saccade.latency), wid = .(subject), within = .(stimulation, leg, direction), type = 3)
