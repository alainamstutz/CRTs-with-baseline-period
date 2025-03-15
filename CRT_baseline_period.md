---
title: "CRT_baseline_period"
author: "A.Amstutz"
date: "2023-10-18"
output:
  html_document:
    keep_md: yes
    toc: yes
    toc_float: yes
    code_folding: hide
  pdf_document:
    toc: yes
---

## Parallel CRT with baseline period
1. "A common enhancement of a simple parallel CRT is to add an assessment of participants’ outcomes in a baseline period (before randomisation). Even if different participants are assessed at baseline and follow-up [i.e. cross-sectional sampling], the fact that they are sampled from the same cluster allows some control for cluster differences." -> https://www.bmj.com/content/360/bmj.k1121
2. This is illustratively shown in the sample size calculator: https://clusterrcts.shinyapps.io/rshinyapp/ (switch between "Parallel" and "Parallel with baseline measure") -> can yield a substantial increase in power! See last chapter below.

Let's explore this further
* The rationale: The more variability there is in the outcome across clusters, the more difficult to identify the effect.
* On the example of an individual RCT: "if we are interested in measuring the impact of an intervention on the quality of life (QOL) across a diverse range of patients, the measurement (which typically ranges from 0 to 1) might vary considerably from person to person, regardless of the intervention. If the intervention has a real but moderate effect of, say, 0.1 points, it could easily get lost if the standard deviation is considerably larger, say 0.25."
* If we collect baseline QOL scores and can “control” for those measurements in some way (by conducting a repeated measures analysis, using ANCOVA, or assessing the difference itself as an outcome), we might be able to reduce the variability and have better chance to pick up the effect. 

Literature:
* https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.5352 
* Hooper Richard et al. Sample size calculation for stepped wedge and other longitudinal cluster randomised trials. Statistics in Medicine. 2016: https://onlinelibrary.wiley.com/doi/10.1002/sim.7028 and 
* Leyrat Clémence et al. Practical considerations for sample size calculation for cluster randomized trials. Journal of Epidemiology and Population Health. 2024: https://www.sciencedirect.com/science/article/pii/S2950433324000090
* Based on: https://www.bmj.com/content/360/bmj.k1121.long

## There are various ways to do it:
1. Analysis of covariance (ANCOVA): Aggregate outcomes at baseline, and adjusts each individual participant at follow-up for the baseline cluster mean
2. Constrained baseline analysis: Treat outcomes collected at baseline and follow-up as longitudinal, and to use a repeated measures analysis to estimate the effect of the intervention being switched on in one of the randomised groups on the second of these occasions, see design matrix in https://clusterrcts.shinyapps.io/rshinyapp/. Unlike a difference of differences analysis, it assumes that there is no systematic difference between the groups at baseline.

## Simulation
### Let's start with a simple individual randomized trial

```r
RNGkind("L'Ecuyer-CMRG")
set.seed(19287)
library(simstudy)
library(ggplot2)
library(lmerTest)
library(parallel)
library(data.table)
library(pwr)
library(gtsummary)
library(paletteer)
library(magrittr)
```

In this examples the overall variance is𝜎=64 (regardless of treatment group) and thus = individual-level𝜎=64 is the only source of variation. The overall effect size 𝛿, which is the difference in average QoL scores across treatment groups, is assumed to be 2.4, a standardized effect size 2.4/8=0.3, i.e., the treatment effect of 2.4 is about 0.3 standard deviations large (Cohen's d), which is a moderate effect size in typical clinical trials.
Let's calculate the sample size for a two-sample t-test (2 groups, independent, only 1 measurement at the end) based on the effect size (Cohen's d) and a power of 80%.

```r
pwr.t.test(d = 0.3, power = 0.80)
```

```
## 
##      Two-sample t test power calculation 
## 
##               n = 175.3847
##               d = 0.3
##       sig.level = 0.05
##           power = 0.8
##     alternative = two.sided
## 
## NOTE: n is number in *each* group
```
We will need 350 participants (175 in each arm) to achieve a power of 80%.

### Let's generate the data

```r
simple_rct <- function(N) {
  
  # data definition for outcome
  
  defS <- defData(varname = "rx", formula = "1;1", dist = "trtAssign") # treatment variable, 1:1 allocation, 50% prob to be in a or b
  defS <- defData(defS, varname = "y", formula = "2.4*rx", variance = 64, dist = "normal") # outcome, continuous, 
  # If rx = 1 (treatment group), the outcome will be 2.4*1 = 2.4.
  # If rx = 0 (control group), the outcome will be 2.4*0 = 0.
  # variance: spread or variability in the scores within each group.
  # dist = "normal" indicates that the outcome variable (y) will follow a normal distribution.
  # The outcome variable has a mean of 2.4 for the treatment group and 0 for the control group, with a variance of 64.
  dd <- genData(N, defS)
  
  dd[]
}

dd <- simple_rct(350) # simulate the trial with 350 participants

ggplot(dd, aes(x = factor(rx), y = y, fill = factor(rx))) + 
  # Violin plot
  geom_violin(trim = FALSE) + 
  # Add the individual data points
  geom_jitter(width = 0.15, size = 1, alpha = 0.6, color = "black") +
  # Add a point estimate (median)
  stat_summary(fun = "median", geom = "point", shape = 18, size = 3, color = "red") +
  # Boxplot to show interquartile range (optional, can be removed if not needed)
  geom_boxplot(width = 0.1, alpha = 0.5, outlier.shape = NA, color = "black") +
  labs(x = "Treatment Group (rx)", y = "Outcome (y)") +
  scale_fill_manual(values = c("skyblue", "salmon")) + 
  theme_minimal() +
  theme(legend.title = element_blank()) +
  ggtitle("Violin Plots with Point Estimate and Individual Data Points")
```

![](CRT_baseline_period_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

### Let's estimate the effect size
A simple linear regression model will do

```r
fit1 <- lm(y ~ rx, data = dd)
tbl_regression(fit1) %>% 
  modify_footnote(ci ~ NA, abbreviation = TRUE)
```

```{=html}
<div id="gpdyzpkjsn" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#gpdyzpkjsn table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#gpdyzpkjsn thead, #gpdyzpkjsn tbody, #gpdyzpkjsn tfoot, #gpdyzpkjsn tr, #gpdyzpkjsn td, #gpdyzpkjsn th {
  border-style: none;
}

#gpdyzpkjsn p {
  margin: 0;
  padding: 0;
}

#gpdyzpkjsn .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#gpdyzpkjsn .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#gpdyzpkjsn .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#gpdyzpkjsn .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#gpdyzpkjsn .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#gpdyzpkjsn .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#gpdyzpkjsn .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#gpdyzpkjsn .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#gpdyzpkjsn .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#gpdyzpkjsn .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#gpdyzpkjsn .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#gpdyzpkjsn .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#gpdyzpkjsn .gt_spanner_row {
  border-bottom-style: hidden;
}

#gpdyzpkjsn .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#gpdyzpkjsn .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#gpdyzpkjsn .gt_from_md > :first-child {
  margin-top: 0;
}

#gpdyzpkjsn .gt_from_md > :last-child {
  margin-bottom: 0;
}

#gpdyzpkjsn .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#gpdyzpkjsn .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#gpdyzpkjsn .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#gpdyzpkjsn .gt_row_group_first td {
  border-top-width: 2px;
}

#gpdyzpkjsn .gt_row_group_first th {
  border-top-width: 2px;
}

#gpdyzpkjsn .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#gpdyzpkjsn .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#gpdyzpkjsn .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#gpdyzpkjsn .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#gpdyzpkjsn .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#gpdyzpkjsn .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#gpdyzpkjsn .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#gpdyzpkjsn .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#gpdyzpkjsn .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#gpdyzpkjsn .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#gpdyzpkjsn .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#gpdyzpkjsn .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#gpdyzpkjsn .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#gpdyzpkjsn .gt_left {
  text-align: left;
}

#gpdyzpkjsn .gt_center {
  text-align: center;
}

#gpdyzpkjsn .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#gpdyzpkjsn .gt_font_normal {
  font-weight: normal;
}

#gpdyzpkjsn .gt_font_bold {
  font-weight: bold;
}

#gpdyzpkjsn .gt_font_italic {
  font-style: italic;
}

#gpdyzpkjsn .gt_super {
  font-size: 65%;
}

#gpdyzpkjsn .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#gpdyzpkjsn .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#gpdyzpkjsn .gt_indent_1 {
  text-indent: 5px;
}

#gpdyzpkjsn .gt_indent_2 {
  text-indent: 10px;
}

#gpdyzpkjsn .gt_indent_3 {
  text-indent: 15px;
}

#gpdyzpkjsn .gt_indent_4 {
  text-indent: 20px;
}

#gpdyzpkjsn .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;Characteristic&lt;/strong&gt;"><strong>Characteristic</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;Beta&lt;/strong&gt;"><strong>Beta</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;95% CI&lt;/strong&gt;"><strong>95% CI</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;p-value&lt;/strong&gt;"><strong>p-value</strong></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">rx</td>
<td headers="estimate" class="gt_row gt_center">2.7</td>
<td headers="ci" class="gt_row gt_center">1.0, 4.4</td>
<td headers="p.value" class="gt_row gt_center">0.002</td></tr>
  </tbody>
  
  
</table>
</div>
```

### Let's confirm the power
We can confirm the power by repeatedly generating data sets and fitting models, recording the p-values for each replication.

```r
replicate <- function() {
  dd <- simple_rct(350)
  fit1 <- lm(y ~  rx, data = dd)
  coef(summary(fit1))["rx", "Pr(>|t|)"]
}

p_values <- mclapply(1:1000, function(x) replicate(), mc.cores = 4) # mclapply() is a parallelized version of the lapply() function that applies the replicate() function 1000 times. mc.cores = 4 means the code will use 4 cores for parallel computation, which speeds up the process by running multiple replications simultaneously.

# Estimated power based on 1000 replications.
# Convert the list of p-values into a vector and calculates the proportion of replications where the p-value is less than 0.05
mean(unlist(p_values) < 0.05)
```

```
## [1] 0.789
```



