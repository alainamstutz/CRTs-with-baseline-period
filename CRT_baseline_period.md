---
title: "CRT_baseline_period"
author: "A.Amstutz"
date: "2025-03-15"
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
### First, start with a simple individual randomized trial

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

In this examples the overall variance is𝜎= 64 (regardless of treatment group) and thus here, the individual-level𝜎= 64 is the only source of variation. The overall effect size 𝛿, which is the difference in average QoL scores across treatment groups, is assumed to be 2.4, a standardized effect size 2.4/8=0.3, i.e., the treatment effect of 2.4 is about 0.3 standard deviations large (Cohen's d), which is a moderate effect size in typical clinical trials.
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

#### Let's generate the data

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

#### Let's estimate the effect size
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

#### Let's confirm the power
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

### Second, now move to a cluster randomized trial, but same context
Parallel cluster randomized trial:
𝑌𝑖𝑗=𝛼+𝛿𝑍𝑗+𝑐𝑗+𝑠𝑖

where 𝑌𝑖𝑗 is a continuous outcome for participant 𝑖 in site 𝑗. 𝑍𝑗 is the treatment indicator for site 𝑗. Again, 𝛿 is the treatment effect. 𝑐𝑗∼𝑁(0,𝜎^2𝑐) is a site level effect, and 𝑠𝑖∼𝑁(0,𝜎^2𝑠) is the participant level effect. The correlation of any two participants in a cluster is 𝜌 (the ICC):

𝜌=𝜎^2𝑐/ (𝜎^2𝑐+𝜎^2𝑠)

This tells us how correlated participants are within the same cluster.
If ρ is close to 0, most of the variability is at the individual level.
If ρ is close to 1, most of the variability is at the site level.

If we have a pre-specified number (𝑛) of participants at each site, we can estimate the sample size required in the CRT applying a design effect 1+(𝑛−1)𝜌 to the sample size of an RCT that has the same overall variance.
We know the overall variance (𝜎^2) is 64 and we assume/know the ICC/p is 0.15. That brings us to (using the ICC formula above)

σ^2 c = 9.6 (site-level variance)
σ^2 s = 54.4 (individual-level variance)
ICC (ρ) = 0.15: 15% of the total variability in the outcome is due to differences between sites, while 85% is due to individual differences within sites.

Since individuals in the same cluster are correlated, we need to adjust the sample size using the design effect.
Design Effect = 1+(n−1)ρ; where, n=30 (number of individuals per site), ρ = 0.15
=> Design effect = 5.35
=> New Total Sample Size = 5.35 × 350 = 1872
Since each site includes 30 participants, the number of required sites is: 62.4 => 64

#### Let's generate the data

```r
simple_crt <- function(nsites, n) {
  # treatment assignment, again, 1:1, 50% prob to be in each group
  defC <- defData(varname = "rx", formula = "1;1", dist = "trtAssign")
  # Define Cluster-Level Data (defC)
  # c (Random Site Effect from Normal Distribution)
  # This represents site-level variation. It follows a normal distribution
  # variance = 9.6 is the site-level variance (see above)
  defC <- defData(defC, varname = "c", formula = "0", variance = 9.6, dist = "normal")  
  # Define Individual-Level Outcome (defS)
  # c: Site-level effect (previously defined).
  # 2.4 * rx: Treatment effect for clusters assigned to intervention, see above
  # s_i \sim N(0, 54.4): Individual-level noise, with variance 54.4.
  defS <- defDataAdd(varname="y", formula="c + 2.4*rx", variance = 54.4, dist="normal")

# site/cluster level data
  # Generates nsites rows (one per cluster)
  # Each row includes: rx (treatment assignment for the site) and c (site-specific random effect).
  dc <- genData(nsites, defC, id = "site")

# individual level data
  # Creates individual-level data by assigning n individuals to each site.
  # Adds outcome (y) based on c, rx, and individual-level variance.
  dd <- genCluster(dc, "site", n, "id")
  dd <- addColumns(defS, dd)
  
  dd[]
}
# generating data for 20 sites, each with 50 individuals.
dd <- simple_crt(20, 50)

# Create violin plot with sites on x-axis and colors for treatment/control
ggplot(dd, aes(x = factor(site), y = y, fill = factor(rx))) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +  # Boxplot for central tendency
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.4) +  # Scatter individual points
  scale_fill_manual(values = c("blue", "red"), labels = c("Control", "Treatment")) +
  labs(x = "Site (Cluster)", y = "Outcome (y)", title = "Outcome Distribution by Site") +
  theme_minimal() +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
```

![](CRT_baseline_period_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

#### Let's estimate the effect size
A mixed effects model is used to estimate the effect size.

```r
dd <- simple_crt(200,100)

fit2 <- lmer(y ~  rx + (1|site), data = dd)
tbl_regression(fit2, tidy_fun = broom.mixed::tidy)  %>% 
  modify_footnote(ci ~ NA, abbreviation = TRUE)
```

```{=html}
<div id="izjetofxgq" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#izjetofxgq table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#izjetofxgq thead, #izjetofxgq tbody, #izjetofxgq tfoot, #izjetofxgq tr, #izjetofxgq td, #izjetofxgq th {
  border-style: none;
}

#izjetofxgq p {
  margin: 0;
  padding: 0;
}

#izjetofxgq .gt_table {
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

#izjetofxgq .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#izjetofxgq .gt_title {
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

#izjetofxgq .gt_subtitle {
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

#izjetofxgq .gt_heading {
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

#izjetofxgq .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#izjetofxgq .gt_col_headings {
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

#izjetofxgq .gt_col_heading {
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

#izjetofxgq .gt_column_spanner_outer {
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

#izjetofxgq .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#izjetofxgq .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#izjetofxgq .gt_column_spanner {
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

#izjetofxgq .gt_spanner_row {
  border-bottom-style: hidden;
}

#izjetofxgq .gt_group_heading {
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

#izjetofxgq .gt_empty_group_heading {
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

#izjetofxgq .gt_from_md > :first-child {
  margin-top: 0;
}

#izjetofxgq .gt_from_md > :last-child {
  margin-bottom: 0;
}

#izjetofxgq .gt_row {
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

#izjetofxgq .gt_stub {
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

#izjetofxgq .gt_stub_row_group {
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

#izjetofxgq .gt_row_group_first td {
  border-top-width: 2px;
}

#izjetofxgq .gt_row_group_first th {
  border-top-width: 2px;
}

#izjetofxgq .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#izjetofxgq .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#izjetofxgq .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#izjetofxgq .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#izjetofxgq .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#izjetofxgq .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#izjetofxgq .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#izjetofxgq .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#izjetofxgq .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#izjetofxgq .gt_footnotes {
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

#izjetofxgq .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#izjetofxgq .gt_sourcenotes {
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

#izjetofxgq .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#izjetofxgq .gt_left {
  text-align: left;
}

#izjetofxgq .gt_center {
  text-align: center;
}

#izjetofxgq .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#izjetofxgq .gt_font_normal {
  font-weight: normal;
}

#izjetofxgq .gt_font_bold {
  font-weight: bold;
}

#izjetofxgq .gt_font_italic {
  font-style: italic;
}

#izjetofxgq .gt_super {
  font-size: 65%;
}

#izjetofxgq .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#izjetofxgq .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#izjetofxgq .gt_indent_1 {
  text-indent: 5px;
}

#izjetofxgq .gt_indent_2 {
  text-indent: 10px;
}

#izjetofxgq .gt_indent_3 {
  text-indent: 15px;
}

#izjetofxgq .gt_indent_4 {
  text-indent: 20px;
}

#izjetofxgq .gt_indent_5 {
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
<td headers="estimate" class="gt_row gt_center">1.2</td>
<td headers="ci" class="gt_row gt_center">0.21, 2.1</td>
<td headers="p.value" class="gt_row gt_center">0.018</td></tr>
    <tr><td headers="label" class="gt_row gt_left">site.sd__(Intercept)</td>
<td headers="estimate" class="gt_row gt_center">3.4</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">Residual.sd__Observation</td>
<td headers="estimate" class="gt_row gt_center">7.4</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
  </tbody>
  
  
</table>
</div>
```

#### Let's confirm the power
Confirm power and go back to the initial assumption of 64 sites with 30 participants per site, which is a massive increase from total sample size 350 (individual RCT) to 1920 in a simple CRT

```r
replicate <- function() {
  dd <- simple_crt(64, 30)
  fit2 <- lmer(y ~  rx + (1|site), data = dd)
  
  coef(summary(fit2))["rx", "Pr(>|t|)"]
}

p_values <- mclapply(1:1000, function(x) replicate(), mc.cores = 4)

mean(unlist(p_values) < 0.05)
```

```
## [1] 0.801
```

```r
# Run the power analysis simulation
p_values <- unlist(p_values)  # Convert list to vector

# Create histogram
ggplot(data.frame(p_values), aes(x = p_values)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Distribution of p-values from 1000 Simulations",
       x = "p-value",
       y = "Frequency") +
  theme_minimal()
```

![](CRT_baseline_period_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

### Third, we move to a CRT with baseline period
The baseline and follow-up measurements can be collected from the same participants (cohort design) or different participants (cross-sectional design), though the impact on the design effect depends on what approach is taken.
The key idea is to measure the same outcome at two different time points to reduce variance and improve statistical power.

𝑌𝑖𝑗𝑘 = 𝛼0 + 𝛼1𝑘 + 𝛿0𝑍𝑗 + 𝛿1𝑘𝑍𝑗 + 𝑐𝑗 + 𝑐𝑝𝑗𝑘 + 𝑠𝑖𝑗 + 𝑠𝑝𝑖𝑗𝑘

where 𝑌𝑖𝑗𝑘 is a continuous outcome measure for individual 𝑖 in site 𝑗 and measurement 𝑘∈{0,1}. 𝑘=0 for baseline measurement, and 𝑘=1 for the follow-up. 
𝑍𝑗 is the treatment status of cluster 𝑗, 𝑍𝑗∈{0,1}. 
𝛼0 is the mean outcome at baseline for participants in the control cluster. Baseline mean outcome in control clusters.
𝛼1 is the change from baseline to follow-up in the control arm. Change over time in the control arm (i.e., how much the outcome changes naturally without intervention).
𝛿0 is the difference at baseline between control and treatment arms (we would expect this to be 0 in a randomized trial)  
𝛿1 is the difference in the change from baseline to follow-up between the two arms. In a randomized trial, since 𝛿0 should be close to 0, 𝛿1 is the treatment effect. Treatment effect, i.e., the difference in change from baseline to follow-up between treatment and control arms.

The model has cluster-specific and individual-specific random effects. 
For both, there can be time-invariant effects and time-varying effects.
Cluster-level effects (between-site variation):
𝑐𝑗∼𝑁(0,𝜎2𝑐) are time invariant site-specific effects. Intrinsic differences between clusters.
𝑐𝑝𝑗𝑘 ∼𝑁(0,𝜎2𝑐𝑝) are the site-specific period (time varying) effects. Changes over time within a site.
Individual-level effects (within-site variation):
𝑠𝑖𝑗∼𝑁(0,𝜎2𝑠) are time invariant individual-level effects, stable individual characteristics.
𝑠𝑝𝑖𝑗𝑘∼𝑁(0,𝜎2𝑠𝑝) are the individual-level period (time varying) effects, measurement error or individual change over time.

Why Does This Help Reduce Sample Size?
1. Since each individual (or site) has two measurements, we can control for their baseline score when estimating the treatment effect.
2. This reduces residual variance, leading to a lower design effect and a smaller required sample size compared to a standard CRT.
3. The benefit depends on the intra-cluster correlation (ICC) and how much baseline values predict follow-up values.

#### Let's generate the data

```r
# Parameters:

# effect: The treatment effect to be used for the formula in the outcome variable.
# nsites: The number of sites (clusters) in the study.
# n: The number of participants at each site.
# s_c, s_cp, s_s, s_sp: These represent the variances of the respective random effects:
#  s_c: Variance of cluster-level (site-level) time-invariant random effects.
#  s_cp: Variance of cluster-level (site-level) time-varying random effects (due to different periods).
#  s_s: Variance of participant-level time-invariant random effects.
#  s_sp: Variance of participant-level time-varying random effects (due to periods).

crt_base <- function(effect, nsites, n, s_c, s_cp, s_s, s_sp) {
  # Variable c = cluster level, with variance s_c. It represents the time-invariant random effects at the cluster level.
  defC <- defData(varname = "c", formula = 0, variance = "..s_c")
  defC <- defData(defC, varname = "rx", formula = "1;1", dist = "trtAssign")
  # Cluster-Level Time-Varying Effects (defCP)
  # This defines time-varying cluster-level random effects c.p with variance s_cp. This effect captures the change over time within clusters, for example, the different responses between baseline and follow-up at each site.
  defCP <- defDataAdd(varname = "c.p", formula = 0, variance = "..s_cp")
  # This defines the individual-level random effects s with variance s_s. This represents individual individual variability that is constant across time.
  defS <- defDataAdd(varname = "s", formula = 0, variance = "..s_s")
  # The treatment effect (effect) for the interaction between the treatment indicator (rx) and the period (period), which gives the treatment effect over time.
  defSP <- defDataAdd(varname = "y",
    formula = "..effect * rx * period + c + c.p + s", 
    variance ="..s_sp")
  
  dc <- genData(nsites, defC, id = "site")
  
  # Add Periods to Cluster Data (dcp)
  # This adds two periods (baseline and follow-up) to each cluster and adds the time-varying random effects (c.p).
  dcp <- addPeriods(dc, 2, "site")
  dcp <- addColumns(defCP, dcp)
  dcp <- dcp[, .(site, period, c.p, timeID)]
  
  # Generate Individual-Level Data (ds)
  ds <- genCluster(dc, "site", n, "id")
  ds <- addColumns(defS, ds)
  
  # Add Periods to Individual-Level Data (dsp)
  # This adds two periods (baseline and follow-up) for each individual, creating the observational ID obsID.
  dsp <- addPeriods(ds, 2)
  setnames(dsp, "timeID", "obsID")
  
  setkey(dsp, site, period)
  setkey(dcp, site, period)
  
  # Merge Cluster-Level and Individual-Level Data (dd)
  dd <- merge(dsp, dcp)
  dd <- addColumns(defSP, dd)
  setkey(dd, site, id, period)
  
  dd[]
}
```

Design effect

In their paper, Teerenstra et al develop a design effect that takes into account the baseline measurement.

The correlation of two participant measurements in the same cluster and same time period is the ICC or 𝜌, and is:

𝜌 = (𝜎2𝑐 + 𝜎2𝑐𝑝) / (𝜎2𝑐 + 𝜎2𝑐𝑝 + 𝜎2𝑠 + 𝜎2𝑠𝑝)

In order to estimate the design effect, we need two more correlations. 
First, the correlation between baseline and follow-up random effects at the cluster level:

𝜌𝑐 = 𝜎2𝑐/ (𝜎2𝑐 + 𝜎2𝑐𝑝)

Second, the correlation between baseline and follow-up random effects at the individual level:

𝜌𝑠 = 𝜎2𝑠/ (𝜎2𝑠+𝜎2𝑠𝑝)

A value 𝑟 is used to estimate the design effect, and is defined as:

𝑟= (𝑛𝜌𝜌𝑐 + (1−𝜌)𝜌𝑠) / (1 + (𝑛−1)𝜌)

If we are able to collect baseline measurements and our focus is on estimating 𝛿1 from the model, the design effect is slightly modified from before:

(1+(𝑛−1)𝜌) * (2(1−𝑟))

-> Basically, adding the second part with rho. r: The correlation factor, which is a weighted combination of the cluster-level and individual-level correlations.
The first part of the formula is similar to the design effect formula used when baseline data isn't available. It adjusts for the clustering of individuals within sites (clusters), accounting for the intra-cluster correlation.
The second part of the formula adjusts for the reduction in variance due to the inclusion of baseline measurements in the analysis. The value of r reflects how strongly the baseline and follow-up measurements correlate at both the cluster and individual levels.
If the correlation between baseline and follow-up measurements is strong (i.e., high ps and high pc), this part of the formula will be small, implying a lower design effect (and thus a lower required sample size).

By collecting baseline measurements, the baseline-to-follow-up correlation reduces the variability in the outcome variable, which allows for more precise estimates of the treatment effect and smaller sample size requirements. Specifically, reducing noise and allowing you to detect smaller treatment effects with the same sample size.
By using the correlation between baseline and follow-up data, you can "shrink" the variability and make the study more efficient. This is crucial if the cost or effort of data collection is high, as it helps you achieve the same statistical power with a smaller sample.


#### Cross-section cohort design
We may not be able to collect two measurements for each participants at a site, but if we can collect measurements on two different cohorts, one at baseline before the intervention is implemented, and one cohort in a second period (either after the intervention has been implemented or not, depending on the randomization assignment of the cluster), we might be able to reduce the number of clusters.

In this case, 𝜎2𝑠= 0 and 𝜌𝑠= 0, so the general model reduces to:

𝑌𝑖𝑗𝑘 = 𝛼0 + 𝛼1𝑘 + 𝛿0𝑍𝑗 + 𝛿1𝑘𝑍𝑗 + 𝑐𝑗 + 𝑐𝑝𝑗𝑘 + 𝑠𝑝𝑖𝑗𝑘

So, we can drop 𝑠𝑖𝑗 because there is no individual-level correlation (independent participants)

The parameters for this simulation are 
𝛿1= 2.4 (treatment effect, see above)
𝜎2𝑐= 6.8 (Cluster-level variance)
𝜎2𝑐𝑝= 2.8 (Cluster-level period variance)
𝜎2𝑠𝑝= 54.4 (Individual-level period variance) 
Total variance:𝜎2𝑐 + 𝜎2𝑐𝑝 + 𝜎2𝑠𝑝 = 6.8 + 2.8 + 54.4 = 64, as used previously.


```r
dd <- crt_base(effect = 2.4, nsites = 20, n = 30, s_c = 6.8, s_cp = 2.8, s_s = 0, s_sp = 54.4)

# Create a new column for period names (baseline and follow-up)
dd[, period_label := ifelse(period == 0, "Baseline", "Follow-up")]

# Create a violin plot, faceted by site
ggplot(dd, aes(x = period_label, y = y, fill = factor(rx))) + 
  geom_violin(trim = FALSE, alpha = 0.1) +  # Make the violin plot transparent
  geom_jitter(aes(color = factor(rx)), width = 0.1, alpha = 0.5) +  # Add individual data points
  geom_boxplot(width = 0.1, color = "black", alpha = 0.3, outlier.shape = NA) +  # Add boxplots
  facet_wrap(~site, ncol = 5) +  # 5 columns for 20 sites
  scale_fill_manual(values = c("blue", "red"), labels = c("Control", "Intervention")) +
  scale_color_manual(values = c("blue", "red")) +
  labs(title = "Violin Plot with Data Points and Boxplots by Site and Treatment Group", 
       x = "Measurement Time", y = "Outcome Value") +
  theme_minimal() + 
  theme(legend.position = "bottom", 
        strip.text = element_text(size = 8), 
        axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust text angles for clarity
```

![](CRT_baseline_period_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

#### Let's estimate the effect size
To estimate the effect size we fit a mixed effect model with cluster-specific effects only (both time invariant and time varying).
Treatment effect under "period*rx" !
Translating the Formula:
Fixed Effects:
Fixed effect for 𝛼1𝑘: This is represented by period in the model, which captures the difference between baseline and follow-up.
Fixed effect for treatment group (rx) 𝛿0𝑍𝑗 : This is represented by rx and captures whether the site received treatment or control.
Interaction of time and treatment (period * rx) 𝛿1𝑘𝑍𝑗 : This interaction term captures how the treatment effect differs between baseline and follow-up.
Random Effects:
(1 | site): 𝑐𝑗 : A random intercept for each site. This accounts for the variability between sites but assumes the variability is constant across time (since it's not tied to the period).
(1 | timeID:site):𝑐𝑝𝑗𝑘: A random intercept for each site in each time period (baseline and follow-up). This accounts for the site-specific changes over time—it captures how variability between sites changes from baseline to follow-up.
s pij is not included in your current model for the cross-sectional design because you are assuming independence between individuals at baseline and follow-up.

```r
dd <- crt_base(effect = 2.4, nsites = 200, n = 100, s_c=6.8, s_cp=2.8, s_s=0, s_sp=54.4)

fit3 <- lmer(y ~ period * rx + (1|timeID:site) + (1 | site), data = dd)
tbl_regression(fit3, tidy_fun = broom.mixed::tidy)  %>% 
  modify_footnote(ci ~ NA, abbreviation = TRUE)
```

```{=html}
<div id="dqyroxmggv" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#dqyroxmggv table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#dqyroxmggv thead, #dqyroxmggv tbody, #dqyroxmggv tfoot, #dqyroxmggv tr, #dqyroxmggv td, #dqyroxmggv th {
  border-style: none;
}

#dqyroxmggv p {
  margin: 0;
  padding: 0;
}

#dqyroxmggv .gt_table {
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

#dqyroxmggv .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#dqyroxmggv .gt_title {
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

#dqyroxmggv .gt_subtitle {
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

#dqyroxmggv .gt_heading {
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

#dqyroxmggv .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#dqyroxmggv .gt_col_headings {
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

#dqyroxmggv .gt_col_heading {
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

#dqyroxmggv .gt_column_spanner_outer {
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

#dqyroxmggv .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#dqyroxmggv .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#dqyroxmggv .gt_column_spanner {
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

#dqyroxmggv .gt_spanner_row {
  border-bottom-style: hidden;
}

#dqyroxmggv .gt_group_heading {
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

#dqyroxmggv .gt_empty_group_heading {
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

#dqyroxmggv .gt_from_md > :first-child {
  margin-top: 0;
}

#dqyroxmggv .gt_from_md > :last-child {
  margin-bottom: 0;
}

#dqyroxmggv .gt_row {
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

#dqyroxmggv .gt_stub {
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

#dqyroxmggv .gt_stub_row_group {
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

#dqyroxmggv .gt_row_group_first td {
  border-top-width: 2px;
}

#dqyroxmggv .gt_row_group_first th {
  border-top-width: 2px;
}

#dqyroxmggv .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#dqyroxmggv .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#dqyroxmggv .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#dqyroxmggv .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#dqyroxmggv .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#dqyroxmggv .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#dqyroxmggv .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#dqyroxmggv .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#dqyroxmggv .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#dqyroxmggv .gt_footnotes {
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

#dqyroxmggv .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#dqyroxmggv .gt_sourcenotes {
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

#dqyroxmggv .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#dqyroxmggv .gt_left {
  text-align: left;
}

#dqyroxmggv .gt_center {
  text-align: center;
}

#dqyroxmggv .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#dqyroxmggv .gt_font_normal {
  font-weight: normal;
}

#dqyroxmggv .gt_font_bold {
  font-weight: bold;
}

#dqyroxmggv .gt_font_italic {
  font-style: italic;
}

#dqyroxmggv .gt_super {
  font-size: 65%;
}

#dqyroxmggv .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#dqyroxmggv .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#dqyroxmggv .gt_indent_1 {
  text-indent: 5px;
}

#dqyroxmggv .gt_indent_2 {
  text-indent: 10px;
}

#dqyroxmggv .gt_indent_3 {
  text-indent: 15px;
}

#dqyroxmggv .gt_indent_4 {
  text-indent: 20px;
}

#dqyroxmggv .gt_indent_5 {
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
    <tr><td headers="label" class="gt_row gt_left">period</td>
<td headers="estimate" class="gt_row gt_center">-0.03</td>
<td headers="ci" class="gt_row gt_center">-0.52, 0.46</td>
<td headers="p.value" class="gt_row gt_center">>0.9</td></tr>
    <tr><td headers="label" class="gt_row gt_left">rx</td>
<td headers="estimate" class="gt_row gt_center">0.17</td>
<td headers="ci" class="gt_row gt_center">-0.78, 1.1</td>
<td headers="p.value" class="gt_row gt_center">0.7</td></tr>
    <tr><td headers="label" class="gt_row gt_left">period * rx</td>
<td headers="estimate" class="gt_row gt_center">2.7</td>
<td headers="ci" class="gt_row gt_center">2.0, 3.4</td>
<td headers="p.value" class="gt_row gt_center"><0.001</td></tr>
    <tr><td headers="label" class="gt_row gt_left">timeID:site.sd__(Intercept)</td>
<td headers="estimate" class="gt_row gt_center">1.6</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">site.sd__(Intercept)</td>
<td headers="estimate" class="gt_row gt_center">2.9</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">Residual.sd__Observation</td>
<td headers="estimate" class="gt_row gt_center">7.4</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
  </tbody>
  
  
</table>
</div>
```

#### Let's update the Design Effect

```r
# Based on the variance assumptions, we can update our design effect:
s_c <- 6.8
s_cp <- 2.8
s_s <- 0
s_sp <- 54.4

rho <- (s_c + s_cp)/(s_c + s_cp + s_s + s_sp)
rho_c <- s_c/(s_c + s_cp)
rho_s <- s_s/(s_s + s_sp)

n <- 30

r <- (n * rho * rho_c + (1-rho) * rho_s) / (1 + (n-1) * rho)

# The design effect for the CRT without any baseline measurement was 5.35. With the two-cohort design, the design effect is reduced slightly:
(des_effect <- (1 + (n - 1) * rho) * 2 * (1 - r))
```

```
## [1] 4.325
```

```r
## [1] 4.3

# and thus leaves us with:
des_effect * 350 / n
```

```
## [1] 50.45833
```

```r
## [1] 50

# The desired number of sites is over 50, so rounding up to the next even number gives us 52
```

#### Let's confirm the power
After calculating the design effect, we run a simulation to confirm the statistical power based on the specified number of sites (52) and participants per site (30). The goal is to check whether the p-value for the treatment effect (period*rx) is statistically significant in at least 80% of simulations.

```r
replicate <- function() {
  dd <- crt_base(2.4, 52, 30, s_c = 6.8, s_cp = 2.8, s_s = 0, s_sp = 54.4)
  fit3 <- lmer(y ~ period * rx + (1|timeID:site) + (1 | site), data = dd)
  
  coef(summary(fit3))["period:rx", "Pr(>|t|)"]
}

p_values <- mclapply(1:1000, function(x) replicate(), mc.cores = 4)

mean(unlist(p_values) < 0.05)
```

```
## [1] 0.805
```

#### Closed cohort design
We can reduce the number of clusters further if instead of measuring one cohort prior to the intervention and another after the intervention, we measure a single cohort twice - once at baseline and once at follow-up. Now we use the full model that decomposes the participant level variance into a time invariant effect (𝑐𝑗) and a time varying effect (𝑐𝑝𝑗𝑘):

𝑌𝑖𝑗𝑘 = 𝛼0 + 𝛼1𝑘 + 𝛿0𝑍𝑗 + 𝛿1𝑘𝑍𝑗 + 𝑐𝑗 + 𝑐𝑝𝑗𝑘 + 𝑠𝑖𝑗 + 𝑠𝑝𝑖𝑗𝑘

#### Let's generate the data
The parameters for this simulation are 
𝛿1 = 2.4 (treatment effect, see above)
𝜎2𝑐 = 6.8 (Cluster-level variance)
𝜎2𝑐𝑝 = 2.8 (Cluster-level period variance)
𝜎𝑠 = 38 (Individual-level variance) 
𝜎2𝑠𝑝 = 16.4 (Individual-level period variance) 
Total variance = 64, as used previously.

```r
dd <- crt_base(effect=2.4, nsites=20, n=30, s_c=6.8, s_cp=2.8, s_s=38, s_sp=16.4)

# Visualizing the data (2 time points: baseline and follow-up)
ggplot(dd, aes(x = period, y = y, group = id, color = as.factor(rx))) +
  geom_line() +
  geom_point(aes(shape = as.factor(rx)), size = 3) +
  facet_wrap(~site, ncol = 5) +
  scale_color_manual(values = c("blue", "red")) +
  labs(title = "Repeated Measurements: Outcome by Site and Time",
       x = "Time Period (0 = Baseline, 1 = Follow-up)",
       y = "Outcome (y)",
       color = "Treatment",
       shape = "Treatment") +
  theme_minimal() +
  theme(legend.position = "bottom")
```

![](CRT_baseline_period_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

#### Let's estimate the effect size
The mixed effect model includes cluster-specific effects only (both time invariant and time varying), as well as subject level effects. 

```r
dd <- crt_base(effect = 2.4, nsites = 200, n = 100, 
  s_c = 6.8, s_cp = 2.8, s_s = 38, s_sp = 16.4)

fit4 <- lmer(y ~ period*rx + (1 | id:site) + (1|timeID:site) + (1 | site), data = dd)
tbl_regression(fit4, tidy_fun = broom.mixed::tidy)  %>% 
  modify_footnote(ci ~ NA, abbreviation = TRUE)
```

```{=html}
<div id="syvravrbas" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#syvravrbas table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#syvravrbas thead, #syvravrbas tbody, #syvravrbas tfoot, #syvravrbas tr, #syvravrbas td, #syvravrbas th {
  border-style: none;
}

#syvravrbas p {
  margin: 0;
  padding: 0;
}

#syvravrbas .gt_table {
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

#syvravrbas .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#syvravrbas .gt_title {
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

#syvravrbas .gt_subtitle {
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

#syvravrbas .gt_heading {
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

#syvravrbas .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#syvravrbas .gt_col_headings {
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

#syvravrbas .gt_col_heading {
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

#syvravrbas .gt_column_spanner_outer {
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

#syvravrbas .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#syvravrbas .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#syvravrbas .gt_column_spanner {
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

#syvravrbas .gt_spanner_row {
  border-bottom-style: hidden;
}

#syvravrbas .gt_group_heading {
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

#syvravrbas .gt_empty_group_heading {
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

#syvravrbas .gt_from_md > :first-child {
  margin-top: 0;
}

#syvravrbas .gt_from_md > :last-child {
  margin-bottom: 0;
}

#syvravrbas .gt_row {
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

#syvravrbas .gt_stub {
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

#syvravrbas .gt_stub_row_group {
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

#syvravrbas .gt_row_group_first td {
  border-top-width: 2px;
}

#syvravrbas .gt_row_group_first th {
  border-top-width: 2px;
}

#syvravrbas .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#syvravrbas .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#syvravrbas .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#syvravrbas .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#syvravrbas .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#syvravrbas .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#syvravrbas .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#syvravrbas .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#syvravrbas .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#syvravrbas .gt_footnotes {
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

#syvravrbas .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#syvravrbas .gt_sourcenotes {
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

#syvravrbas .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#syvravrbas .gt_left {
  text-align: left;
}

#syvravrbas .gt_center {
  text-align: center;
}

#syvravrbas .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#syvravrbas .gt_font_normal {
  font-weight: normal;
}

#syvravrbas .gt_font_bold {
  font-weight: bold;
}

#syvravrbas .gt_font_italic {
  font-style: italic;
}

#syvravrbas .gt_super {
  font-size: 65%;
}

#syvravrbas .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#syvravrbas .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#syvravrbas .gt_indent_1 {
  text-indent: 5px;
}

#syvravrbas .gt_indent_2 {
  text-indent: 10px;
}

#syvravrbas .gt_indent_3 {
  text-indent: 15px;
}

#syvravrbas .gt_indent_4 {
  text-indent: 20px;
}

#syvravrbas .gt_indent_5 {
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
    <tr><td headers="label" class="gt_row gt_left">period</td>
<td headers="estimate" class="gt_row gt_center">-0.21</td>
<td headers="ci" class="gt_row gt_center">-0.73, 0.31</td>
<td headers="p.value" class="gt_row gt_center">0.4</td></tr>
    <tr><td headers="label" class="gt_row gt_left">rx</td>
<td headers="estimate" class="gt_row gt_center">-0.19</td>
<td headers="ci" class="gt_row gt_center">-1.1, 0.73</td>
<td headers="p.value" class="gt_row gt_center">0.7</td></tr>
    <tr><td headers="label" class="gt_row gt_left">period * rx</td>
<td headers="estimate" class="gt_row gt_center">2.4</td>
<td headers="ci" class="gt_row gt_center">1.7, 3.2</td>
<td headers="p.value" class="gt_row gt_center"><0.001</td></tr>
    <tr><td headers="label" class="gt_row gt_left">id:site.sd__(Intercept)</td>
<td headers="estimate" class="gt_row gt_center">6.2</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">timeID:site.sd__(Intercept)</td>
<td headers="estimate" class="gt_row gt_center">1.8</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">site.sd__(Intercept)</td>
<td headers="estimate" class="gt_row gt_center">2.7</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">Residual.sd__Observation</td>
<td headers="estimate" class="gt_row gt_center">4.1</td>
<td headers="ci" class="gt_row gt_center"></td>
<td headers="p.value" class="gt_row gt_center"></td></tr>
  </tbody>
  
  
</table>
</div>
```

#### Let's update the Design Effect

```r
# Based on the variance assumptions, we can update our design effect a second time:
s_c <- 6.8
s_cp <- 2.8
s_s <- 38
s_sp <- 16.4

rho <- (s_c + s_cp)/(s_c + s_cp + s_s + s_sp)
rho_c <- s_c/(s_c + s_cp)
rho_s <- s_s/(s_s + s_sp)

n <- 30

r <- (n * rho * rho_c + (1-rho) * rho_s) / (1 + (n-1) * rho)

# And again, the design effect (and sample size requirement) is reduced:
(des_effect <- (1 + (n - 1) * rho) * 2 * (1 - r))
```

```
## [1] 3.1375
```

```r
## [1] 3.1
des_effect * 350 / n
```

```
## [1] 36.60417
```

```r
## [1] 37

# The desired number of sites is over 36, so I will round up to 38
```

#### Let's confirm the power
After calculating the design effect, we run a simulation to confirm the statistical power based on the specified number of sites (52) and participants per site (30). The goal is to check whether the p-value for the treatment effect (period*rx) is statistically significant in at least 80% of simulations.

```r
replicate <- function() {
  dd <- crt_base(2.4, 38, 30, s_c = 6.8, s_cp = 2.8, s_s = 38, s_sp = 16.4)
  fit4 <-  lmer(y ~ period*rx + (1 | id:site) + (1|timeID:site) + (1 | site), data = dd)
  
  coef(summary(fit4))["period:rx", "Pr(>|t|)"]
}

p_values <- mclapply(1:50, function(x) replicate(), mc.cores = 4) # increase once more computational capacity

mean(unlist(p_values) < 0.05)
```

```
## [1] 0.74
```

```r
## [1] 0.79
```


