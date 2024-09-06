# Title: CT Power Analyses
# Author: Walter Campbell
# Initial Date: 08/21/2024
# Last Edited Date: 08/21/2024
# Purpose: To produce a couple different estimates of power for the 
# CT Drug Sentencing Proposal, 960000-1701-000-01207.

# Clear workspace
rm(list=ls())

# Install and load packages - More than I ended up using, but can't hurt to have all
packages <- c(
  "arrow", "dataReporter", "data.table", "fs", "Hmisc", "here", "haven",
  "janitor", "lubridate", "magrittr", "purrr", "readr", "skimr", "stringr",
  "tidyr", "tidyselect", "labelled", "dplyr", "DescTools", "broom", 
  "janitor", "gmodels", "crosstable", "flextable", "pwr", "ggplot2", "SimEngine", 
  "tidyverse", "pwrss"
)
install_missing_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) install.packages(new_packages)
}
install_missing_packages(packages)
invisible(lapply(packages, library, character.only = TRUE))

## Determining a Reasonable MDEs ##
#Prop1 are those proportions adjusted for what's observed in the CA paper
prop1 <- c(0.200, 0.330)
prop2 <- c(0.174, 0.304)
hs <-data.frame(prop1, prop2)
hs <- hs %>%
  mutate(MDEs = ES.h(p1=prop1, p2=prop2))
hs

## Classic power estimate ##
library(pwr)
#Varying Sample Size to Get MDES
#Set the power analysis parameters to obtain the MDEs for a test of two proportions
samplesize <- c(23654, 21288, 17740, 14192, 11827, 9461)
power <- (0.80)
siglev <- (0.05)
alternative <- ("two.sided")
#Run the power analysis
mdes <- sapply(samplesize, 
                FUN = function(x) {
                  pwr.2p.test(n = x, 
                             power = power, sig.level = siglev, 
                             alternative = alternative)$h})
plot_df <- data.frame(samplesize, mdes)
library(ggplot2)
ggplot(plot_df, aes(x=samplesize,
                    y=mdes))+geom_point()+geom_line()

#Now assess the sample size needed with a logistic regression which is a much better fit
#Assume that other covariates explain a third of the variance
#For a 20% recidivism rate
pwrss.z.logreg(p0 = 0.200, p1 = 0.174, r2.other.x = 0.33,
               power = 0.80, alpha = 0.05, 
               dist = "normal")

#For a 54% recidivism rate
pwrss.z.logreg(p0 = 0.330, p1 = 0.304, r2.other.x = 0.33,
               power = 0.80, alpha = 0.05, 
               dist = "normal")

## Simulation power estimate ##
#NOTE: With little information from CT OPM, it is difficult to conduct a simulation 
#that includes things that relevant covariates, year dummies, and county dummies. 
#This exercise will be much more useful when CT OPM can either send us a sample 
#file or has time to pull detailed information. With that in mind, simulate 
#data that only includes our time and treatment indicators and recidivism. This
#This limits what we glean from these analyses other than confirming what we
#see in more traditional power analyses. Thinking of this as an initial estimate
#with the best information available and starting point for a more detailed function
#that incorporate more covariates and a more complex model at the design planning phase.

#Create a function that will simulate data with different sizes and baseline proportions
my_power_function <- function(baseline_prop, sample_size) {
  #Store results here
  sig_results <- c()
  #Simulate data 1000 times
  for (i in 1:1000) {
    #Create the pre treatment group
    treatment_pre <- data.frame(condition=rep(1, sample_size), time=rep(0, sample_size), recid=rbinom(n=sample_size, size=1, prob=baseline_prop))
    #Create the post treatment group
    treatment_post <- data.frame(condition=rep(1, sample_size), time=rep(1, sample_size), recid=rbinom(n=sample_size, size=1, prob=baseline_prop-0.026))
    #Create the pre comparison group
    comparison_pre <- data.frame(condition=rep(0, sample_size), time=rep(0, sample_size), recid=rbinom(n=sample_size, size=1, prob=baseline_prop))
    #Create the post comparison group
    comparison_post <- data.frame(condition=rep(0, sample_size), time=rep(1, sample_size), recid=rbinom(n=sample_size, size=1, prob=baseline_prop))
    #Combine all conditions and time periods
    sample <- rbind(treatment_pre, treatment_post, comparison_pre, comparison_post)
    #Run a rudimentary model
    mylogit <- glm(recid ~ condition + time + condition:time, data = sample, family = "binomial")
    summary(mylogit)
    #Store whether the impact of the legislation is significant at an 0.05 level
    sig_results[i]  <- tidy(mylogit)$p.value[4] <= .05
  }
  #Store mean of significance
  sig_results %>%
    mean() %>%
    return()
}

#Try different sample sizes for each condition
sample_sizes_to_try <- c(23654, 21288, 17740, 14192, 11827, 9461)

#Store the findings for the low (20%) recidivism rate
power_levels_low <- c()

#Run for 20% recidivism rate
for (i in 1:6) {
    power_levels_low[i] <- my_power_function(0.20, sample_sizes_to_try[i])
}
power_levels_low

# Where do we cross 80%?
power_results <- tibble(sample = sample_sizes_to_try,
                        power = power_levels_low)
power_results

#Plot the results
ggplot(power_results, 
       aes(x = sample, y = power)) +
  geom_line(color = 'red', linewidth = 1.5) + 
  #Add a horizontal line at 80%
  geom_hline(aes(yintercept = .8), linetype = 'dashed') + 
  #Make it look nice
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Sample Size', y = 'Power') +
  ggtitle("20% Recidivism")

#Store the findings for the high (33%) recidivism rate
power_levels_high <- c()

#Run for 33% recidivism rate
for (i in 1:6) {
  power_levels_high[i] <- my_power_function(0.33, sample_sizes_to_try[i])
}
power_levels_high

# Where do we cross 80%?
power_results <- tibble(sample = sample_sizes_to_try,
                        power = power_levels_high)
power_results

#Plot the results
ggplot(power_results, 
       aes(x = sample, y = power)) +
  geom_line(color = 'red', size = 1.5) + 
  #Add a horizontal line at 80%
  geom_hline(aes(yintercept = .8), linetype = 'dashed') + 
  #Make it look nice
  theme_minimal() + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = 'Sample Size', y = 'Power') +
  ggtitle("33% Recidivism")

knitr::stitch('CT Power Analyses.r')
