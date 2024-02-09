# Analysis on domain V reversion

#### Load libraries ####

#### Import Data ####
# GC Subgroup data
data <- read.csv("reversion_data.csv")
dim(data)
mean(data$Days)
median(data$Days)

# POLIS virus data

# POLIS SIA data

#### Clean data ####
# Remove negative SIA1 dates

# Remove negative SIA2 dates

# Create indicator for Dom V reversion (classification 1 or 2)


#### Exclude non-index cVDPV2 detections ####
# Assign emergence group using EPID

# Assign index isolate


#### Cross-check SIA1 field against POLIS ####
# Create new field for SIA in nearby district

#### Calculate expected reversions ####
# Given values
rate_per_day <- 154e-3
lower_bound_rate <- 70e-3
upper_bound_rate <- 297e-3

# Days
days <- c(1, 10, 15, 20, 29)
days <- c(15, 15, 15)
days <- c(1, 1, 1, 29, 29, 29)
days <- c(1, 1, 1, 29*3)
days <- data$Days

# Function to calculate the probability of observing at least one event
calc_prob_at_least_one_event <- function(rate_per_day, days) {
  lambda_value <- rate_per_day * days
  prob_zero_events <- dpois(0, lambda_value)
  prob_at_least_one_event <- 1 - prob_zero_events
  return(prob_at_least_one_event)
}

# Calculating the probability of observing at least one event for point estimate, lower and upper bounds
prob_at_least_one_event <- mean(calc_prob_at_least_one_event(rate_per_day, days))
prob_at_least_one_event_lower_bound <- mean(calc_prob_at_least_one_event(lower_bound_rate, days))
prob_at_least_one_event_upper_bound <- mean(calc_prob_at_least_one_event(upper_bound_rate, days))

# Printing the results
print(paste("Point estimate:", round(prob_at_least_one_event * 100, 2), "%"))
print(paste("Lower bound:", round(prob_at_least_one_event_lower_bound * 100, 2), "%"))
print(paste("Upper bound:", round(prob_at_least_one_event_upper_bound * 100, 2), "%"))

# Hil's simpler method
mean(1-exp(-rate_per_day*days))
mean(1-exp(-lower_bound_rate*days))
mean(1-exp(-upper_bound_rate*days))


#### Hil's experiment into heterogeneity ####
library(tidyverse)

target = 0.73
lambda = 154e-3
t_sia = qexp(target, lambda)
t_sia

# simulated
md = 10
mn = 20
n = 10000

t.s = rlnorm(n, meanlog = log(md), sdlog = sqrt(2*(log(mn) - log(md))))

p = mean(1-exp(-lambda*t.s))
p

df = expand_grid(md = seq(7, 20, 1), mn_mult = c(1.1, 1.5, 2, 4)) |>
  mutate(mn = md*mn_mult)

df = df |> mutate(p = map2_dbl(md, mn,
                               ~mean(
                                 1-exp(-lambda * 
                                         rlnorm(n, meanlog = log(.x), sdlog = sqrt(2*(log(.y) - log(.x))))))))

ggplot(df) + geom_point(aes(x=md, y = p, col = factor(mn_mult)))
