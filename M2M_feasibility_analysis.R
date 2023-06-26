# M2M feasibility analysis 


# Shift to # rounds to get to 80% immunity, and 90% immunity, instead of R_0
# Consider changing the Y axis
# More recent references?
# isolines at transition, labeled with #rounds (eg, 4, 9, 12, 25, 50)
# use distrinct color scheme from Arie's figures (to avoid confusion)


#### Libraries ####
library(ggplot2)
library(magrittr)
library(tidyverse)

#### Functions ####
#x is coverage
#y is fraction of missingness that is random
#z is number of rounds
#R0 is basic repro number
#Ve is per-dose vaccine efficacy against type 1

Re <- function(x, y, z, R0, Ve){
  
  frac_rep_miss <- (1-x)*(1-y)
  x_effective <- x/(1-frac_rep_miss)
  
  R0*frac_rep_miss +
    R0*(1-frac_rep_miss)*(1 - x_effective*Ve)^z
}

Re_asymp <- function(Re){  #What is Re at the asymptote (err, last obs)
  min(Re)
}

Re1 <- function(Re, target=1){  #time when Re reaches 1. 
  # if(min(Re) < 1){
    min(which(Re<target), 999)
  # } else{"Inf"}
}

#### Sample Calc ####
Re(x = 0.5, y = 0.68, z = seq(from=1, to=100), R0 = 5, Ve = 0.4)
Re(x = 0.5, y = 0.68, z = seq(from=1, to=100), R0 = 5, Ve = 0.4) %>% Re_asymp
Re(x = 0.5, y = 0.68, z = seq(from=1, to=100), R0 = 5, Ve = 0.4) %>% Re1

#### Construct Data ####
x_inputs <- seq(from = 0.0025, to = 1, by = 0.0025)
y_inputs <- seq(from = 0, to = 1, by = 0.0025)
z_input <- 100
R0_input = 5
Ve = 0.4 # 0.5 for bOPV, 0.4 for tOPV

data <- data.frame(x_in = rep(x_inputs, times = length(y_inputs)),
                  y_in = rep(y_inputs, each = length(x_inputs)),
                  z_in = z_input,
                  R0_in = R0_input,
                  Ve_in = Ve)
head(data)

data <- data %>% rowwise %>% 
  mutate(Re_asymp = Re(x = x_in, y = y_in, z = 1:z_in, R0 = R0_in, Ve = Ve_in) %>% Re_asymp) %>%
  rowwise %>% 
  mutate(Re1 = Re(x = x_in, y = y_in, z = 1:z_in, R0 = R0_in, Ve = Ve_in) %>% 
           Re1 %>% na_if(999))

#### Plot ####
ggplot(data, aes(x = x_in, y = y_in)) +
  theme_bw() + theme(panel.grid = element_blank())+
  geom_tile(aes( fill = Re_asymp < 1)) +
  ylab("Fraction of missingness that is random") +
  xlab("Overall Coverage") +
  scale_fill_discrete(name = "Herd Immunity Possible")

ggplot() +
  theme_bw() + theme(panel.grid = element_blank())+
  geom_contour_filled(data = drop_na(data), aes(x = x_in, y = y_in, z = Re1), 
                      breaks = c(4, 6, 9, 12, 25, 50, 100)) +
  geom_point(aes(x = 0.5, y = 0.68)) +
  geom_point(aes(x = 0.7, y = 0.68)) +
  geom_point(aes(x = 0.5, y = 0.9)) +
  labs(y = "Fraction of missingness that is random", x = "Overall Coverage",
       fill = "Number of Rounds to\nreach 80% Immunity") 

ggplot() +
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0.8, color = "grey") +
  geom_hline(yintercept = 0.9, color = "grey") +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 1, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 0.68, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", linetype = "longdash", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 0.5, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", linetype = "dotted", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 1, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 0.68, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", linetype = "longdash", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 0.5, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", linetype = "dotted", size = 1.5) +
  scale_x_continuous(name = "Number of Rounds", breaks = seq(0,12,2)) +
  scale_y_continuous(name = "Immunity", breaks = seq(0, 1, 0.2))

## REDO TWICE. ONE WITH 100, 75, 50, 25, 0 AND AGAIN WITH 100, 50, 0.
ggplot() +
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0.8, color = "grey") +
  geom_hline(yintercept = 0.9, color = "grey") +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 1, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 0.5, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", linetype = "longdash", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 0, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", linetype = "dotted", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 1, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 0.5, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", linetype = "longdash", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 0, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", linetype = "dotted", size = 1.5) +
  scale_x_continuous(name = "Number of Rounds", breaks = seq(0,12,2)) +
  scale_y_continuous(name = "Immunity", breaks = seq(0, 1, 0.2))

ggplot() +
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_hline(yintercept = 0.8, color = "grey") +
  geom_hline(yintercept = 0.9, color = "grey") +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 1, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 0.75, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", linetype = "longdash", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 0.5, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", linetype = "dotdash", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 0.25, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", linetype = "dotted", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.5, y = 0, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "red", linetype = "solid", size = .5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 1, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 0.75, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", linetype = "longdash", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 0.5, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", linetype = "dotdash", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 0.25, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", linetype = "dotted", size = 1.5) +
  geom_line(aes(x = 0:12, y = 1-Re(x = 0.8, y = 0, z = 0:12, R0 = 5, Ve = 0.4)/5),
            color = "blue", linetype = "solid", size = .5) +
  scale_x_continuous(name = "Number of Rounds", breaks = seq(0,12,2)) +
  scale_y_continuous(name = "Immunity", breaks = seq(0, 1, 0.2))

#### Plot by repeated missing ####
# Problem is rep_miss can appear multiple times with different combinations of coverage and fraction of missingness that is random

data$rep_miss <- (1-data$x_in)*(data$y_in)

ggplot() +
  theme_bw() + theme(panel.grid = element_blank())+
  geom_tile(data = data, aes(x = x_in, y = y_in, fill = Re1)) +
  # geom_contour_filled(data = drop_na(data), aes(x = x_in, y = rep_miss, z = Re1), 
  #                     breaks = c(2, 4, 6, 9, 12, 25, 50, 100)) +
  labs(y = "Fraction of population\nrepeatedly missed", x = "Overall Coverage",
       fill = "Number of Rounds to\nreach 80% Immunity") 

View(data)