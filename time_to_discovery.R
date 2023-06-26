
library(tidyverse)

ttd_data <- load("C:/Users/coreype/OneDrive - Bill & Melinda Gates Foundation/Documents/GitHub/polio-immunity-mapping/emergence_fit_for_Corey_20230227.rda")


f.l = function(l, m, s) dlnorm(l, m, s)
f.l_composite = function(l, m, s, m_es, s_es)
  dlnorm(l, m, s)*plnorm(l, m_es, s_es, lower.tail = F) +
  dlnorm(l, m_es, s_es)*plnorm(l, m, s, lower.tail = F)
time_to_discovery_w_change_es2_dfonly = function(log_scale, m, s, m_es, s_es, rate = NA){
  tibble(x = seq(0,36,length.out=500)) %>%
    expand_grid(Vaccine = factor(c("Sabin 2", "nOPV2"), levels = c("Sabin 2", "nOPV2"))) %>%
    mutate(val_add = if_else(Vaccine == "Sabin 2", 0, log_scale)) %>%
    mutate(ES = f.l(x, m_es + val_add, s_es),
           AFP = f.l(x, m + val_add, s),
           Composite = f.l_composite(x, m+val_add, s, m_es+val_add, s_es)) %>%
    select(-val_add) %>%
    pivot_longer(cols = -c(x, Vaccine), names_to = "Surv.") %>%
    mutate(name = factor(`Surv.`, levels = c("AFP","ES","Composite"), labels = c("AFP only", "ES", "with ES"))) %>%
    filter(name!="ES")
}

ttd_data = bind_rows(
  time_to_discovery_w_change_es2_dfonly(log(1), m = m, s = s,
                                        m_es = m_es, s_es = s_es, "1") %>%
    filter(Vaccine!="nOPV2"),
  time_to_discovery_w_change_es2_dfonly(log(1.25), m = m, s = s,
                                        m_es = m_es, s_es = s_es, "1") %>%
    filter(Vaccine == "nOPV2") %>%
    mutate(Vaccine = "nOPV2 (1.25x)") %>%
    filter(name=="AFP only"),
  time_to_discovery_w_change_es2_dfonly(log(1.5), m = m, s = s,
                                        m_es = m_es, s_es = s_es, "1") %>%
    filter(Vaccine == "nOPV2") %>%
    mutate(Vaccine = "nOPV2 (1.5x)") %>%
    filter(name=="AFP only")
) %>%
  mutate(Vaccine = factor(Vaccine, levels = c("Sabin 2", "nOPV2 (1.25x)", "nOPV2 (1.5x)")))

ggplot(ttd_data %>% mutate(grp = paste0(Vaccine, name)), aes(x=x, y=value, col=name, lty = Vaccine, group = grp)) + geom_line() +
  xlab('Months')+ ylab(NULL)+
  scale_x_continuous(breaks = seq(0,36,6)) + 
  scale_color_discrete("Surveillance") +
  ggtitle(paste0("Fitted time to discovery w/scenarios")) + theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave(f_name("ttd_scenarios"), height = 4, width = 5)