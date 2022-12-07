# in this script, we demonstrate how to model
# acute and chronic loads separately
# in an overly simplistic approach
# we will not account for non-linearity or other complexities here
library(tidyverse)
library(slider) # for moving EWMA and sums iteratively from one day to the next

# read data
folder_base = "O://Prosjekter//Bache-Mathiesen-004-modelling-recurrent-injury//Data//qsl//"
l_mids_qsl = readRDS(paste0(folder_base, "qsl_imputed.rds"))

d_daily_1 = l_mids_qsl[[1]]

# add day-number per football player
d_daily_1 = d_daily_1 %>% group_by(player_id_num) %>% 
  mutate(day_nr = 1:n()) %>% ungroup()

#--------------------------------------- demonstrating a simplistic approach

# to show how acute and chronic loads can be modelled separately
# in any kind of statistical framework
# we demonstrate in a much simpler approach than the main analysis
# here, acute will be considered 1 week of training, and chronic the 3 concurrent weeks before that week
# essentially the same timeframe as the 7:21 uncoupled ACWR

# Lag the exposure data  7 days, so we're only looking at day 8:21 for chronic load
d_daily_1 = d_daily_1 %>% mutate(exposure_daily_lag_chronic = lag(exposure_daily, 7),
                                 exposure_daily_lag_imputed_chronic = 
                                   ifelse(is.na(exposure_daily_lag_chronic), 
                                          mean(exposure_daily, na.rm = TRUE), exposure_daily_lag_chronic))

# calculate number of players included in analysis
nsub = nrow(d_daily_1 %>% distinct(player_id_num))

# function for calculating exponentially waited moving averages
# using similar syntax as the RA-function
# an exponential smoothing ratio of 2/(n+1)
# same as in williams et al. 2016
ewma = function(x, n_days){
  TTR::EMA(x, n = n_days, wilder = FALSE)
}

# function calculates ewma on a sliding window of lag_max +1 days
# slide is from the slider package
# lag_max + 1 = 21 day window-size
# .before = lag_max, means 20 days prior to current day included in calculation. 
# The current day is included unless otherwise specified, hence, we don't add +1
# step = 1 means the calculation moves 1 day at a time
# complete = TRUE means the calculation waits until the first 21 days of data before first calculation 
# meaning it does not calculate on incomplete windows
slide_ewma = function(x, lag_max){
  l = slide(x, ~ewma(., lag_max+1), .before = lag_max, step = 1, .complete = TRUE) %>% map(last)
  l = compact(l)
  l = unlist(l)
  l
}

# Now, we run the ewma-function per player in the dataset

# EWMA will be calculated on 21 days prior to the acute week
# we will lag the training load data by 7 days, and so, minimum lag is still 0
# lag_max will be 20, representing day nr. 21 (since the counter starts at 0)
lag_max = 20

nested_list = d_daily_1 %>% group_by(player_id_num) %>% nest()
nested_list$data = nested_list$data %>% map(., ~slide_ewma(.$exposure_daily_lag_imputed_chronic, lag_max))
d_chronic_ewma = unnest(nested_list, cols = c(data)) %>% 
  group_by(player_id_num) %>% 
  mutate(day = 1:n()) %>% 
  ungroup() %>% 
  rename(ewma = data)

# we do the same procedure for the weekly sums
slide_sum = function(x){
  l = slide(x, ~sum(.), .before = 6, step = 1, .complete = TRUE)
  l = compact(l)
  l = unlist(l)
  l
}

nested_list = d_daily_1 %>% group_by(player_id_num) %>% nest()
nested_list$data = nested_list$data %>% map(., ~slide_sum(.$exposure_daily))
d_acute_sum = unnest(nested_list, cols = c(data)) %>% 
  group_by(player_id_num) %>% 
  mutate(day = 1:n()) %>% 
  ungroup() %>% 
  rename(weekly_sum = data)

# remove the first 21 rows for to be able to join the EWMA data
# since it does not have the first 21 days
d_daily_mods = d_daily_1 %>% group_by(player_id_num) %>% filter(day_nr >= lag_max+1)
d_daily_mods = d_daily_mods %>% 
  left_join(d_chronic_ewma, by = c("player_id_num", "day_nr" = "day")) %>% 
  left_join(d_acute_sum, by = c("player_id_num", "day_nr" = "day")) %>% ungroup()

# remove weeks where players are not at risk of injury
d_at_risk = d_daily_mods %>% filter(weekly_sum != 0)

# check sample size after these calculations
d_at_risk %>% summarize(sum(injury == 1, na.rm = TRUE))
nrow(d_at_risk)

# run the model
fit_ewma_simple = glm(injury ~ weekly_sum + ewma, 
                      family = "binomial", 
                      data = d_at_risk)

# wald is only for faster computation
# use bootstrap for final model
parameters = parameters::parameters(fit_ewma_simple, ci_method="wald", exponentiate = TRUE)
write_excel_csv(parameters, paste0("simple_approach.csv"), delim = ";", na = "")
#-------------------------------Figures

# make figures
library(directlabels)
library(lmisc) # homemade package
text_size = 16
ostrc_theme =  theme(panel.border = element_blank(), 
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     axis.line = element_line(color = nih_distinct[4]),
                     strip.background = element_blank(),
                     strip.text.x = element_text(size = text_size, 
                                                 family="Trebuchet MS", 
                                                 colour="black", 
                                                 face = "bold", hjust = -0.01),
                     axis.ticks = element_line(color = nih_distinct[4]))


library(ggeffects)
# use quantiles as example levels
hist(d_daily_mods$weekly_sum)
quantile(d_daily_mods$ewma, na.rm = TRUE)

d_preds_ewma_simple = ggpredict(fit_ewma_simple, terms = 
                                  c("weekly_sum [all]", 
                                    "ewma [0, 25, 50, 100]"))

preds_ewma_acute = plot(d_preds_ewma_simple, ci = TRUE) +
  theme_line(text_size) +
  ostrc_theme +
  ylab("Probability of injury") +
  xlab("Minutes in activity per week (acute load)") +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(0, 1000))  +
  scale_color_manual(values = nih_distinct[1:4]) +
  scale_fill_manual(values = nih_distinct[1:4]) +
  ggtitle(NULL)

direct.label.ggplot(preds_ewma_acute, method = list("top.bumpup", cex = 1.2))

hist(d_daily_mods$ewma)
quantile(d_daily_mods$weekly_sum, na.rm = TRUE)
d_preds_ewma_simple_chronic = ggpredict(fit_ewma_simple, terms = 
                                          c("ewma [all]", "weekly_sum [0, 100, 500, 1000]"))

preds_ewma_chronic = plot(d_preds_ewma_simple_chronic, ci = TRUE) +
  scale_color_manual(values = nih_distinct[1:4]) +
  scale_fill_manual(values = nih_distinct[1:4]) +
  theme_line(text_size) +
  ostrc_theme +
  ylab("Probability of injury") +
  scale_x_continuous(breaks = scales::breaks_width(50, 0), limits = c(0, 150))  +
  xlab("Average daily minutes in activity prior to acute week (chronic load)") +
  theme(legend.position = "bottom") +
  ggtitle(NULL) +
  scale_y_continuous(limits = c(0, 0.08), labels = axis_percent)

direct.label.ggplot(preds_ewma_chronic, method = list("top.bumpup", cex = 1.2))

# to save figures
cairo_pdf("plot_simple_acute.pdf", width = 6, height = 6)
direct.label.ggplot(preds_ewma_acute, method = list("top.bumpup", cex = 1.2))
dev.off()

cairo_pdf("plot_simple_chronic.pdf", width = 6, height = 6)
direct.label.ggplot(preds_ewma_chronic, method = list("top.bumpup", cex = 1.2))
dev.off()