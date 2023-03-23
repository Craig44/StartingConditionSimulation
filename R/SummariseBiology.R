#'
#'
#' Summarise Biology between life-histories
#'
#'
source("AuxillaryFunctions.R")
library(dplyr)
library(ggplot2)
library(reshape2)
#-------------------------
### short (flat fish) growing biology
#-------------------------
# Based Life-history params for Long-lived Rock fish from Wetzel and Punt
# what about scale parameters like R0?
fast_bio = list(
  ages = 0:50,
  L_inf = 58,
  K = 0.133,
  t0 = 0,
  M = 0.15,
  a = 2.08e-6,
  b = 3.50,
  m_a50 = 6.4,
  m_ato95 = 1.1,
  f_a50 = 7,
  s_ato95 = 2,
  s_a50 = 5,
  f_ato95 = 2,
  sigma = 0.6,
  h = 0.85,
  sigma_r = 0.6,
  R0 = 2383000
)
#-------------------------
### short (flat fish) growing biology
#-------------------------
# Based Life-history params for Long-lived Rock fish from Wetzel and Punt
# what about scale parameters like R0?
medium_bio = list(
  ages = 0:50,
  L_inf = 34,
  K = 0.115,
  t0 = 0,
  M = 0.08,
  a = 2.08e-6,
  b = 3.17,
  m_a50 = 8.7,
  m_ato95 = 3.1,
  f_a50 = 7,
  s_ato95 = 5,
  s_a50 = 3,
  f_ato95 = 2,
  sigma = 0.6,
  h = 0.65,
  sigma_r = 0.6,
  R0 = 2383000
)
#-------------------------
### Slow growing biology
#-------------------------
# Based Life-history params for Long-lived Rock fish from Wetzel and Punt
# what about scale parameters like R0?
slow_bio = list(
  ages = 0:100,
  L_inf = 64,
  K = 0.047,
  t0 = 0,
  M = 0.05,
  a = 9.76e-6,
  b = 3.17,
  m_a50 = 19.6,
  m_ato95 = 6.1,
  f_a50 = 15,
  s_ato95 = 7,
  s_a50 = 10,
  f_ato95 = 7,
  sigma = 0.6,
  h = 0.5,
  sigma_r = 0.6,
  R0 = 5234132
)


###
# Save these objects
# for later reference and R-scripts
###
saveRDS(slow_bio, file = file.path(DIR$data, "Slow_biology.RDS"))
saveRDS(medium_bio, file = file.path(DIR$data, "Medium_biology.RDS"))
saveRDS(fast_bio, file = file.path(DIR$data, "Fast_biology.RDS"))

###
# Plot mean length at age
###
slow_length_at_age = data.frame(age = slow_bio$ages, length_at_age = vonbert(slow_bio$ages, slow_bio$K, L_inf = slow_bio$L_inf, t0 = slow_bio$t0), life_history = "Slow")
medium_length_at_age = data.frame(age = medium_bio$ages, length_at_age = vonbert(medium_bio$ages, medium_bio$K, L_inf = medium_bio$L_inf, t0 = medium_bio$t0), life_history = "Medium")
fast_length_at_age = data.frame(age = fast_bio$ages, length_at_age = vonbert(fast_bio$ages, fast_bio$K, L_inf = fast_bio$L_inf, t0 = fast_bio$t0), life_history = "Fast")
# combine them
length_at_age_df = rbind(slow_length_at_age, medium_length_at_age, fast_length_at_age)

# visualize them
ggplot(length_at_age_df, aes(x = age, y = length_at_age, col = life_history, linetype = life_history)) +
  geom_line(linewidth = 1.1) +
  labs(x = "Age", y = "Length (cm)", col = "Life\nHistory", linetype = "Life\nHistory") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
  )
ggsave(filename = file.path(DIR$fig, "length_at_age.png"), width = 7, height = 6)


###
# Plot Maturity at age
###
slow_maturity = data.frame(age = slow_bio$ages, length = slow_length_at_age$length_at_age, maturity = logis_sel(slow_bio$ages, slow_bio$m_a50, slow_bio$m_ato95), life_history = "Slow")
medium_maturity  = data.frame(age = medium_bio$ages, length = medium_length_at_age$length_at_age, maturity = logis_sel(medium_bio$ages, medium_bio$m_a50, medium_bio$m_ato95), life_history = "Medium")
fast_maturity  = data.frame(age = fast_bio$ages, length = fast_length_at_age$length_at_age, maturity = logis_sel(fast_bio$ages, fast_bio$m_a50, fast_bio$m_ato95), life_history = "Fast")
# combine them
maturity_at_age_df = rbind(slow_maturity, medium_maturity, fast_maturity)

# visualize it
ggplot(maturity_at_age_df, aes(x = age, y = maturity, col = life_history, linetype = life_history)) +
  geom_line(linewidth = 1.1) +
  labs(x = "Age", y = "Proportion Mature", col = "Life\nHistory", linetype = "Life\nHistory") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
  )
ggsave(filename = file.path(DIR$fig, "proportion_mature.png"), width = 7, height = 6)


## compare maturity assumptions with Original Wetzel & Punt 
length_midpoints = seq(from = 5, to = 70, by = 1)

slow_maturity_alt = data.frame(length = length_midpoints, maturity = logistic_sel_alternative(length_midpoints, 38, -0.44), life_history = "Slow")
medium_maturity_alt  = data.frame(length = length_midpoints, maturity = logistic_sel_alternative(length_midpoints, 21, -0.67), life_history = "Medium")
fast_maturity_alt  = data.frame(length = length_midpoints, maturity = logistic_sel_alternative(length_midpoints, 33, -0.75), life_history = "Fast")
maturity_at_age_alt_df = rbind(slow_maturity_alt, medium_maturity_alt, fast_maturity_alt)
maturity_at_age_alt_df$type = "Original"
maturity_at_age_df$type = "Current"

full_maturity = rbind(maturity_at_age_alt_df, maturity_at_age_df %>% select(-age))

# visualize maturity by length
ggplot(full_maturity, aes(x = length, y = maturity, col = type, linetype = type)) +
  geom_line(linewidth = 1.1) +
  labs(x = "Length (cm)", y = "Proportion Mature", col = "", linetype = "") +
  facet_wrap(~life_history)+
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
  )

ggsave(filename = file.path(DIR$fig, "proportion_mature_by_length.png"), width = 11, height = 6)

###
# Plot Survey selectivity at age
###
slow_sel = data.frame(age = slow_bio$ages, length = slow_length_at_age$length_at_age, maturity = logis_sel(slow_bio$ages, slow_bio$m_a50, slow_bio$m_ato95), survey = logis_sel(slow_bio$ages, slow_bio$s_a50, slow_bio$s_ato95), fishery = logis_sel(slow_bio$ages, slow_bio$f_a50, slow_bio$f_ato95), life_history = "Slow")
medium_sel  = data.frame(age = medium_bio$ages, length = medium_length_at_age$length_at_age, maturity = logis_sel(medium_bio$ages, medium_bio$m_a50, medium_bio$m_ato95), survey = logis_sel(medium_bio$ages, medium_bio$s_a50, medium_bio$s_ato95), fishery = logis_sel(medium_bio$ages, medium_bio$f_a50, medium_bio$f_ato95), life_history = "Medium")
fast_sel  = data.frame(age = fast_bio$ages, length = fast_length_at_age$length_at_age, maturity = logis_sel(fast_bio$ages, fast_bio$m_a50, fast_bio$m_ato95), survey = logis_sel(fast_bio$ages, fast_bio$s_a50, fast_bio$s_ato95), fishery = logis_sel(fast_bio$ages, fast_bio$f_a50, fast_bio$f_ato95), life_history = "Fast")
# combine them
selectivity_at_age_df = rbind(slow_sel, medium_sel, fast_sel)

# visualize it
ggplot(selectivity_at_age_df %>% filter(life_history == "Slow")) +
  geom_line(aes(x = age, y = maturity, col = "Maturity", linetype = "Maturity"), linewidth = 1.1) +
  geom_line(aes(x = age, y = fishery, col = "Fishery", linetype = "Fishery"), linewidth = 1.1) +
  geom_line(aes(x = age, y = survey, col = "Survey", linetype = "Survey"), linewidth = 1.1) +
  labs(x = "Age", y = "Ogive", col = "", linetype = "", title = "Slow life-history") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave(filename = file.path(DIR$fig, "Slow_sel.png"), width = 7, height = 6)


ggplot(selectivity_at_age_df %>% filter(life_history == "Medium")) +
  geom_line(aes(x = age, y = maturity, col = "Maturity", linetype = "Maturity"), linewidth = 1.1) +
  geom_line(aes(x = age, y = fishery, col = "Fishery", linetype = "Fishery"), linewidth = 1.1) +
  geom_line(aes(x = age, y = survey, col = "Survey", linetype = "Survey"), linewidth = 1.1) +
  labs(x = "Age", y = "Ogive", col = "", linetype = "", title = "Medium life-history") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave(filename = file.path(DIR$fig, "Medium_sel.png"), width = 7, height = 6)


ggplot(selectivity_at_age_df %>% filter(life_history == "Fast")) +
  geom_line(aes(x = age, y = maturity, col = "Maturity", linetype = "Maturity"), linewidth = 1.1) +
  geom_line(aes(x = age, y = fishery, col = "Fishery", linetype = "Fishery"), linewidth = 1.1) +
  geom_line(aes(x = age, y = survey, col = "Survey", linetype = "Survey"), linewidth = 1.1) +
  labs(x = "Age", y = "Ogive", col = "", linetype = "", title = "Fast life-history") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave(filename = file.path(DIR$fig, "Fast_sel.png"), width = 7, height = 6)








