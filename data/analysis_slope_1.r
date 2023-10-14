rm(list=ls(all=TRUE))  ## efface les données

source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/delay/data/')

## * preparation
## load data
data <- read_csv('slope_data_norm_exp1_onset_tukey_60ms.csv')
data <- data %>%
    mutate(subject =row_number()) %>%
  pivot_longer(
    cols = -subject,  # Replace 'subject_id' with the actual subject identifier column
    names_to = c("delay", "coh"),  # Specify the new column names
    names_sep = "s_coh",
    values_to = "slope"
  )


## define variables and  contrasts
data$delay <- as.factor(data$delay)
contrasts(data$delay) <- contr.sdif(4)
data$coh <- factor(data$coh, levels = c("2", "11", "40"))
contrasts(data$coh) <- contr.sdif(3)



## * preliminary analysis
d <- data  %>%
  group_by(coh,delay) %>%
  summarise(slope = mean(slope))

plot.slope <- ggplot(data = d, aes(x = coh, y = slope, color = delay)) +
    geom_point() +
    geom_line(aes(group = delay))

plot.slope
ggsave('slope_1.jpg', plot.slope)

## * Frequentist
l.slope <- lmer_alt(slope  ~ coh * delay + (1 + coh * delay   || subject),
                  data = data)
summary(l.slope)

tab_model(l.slope, file = "slope_1.html")

predict <- ggemmeans(l.slope, c('delay','coh'))
plot(predict)
plot <- plot(predict) + 
  labs(x = "Delay (s)", 
       y = "slope", 
       title = "Delay:Coherence interaction") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave('l_slope_1.jpeg', plot)
