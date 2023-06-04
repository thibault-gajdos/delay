rm(list=ls(all=TRUE))  ## efface les données

source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/delay/data/')

## * preparation
## load data
data <- read_csv('data_exp1_johan.csv')
data <- data %>%
    rename(subject_id = `# Subject_id`) %>%
    filter(
        RT>.150,
        MT>0,
        RT<5)

## define variables and  contrasts
data$MT_centered  <- data$MT -mean(data$MT, na.rm = TRUE)
data$PMT_centered  <- data$PMT -mean(data$PMT, na.rm = TRUE)
data$RT_centered  <- data$RT -mean(data$RT, na.rm = TRUE)
data$Accuracy <- as.factor(data$Accuracy)
contrasts(data$Accuracy) <- - contr.sum(2) ## erreur: -1; correct: 1
data$delay <- as.factor(data$delay)
contrasts(data$delay) <- contr.sdif(4)
data$Coherence <- as.factor(data$Coherence)
contrasts(data$Coherence) <- contr.sdif(3)


## * preliminary analysis
p  <- ggplot(data = data, aes(RT)) +
  geom_histogram()  +
  facet_wrap( ~ subject_id) +
  ggtitle('RT')
p


d <- data  %>%
  group_by(Coherence,delay) %>%
  summarise(accuracy = mean(as.numeric(Accuracy)), RT = mean(RT), PMT = mean(PMT), MT = mean(MT))
kable(d, digits = 3)

## |Coherence |delay | accuracy|    RT|   PMT|    MT|
## |:---------|:-----|--------:|-----:|-----:|-----:|
## |0.02      |0     |    1.536| 1.411| 1.230| 0.181|
## |0.02      |3     |    1.560| 0.611| 0.441| 0.170|
## |0.02      |5     |    1.562| 0.574| 0.402| 0.172|
## |0.02      |7     |    1.546| 0.588| 0.410| 0.178|
## |0.11      |0     |    1.762| 1.194| 1.022| 0.171|
## |0.11      |3     |    1.865| 0.615| 0.445| 0.170|
## |0.11      |5     |    1.877| 0.563| 0.390| 0.173|
## |0.11      |7     |    1.861| 0.573| 0.399| 0.174|
## |0.4       |0     |    1.957| 0.806| 0.654| 0.153|
## |0.4       |3     |    1.982| 0.614| 0.445| 0.169|
## |0.4       |5     |    1.982| 0.564| 0.392| 0.173|
## |0.4       |7     |    1.983| 0.568| 0.392| 0.176|

plot.accuracy <- ggplot(data = d, aes(x = Coherence, y = accuracy, color = delay)) +
    geom_point() +
    geom_line(aes(group = delay)) 
plot.accuracy
ggsave('des_accuracy.jpeg', plot.accuracy)

plot.PMT <- ggplot(data = d, aes(x = Coherence, y = PMT, color = delay)) +
    geom_point() +
    geom_line(aes(group = delay)) 
plot.PMT
ggsave('des_PMT.jpeg', plot.PMT)

plot.MT <- ggplot(data = d, aes(x = Coherence, y = MT, color = delay)) +
    geom_point() +
    geom_line(aes(group = delay)) 
plot.MT
ggsave('des_MT.jpeg', plot.MT)

## * Accuracy
l.acc <- lmer_alt(Accuracy  ~ Coherence * delay * RT + (1 + Coherence * delay * RT_centered  || subject_id),
                  family = binomial(link = "logit"),
                  data = data)

## ** frequentist
## Planned model
l.acc <- lmer_alt(accuracy_gabor  ~ size * position * rt_gabor_centered + (1 + size * position * rt_gabor_centered +  || subject_id),
                 family = binomial(link = "logit"),
                 data = data)
summary(l.acc)
summary(rePCA(l.acc))

## remove interactions
l.acc <- lmer_alt(accuracy_gabor  ~ size * position * rt_gabor_centered + (1 + size + position + rt_gabor_centered  || subject_id),
                 family = binomial(link = "logit"),
                 data = data)
summary(l.acc)

## remove: position
l.acc <- lmer_alt(accuracy_gabor  ~ size * position * rt_gabor_centered  + OOZ_time_centered  + (1 + size  + rt_gabor_centered +  OOZ_time_centered || subject_id),
                 family = binomial(link = "logit"),
                 data = data)
summary(l.acc)

save(l.acc, file = 'fit_acc_MG5.rdata')
tab_model(l.acc, file = "acc_MG5.html")

## interaction: size*position*rt
predict <- ggemmeans(l.acc, c('rt_gabor_centered','size','position'))
plot(predict)
plot <- plot(predict) + 
  labs(x = "Centered Response Time (ms)", 
       y = "Accuracy", 
       title = "Response Time x Size x Position interaction on Accuracy") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave('acc_MG5.jpeg', plot)

## ** Bayesian lmer

fit_acc <- brm(accuracy_gabor  ~ size * position * rt_gabor_centered + (1 + size * position * rt_gabor_centered || subject_id),
               family = bernoulli(link = "logit"),
           data = data,
           prior = c(set_prior("normal(0,1)", class = "b")),
           cores = 4, chains = 4,
           control = list(adapt_delta = .98,  max_treedepth = 12),
           iter = 6000,  warmup = 4000, seed = 123,
           save_model = 'acc.stan',
           save_pars = save_pars(all = TRUE)
           )
summary(fit_acc)
save(fit_acc, file ='fit_acc_bayes_MG5.rdata')
tab_model(fit_acc, file = "acc_bayes_MG5.html")

## * Confidence

## ** frequentist

## planned model

l.conf <- clmm(conf_ord ~ (accuracy_gabor * size * position) * (rt_gabor_centered +  OOZ_time_centered)  + (1 * size * position * rt_gabor_centered +  OOZ_time_centered | subject_id),
                  data = data,
                  link = c("probit"))
summary(l.conf)

save(l.conf, file = "fit_conf_MG5.rdata")
tab_model(l.conf, file = "conf_MG5.html")

## plot
predict <- ggemmeans(l.conf, c('size','accuracy_gabor'))
plot(predict)
ggsave('conf_MG5_1.jpeg', plot)

predict <- ggemmeans(l.conf, c('size'))
predict <- as.data.frame(predict) %>%
    rename( confidence = response.level, size = x)

plot <- ggplot(data = predict, aes(x = confidence, y = predicted, colour = size)) +
    geom_point(position = position_dodge(width = .5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = .5, position = "dodge") + 
    labs(y = "Probabilty", 
        title = "Effect of Size on Confidence") +
    theme(plot.title = element_text(hjust = 0.5))
print(plot)
ggsave('conf_MG5_2.jpeg', plot)


## ** bayesian
fit_conf <- brm(conf ~ accuracy_gabor * size * position * rt_gabor_centered + (1 * size * position * rt_gabor_centered | subject_id),
                init_r = 0.05,
        data = data,
    family=cumulative("probit"),
    prior = c(set_prior("normal(0,1)", class = "b")),
    cores = 4, chains = 4,
    control = list(adapt_delta = .95,  max_treedepth = 12),
    iter = 6000,  warmup = 4000, seed = 123,
    save_model = 'conf.stan',
    save_pars = save_pars(all = TRUE)
)
save(fit_conf, file = 'conf.rdata')
tab_model(fit_conf, file = "conf_bayes_MG5.html")


l.time <- lmer_alt(MT_centered  ~  OOZ_time_centered  + (1 | subject_id),
                 data = data)
