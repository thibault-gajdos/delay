rm(list=ls(all=TRUE))  ## efface les données

source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/delay/data/')

## * preparation
## load data
data <- read_csv('data_exp1_johan.csv')
data <- data %>%
    rename(subject_id = `# Subject_id`) %>%
    filter(
        RT>.20,
        MT>0,
        RT<5,
        delay>0,
        response_moment == 1,
        Trial_Type == 5
    )

## define variables and  contrasts
data$MT_centered  <- data$MT -mean(data$MT, na.rm = TRUE)
data$PMT_centered  <- data$PMT -mean(data$PMT, na.rm = TRUE)
data$RT_centered  <- data$RT -mean(data$RT, na.rm = TRUE)
data$Accuracy <- as.factor(data$Accuracy)
contrasts(data$Accuracy) <- - contr.sum(2) ## erreur: -1; correct: 1
data$delay <- as.factor(data$delay)
contrasts(data$delay) <- contr.sdif(3)
data$Coherence <- as.factor(data$Coherence)
contrasts(data$Coherence) <- contr.sdif(3)

## * preliminary analysis
p  <- ggplot(data = data, aes(RT)) +
    geom_histogram()  +
    facet_wrap( ~ subject_id) +
    ggtitle('RT')


p  <- ggplot(data = data, aes(PMT)) +
  geom_histogram()  +
  facet_wrap( ~ subject_id) +
  ggtitle('PMT')

p  <- ggplot(data = data, aes(MT)) +
  geom_histogram()  +
  facet_wrap( ~ subject_id) +
  ggtitle('MT')


d <- data  %>%
  group_by(Coherence,delay) %>%
  summarise(accuracy = mean(as.numeric(Accuracy)), RT = mean(RT), PMT = mean(PMT), MT = mean(MT))


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
## ** Frequentist
l.acc <- lmer_alt(Accuracy  ~ Coherence * delay * RT_centered + (1 + Coherence * delay * RT_centered  || subject_id),
                  family = binomial(link = "logit"),
                  data = data)
summary(l.acc)
## singular
l.acc <- lmer_alt(Accuracy  ~ Coherence * delay * RT_centered + (1 + Coherence * delay + RT_centered  || subject_id),
                  family = binomial(link = "logit"),
                  data = data)
summary(l.acc)

save(l.acc, file = 'reg_accuracy.rdata')
tab_model(l.acc, file = "accuracy.html")

predict <- ggemmeans(l.acc, c('delay','Coherence'))
plot(predict)
plot <- plot(predict) + 
  labs(x = "Delay (s)", 
       y = "Accuracy", 
       title = "Delay:Coherence interaction") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave('acc.jpeg', plot)

## ** Bayesian
fit_acc <- brm(Accuracy  ~ Coherence * delay * RT_centered + (1 + Coherence + delay + RT_centered  || subject_id),
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
save(fit_acc, file ='fit_acc_bayes.rdata')
tab_model(fit_acc, file = "acc_bayes.html")


predict <- ggpredict(fit_acc, c('RT_centered','Coherence','delay'))
plot(predict)
plot <- plot(predict) + 
  labs(x = "Delay (s)", 
       y = "Accuracy", 
       title = "Delay:Coherence:RT interaction") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave('acc_bayes.jpeg', plot)

## * MT
## ** frequantist
l.MT <- lmer_alt(MT  ~ delay * Coherence * Accuracy + (1 + delay + Coherence   || subject_id),
                 data = data)
summary(l.MT)
tab_model(l.MT, file = "MT.html")

ems <- emmeans(l.MT, pairwise ~ Coherence | c(delay,Accuracy))
ems

predict <- ggemmeans(l.MT, c('delay','Coherence','Accuracy')) %>%
    rename(Coherence = group, delay = x, Accuracy = facet)
plot.MT <- ggplot(data = predict, aes(x = delay, y = predicted, colour = Coherence)) +
    geom_point(position = position_dodge(width = .5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = .5, position = "dodge") + 
     labs(y = "predicted MT", 
         title = "MT") +
         theme(plot.title = element_text(hjust = 0.5))+
    facet_wrap(~ Accuracy)
plot.MT
ggsave('MT.jpeg', plot.MT)

## ** Bayesian
fit_MT <- brm(MT  ~ delay * Coherence * Accuracy + (1 + delay + Coherence    || subject_id),
           data = data,
           prior = c(set_prior("normal(0,1)", class = "b")),
           cores = 4, chains = 4,
           control = list(adapt_delta = .9,  max_treedepth = 12),
           iter = 6000,  warmup = 4000, seed = 123,
           save_model = 'MT.stan',
           save_pars = save_pars(all = TRUE)
           )
summary(fit_MT)
tab_model(fit_MT, file = "MT_bayes.html")
save(fit_MT, file ='fit_MT_bayes.rdata')

load('fit_MT_bayes.rdata')
ems <- emmeans(fit_MT, pairwise ~ Coherence | c(delay,Accuracy))
ems

predict <- ggpredict(fit_MT, c('delay','Coherence','Accuracy')) %>%
    rename(Coherence = group, delay = x, Accuracy = facet)
plot.MT <- ggplot(data = predict, aes(x = delay, y = predicted, colour = Coherence)) +
    geom_point(position = position_dodge(width = .5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = .5, position = "dodge") + 
     labs(y = "predicted MT", 
         title = "MT") +
         theme(plot.title = element_text(hjust = 0.5))+
    facet_wrap(~ Accuracy)
plot.MT
ggsave('MT_bayes.jpeg', plot.MT)


## * PMT
## ** frequentist
l.PMT <- lmer_alt(PMT  ~ delay * Coherence * Accuracy + (1 + delay + Coherence   || subject_id),
                 data = data)
summary(l.PMT)
tab_model(l.PMT, file = "PMT.html")

ems <- emmeans(l.PMT, pairwise ~ Coherence | c(delay,Accuracy))
ems


predict <- ggemmeans(l.PMT, c('delay','Coherence','Accuracy')) %>%
    rename(Coherence = group, delay = x, Accuracy = facet)
plot.PMT <- ggplot(data = predict, aes(x = delay, y = predicted, colour = Coherence)) +
    geom_point(position = position_dodge(width = .5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = .5, position = "dodge") + 
     labs(y = "predicted PMT", 
         title = "PMT") +
         theme(plot.title = element_text(hjust = 0.5))+
    facet_wrap(~ Accuracy)
plot.PMT
ggsave('PMT.jpeg', plot.PMT)

## ** Bayesian
fit_PMT <- brm(PMT  ~ delay * Coherence * Accuracy + (1 + delay + Coherence    || subject_id),
           data = data,
           prior = c(set_prior("normal(0,1)", class = "b")),
           cores = 4, chains = 4,
           control = list(adapt_delta = .9,  max_treedepth = 12),
           iter = 6000,  warmup = 4000, seed = 123,
           save_model = 'PMT.stan',
           save_pars = save_pars(all = TRUE)
           )
summary(fit_PMT)
save(fit_PMT, file ='fit_PMT_bayes.rdata')
tab_model(fit_PMT, file = "PMT_bayes.html")

load('fit_PMT_bayes.rdata')

ems <- emmeans(fit_PMT, pairwise ~ Coherence | c(delay,Accuracy))
ems

predict <- ggpredict(fit_PMT, c('delay','Coherence','Accuracy')) %>%
    rename(Coherence = group, delay = x, Accuracy = facet)
plot.PMT <- ggplot(data = predict, aes(x = delay, y = predicted, colour = Coherence)) +
    geom_point(position = position_dodge(width = .5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = .5, position = "dodge") + 
     labs(y = "predicted PMT", 
         title = "PMT") +
         theme(plot.title = element_text(hjust = 0.5))+
    facet_wrap(~ Accuracy)
plot.PMT
ggsave('PMT_bayes.jpeg', plot.PMT)
