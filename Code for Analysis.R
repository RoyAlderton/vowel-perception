### SPEAKER GENDER AND SALIENCE IN SOCIOLINGUISTIC SPEECH PERCEPTION: GOOSE-FRONTING IN STANDARD SOUTHERN BRITISH ENGLISH
### Author: Roy Alderton

### R CODE FOR STATISTICAL ANALYSIS

# Load packages (install them beforehand if needed)
library(tidyverse)
library(lme4)
library(lmerTest)
library(effects)
library(sjPlot)

setwd("...") # change to file location as appropriate

goose <- read.csv("Main Data File.csv")

## VARIABLES
# participant - participant ID (52 individuals)
# condition - experimental condition (maleFace or femaleFace)
# gender - participant gender (male or female)
# age - participant age in years (ranging from 19 to 30)
# region - participant's UK region of origin (northern, southern, midlands, other)
# ('other' includes Scotland, Northern Ireland and those who grew up in various places)
# distanceLondon - distance in miles between central London and the participant's home town
# (for those who grew up in various places, their main care-giver's home town is used)
# stimulus - FLEECE-GOOSE continuum step (1-11; 1 = most FLEECE-like, 11 = most GOOSE-like)
# response - participant's response in experiment (bee or boo)

# Standardise continuous independent variables
goose_z <- goose %>%
  mutate(age_z = scale(age),
         stimulus_z = scale(stimulus),
         distanceLondon_z = scale(distanceLondon),
         reactionTime_z = scale(reactionTime))

##### Confirmatory model
# Tests the main hypothesis - that condition will affect perception - with stimulus as a control variable
conf.model = glmer(response ~ (1 + stimulus_z | participant) +
                     stimulus_z + condition, 
                   goose_z, family = binomial)
summary(conf.model)
plot(allEffects(conf.model))

# Clearly, condition is not significant (p = 0.52), suggesting no priming effect.
# Now let's explore the data to see if condition is mediated by social characteristics of the listeners

#### Exploratory analysis

# There are two ways we can model listener region of origin - 'region' or 'distanceLondon_z'.
# We don't want both in the model at the same time, so let's do two separate models.

## region model
exp.model.reg = glmer(response ~ (1 + stimulus_z | participant) +
                         stimulus_z + condition + gender + age_z + 
                         region +
                         condition : gender +
                         condition : age_z +
                         condition : region, 
                       goose_z, family = binomial)
summary(exp.model.reg)
plot(allEffects(exp.model.reg))

# There are some possible interactions here between condition*gender and condition*region
# But first, let's check the collinearity of each variable using variance inflation factors (VIFs)
# See https://hlplab.wordpress.com/2011/02/24/diagnosing-collinearity-in-lme4/ for info

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

vif.mer(exp.model.reg)
# Anything above 3 is a bit of a concern (Zuur et al 2010).
# Several of the region levels are a bit high (6 or 7).
# This might be because the cells are not balanced. Let's have a look.
reg <- goose %>%
  group_by(participant) %>%
  count(condition, region) %>%
  group_by(condition, region) %>%
  count(region)
reg
# The balance isn't terrible but it's not as good as it could be.
# The north-south balance is better in the maleFace condition than the femaleFace condition.
# Now let's compare this to a model with distance from London instead.

## distanceLondon_z model
exp.model.dist = glmer(response ~ (1 + stimulus_z | participant) +
                      stimulus_z + condition + gender + age_z + 
                      distanceLondon_z +
                      condition : gender +
                      condition : age_z +
                      condition : distanceLondon_z, 
                    goose_z, family = binomial)
summary(exp.model.dist)
plot(allEffects(exp.model.dist))
# There is a significant interaction between condition and gender.

# Let's check for collinearity.
vif.mer(exp.model.dist)
# Three of the effects are between 3 and 4, suggesting less collinearity than the region model.
# Now let's check the balance across conditions.
lon <- goose %>%
  group_by(condition) %>%
  summarise(mean = mean(distanceLondon, na.rm = T))
lon
# Those in the maleFace condition are from places, on average, that are 30 miles closer to London.
# This reflects the region results to some extent but there are various advantages in choosing this option:
# - Less collinearity than the region model
# - Continuous variables are less likely to cause convergence errors
# - Distance from London makes theoretical sense given that GOOSE-fronting is said to be spreading from there through SE England and beyond.
# Hence, we will go for the distanceLondon_z model.

# exp.model.dist doesn't converge, though, so we need to remove some fixed effects or simplify the random-effects stucture.
# Let's do the former, since the model fit might also improve if we remove non-significant effects.
# So let's remove condition*distanceLondon_z first since it has the highest p-value.

exp.model1 = glmer(response ~ (1 + stimulus_z | participant) +
                    stimulus_z + condition + gender + age_z + 
                    distanceLondon_z +
                    condition : gender +
                    condition : age_z, 
                  goose_z, family = binomial)
summary(exp.model1)
plot(allEffects(exp.model1))

# Now let's do a model comparison
anova(exp.model, exp.model1)
# The model with the interaction is not a better fit than the one without, so let's proceed with exp.model1
# The model still doesn't converge, so let's remove distanceLondon_z, as it's nowhere near significant

exp.model2 = glmer(response ~ (1 + stimulus_z | participant) +
                     stimulus_z + condition + gender + age_z + 
                     condition : gender +
                     condition : age_z, 
                   goose_z, family = binomial)
summary(exp.model2)
plot(allEffects(exp.model2))

# Model comparison
anova(exp.model1, exp.model2)
# As before, including the effect does not improve the model. 
# It also still doesn't converge, so let's remove another insignificant predictor.

exp.model3 = glmer(response ~ (1 + stimulus_z | participant) +
                     stimulus_z + condition + gender + age_z + 
                     condition : gender, 
                   goose_z, family = binomial)
summary(exp.model3)
plot(allEffects(exp.model3))

# Model comparison
anova(exp.model2, exp.model3)

# Again, the interaction doesn't improve the model and there are still convergence issues.
# Let's remove the insignificant effect of age.

exp.model4 = glmer(response ~ (1 + stimulus_z | participant) +
                     stimulus_z + condition + gender + 
                     condition : gender, 
                   goose_z, family = binomial)
summary(exp.model4)
plot(allEffects(exp.model4))

# Model comparison
anova(exp.model3, exp.model4)
# Age as a fixed effect does not improve the model fit.
# However, the model converges this time! Looks like we have a good model.

# The p-value for the interaction is significant but very brittle (0.0499), so let's check it with more decimal places
tab_model(exp.model4, show.adj.icc = TRUE, show.est = TRUE, show.se = TRUE, digits.p = 8)
# Interaction p-value: 0.4994772

# Now let's make a model without the condition*gender interaction just for comparison's sake.
exp.model5 = glmer(response ~ (1 + stimulus_z | participant) +
                     stimulus_z + condition + gender, 
                   goose_z, family = binomial)
summary(exp.model5)
plot(allEffects(exp.model5))

# Model comparison
anova(exp.model4, exp.model5)
# None of the predictors are significant in this model other than stimulus_z.
# The model comparison shows that the interaction does not significantly improve the model fit, but it is very marginal (p = 0.051).
# The interaction from exp.model4 is clearly very brittle but it is worth reporting as an interesting possible effect in the exploratory analysis as long as we use some caveats to interpet it.

# Let's check the collinearity of exp.model4 just to see if things have improved after removing variables.
vif.mer(exp.model4)
# The condition*gender interaction has a VIF slightly above 3 (3.04), but we know that condition and gender are balanced anyway (13 listeners in each cell), so this isn't really a concern.

## INTERPRETATION
# So the male listeners in the female face condition were responding differently to the others.
# Could that be because their regional make-up is different to those in the other cells?
# Let's have a look. Here, region is useful as an interpretative variable as it shows more fine-grained regional differences than distanceLondon.
reg.gen <- goose %>%
  group_by(participant) %>%
  count(condition, gender, region) %>%
  group_by(condition, gender, region) %>%
  count(region)
reg.gen
# It doesn't look like there's very much that's different about the male listeners in the female face condition.
# It has the most northerners (9), but not really much higher than most of the other cells (8/4/7).
# There are some other disadvantages with the balance - e.g.:
# - there are no southern women in the femaleFace condition. 
# - most cells have a majority of northerners apart from women in the male face condition.
# However, I don't think these do a lot to explain the male listeners' behaviour in the female face condition.

### BOUNDARY POINT MODEL, AKA CROSSOVER POINT MODEL

# Sometimes it can be helpful to model each participant's FLEECE-GOOSE boundary point (i.e. where their bee and boo responses cross over at 50%), as in previous studies (eg Drager 2011).
# So let's get each participant's boundary  / crossover point and make a similar model to the exploratory one above.
# Taken from https://stackoverflow.com/questions/34374244/r-how-to-quickly-get-decision-boundary-for-logistic-regression
# We only need to do this once.
# Change the participant ID for each participant (f1-26, m1-26).

gp = glm(response ~ stimulus, data = goose, subset = participant == "f1", family = binomial) # Change participant ID as needed
deviationFromZero <- function(y) abs(predict(gp, data.frame(stimulus = y)))
boundary <- optimise(f = deviationFromZero, interval = range(goose$stimulus))
boundary # the 'minimum' figure is the estimated bee-boo boundary

# I tried to automate the process using the code below but I couldn't get it to work.
# So instead, I had to do it manually for each participant using the code above and save it in Excel...


## Automated version - doesn't work at the moment
# make new data frame
goose.boundaries <- data.frame(participant = unique(goose$participant), boundary = NA, condition = NA, gender = NA)
# run script to get boundary for each unique participant in 'goose$participant' column and add to 'goose.boundaries$boundary' column
for(i in unique(goose$participant))
{
  spk = goose$participant == i  
  # fit model
  model = glm(response ~ stimulus, data = goose, subset = participant == i, family = binomial)
  # calculate deviation from zero
  deviationFromZero <- function (y) abs(predict(model, data.frame(stimulus = y)))
  # get boundary and add to goose.boundaries object
  goose.boundaries$boundary[i] <- optimize(f = deviationFromZero, interval = range(goose$stimulus))
  # add condition and gender info
  goose.boundaries$gender[i] <- goose$gender[i]
  goose.boundaries$condition[i] <- goose$condition[i]
}

# In any case, after it's done once, it can be saved as a CSV and then reloaded (see the other .csv file):
cross <- read.csv("Crossover Point Data File.csv")

# standardise continuous independent variables
cross_z <- cross %>%
  mutate(age_z = scale(age),
         distanceLondon_z = scale(distanceLondon))

# We can model this in the same way as above, but with an ordinary linear model (lm)
# No random effects are needed since there is only one crossover value per participant
# A glm is not needed as the crossover point is a continuous variable, not discrete
# The stimulus variable is not needed as the dependent variable here aggregates the stimuli

cross.model <- lm(goose.boundary ~ condition + gender + age_z + distanceLondon +
                    condition : gender +
                    condition : age_z +
                    condition : distanceLondon,
                  cross_z)
summary(cross.model)
plot(allEffects(cross.model))

# Remove variables according to backwards selection, as in the main exploratory analysis above, resulting in:
cross.model1 <- lm(goose.boundary ~ condition + gender +
                  condition : gender,
                  cross_z)
summary(cross.model1)
plot(allEffects(cross.model1))
# The condition*gender interaction is significant here as well, resulting in a very similar model to the exploratory one we made earlier.
