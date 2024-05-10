## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
### Install required packages
# install.packages(c("devtools","bsseq","glmnet", "mice"))


### install the madrid package
# github_install("andyamadrid/madrid")



## ----setup--------------------------------------------------------------------
# load in the madrid pacakge
suppressPackageStartupMessages(library(madrid))



## -----------------------------------------------------------------------------
# get toy example
data(bs)



## -----------------------------------------------------------------------------
# summary of bsseq object
bs

# summary of phenotypic data associated with bsseq object
pData(bs)



## -----------------------------------------------------------------------------
# guess sample sex from bsseq object
guessed_sex <- guessSex(bs)

# look at the output
guessed_sex



## -----------------------------------------------------------------------------
# predict age from bsseq object
predicted_age <- predictAge(bs)

# look at the output
predicted_age



## -----------------------------------------------------------------------------
# assess the correlation between actual and estimated
cor(predicted_age$predictedAge, pData(bs)$Age)



## -----------------------------------------------------------------------------
# test for differences (accelerated/decelerated) in actual and estimated age
predicted_age$difference <- predicted_age$predictedAge - pData(bs)$Age

t.test(predicted_age$difference ~ factor(pData(bs)$Diagnosis))



## -----------------------------------------------------------------------------
# test for differences in the resduals between actual and estimated age
predicted_age$residuals <- resid(lm(predicted_age$predictedAge ~ pData(bs)$Age))
t.test(predicted_age$residuals ~ factor(pData(bs)$Diagnosis))



