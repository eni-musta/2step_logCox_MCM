# 2step_logCox_MCM


The logistic/Cox mixture cure model is commonly used to handle survival data in the presence of a cure fraction. It assumes that the population consists of two subpopulations: the cured and the noncured ones. Logistic regression is used to model the incidence and a Cox proportional hazards model is used for the latency. Maximum likelihood estimators are computed via the Expectation-Maximization algorithm and are implemented in the package smcure (Cai et al. (2012)). 

We focus on the cure probabilities and propose a two-step procedure to improve upon an initially available estimator, that can for example be the smcure estimator, for small and moderate sample sizes. The initial estimator is used to construct a one-dimensional covariate, conditional on which we compute a nonparametric estimator of the cure probabilities. Afterwards, the nonparametric estimator is projected on the desired parametric class (for example logistic). This allows direct estimation of the parametric incidence component despite the latent cure status. Once the incidence component is estimated, one can also fit a semiparametric model to the latency component as proposed in Musta, Patilea and Van Keilegom (2022a). For more details about the method, the asymptotic properties of the estimators and the finite sample behaviour see Musta, Patilea and Van Keilegom (2022b).

Here we demonstrate the method for the melanoma data (ECOG phase III clinical trial e1684) from the smcure package (Cai et al. (2012)). 
The purpose of this study was to evaluate the effect of treatment (high dose interferon alpha-2b regimen) as the postoperative adjuvant therapy. The event time is the time from initial treatment to recurrence of melanoma and three covariates have been considered: age (continuous variable centered to the mean), gender (0=male and 1=female) and treatment (0=control and 1=treatment). The data consists of 284 observations (after deleting missing data) out of which 196 had recurrence of the melanoma cancer (around 30% censoring). The Kaplan Meier estimator of the survival function can be found in 'KM_estimator.png'.

One can run the example via '2step_logCox-MCM.R'. This produces estimators for the parameters of the incidence, the regression coefficients and the baseline survival of the Cox component (latency). The bootstrap procedure is used to estimate the variance of the parameter estimates and the resulting p-values. The R file uses the packages 'smcure' for the data and 'np' for bandwidth selection. Additional functions are defined in the file 'functions_2step_logCox-MCM.R'.


References

Cai, C., Y. Zou, Y. Peng, and J. Zhang (2012). smcure: An r-package for estimating semiparametric
mixture cure models. Comput. Meth. Prog. Bio. 108 (3), 1255-1260.

Musta, E., Patilea, V., & Van Keilegom, I. (2022a). A presmoothing approach for estimation in the semiparametric Cox mixture cure model. Bernoulli, 28(4), 2689-2715.

Musta, E., Patilea, V., & Van Keilegom, I. (2022b). A 2-step estimation procedure for semiparametric mixture cure models. arXiv preprint arXiv:2207.08237.
