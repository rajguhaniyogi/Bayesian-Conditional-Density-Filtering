# Bayesian-Conditional-Density-Filtering
Conditional Density Filtering (C-DF) algorithm is devised for efficient online Bayesian inference. C-DF adapts MCMC sampling to the online setting, sampling from approximations to conditional posterior distributions obtained by propagating surrogate conditional sufficient statistics (a function of data and parameter estimates) as new data arrive. These quantities eliminate the need to store or process the entire dataset simultaneously and offer a number of desirable features. Often, these include a reduction in memory requirements and runtime and improved mixing, along with state-of-the-art parameter inference and prediction. These improvements are demonstrated through several illustrative examples. Finally, we show that C-DF samples converge to the target posterior distribution asymptotically as sampling proceeds and more data arrives. The article can be found at https://arxiv.org/abs/1401.3632
The codes are jointly developed by the joint first authors.

lm_code_CDF.R: C-DF code for the ordinary linear regression.

dlm_code_CDF.R: C-DF code for the dynamic linear model.

anova_code_CDF.R: C-DF code the anova model.

probit_code_CDF.R: C-DF code for the high dimensional probit regression model.
