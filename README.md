# Learning Granger causality for multivariate time series using state-space models

This work was the part of 2102499 Electrical Enginerring Project, academic year 2019, by Anawat Nartkulpat. In this work we studied and applied permutation test in testing the zero and non-zero Granger causality between variables in time series data.


## Abstract

This project aims to develop a scheme to learn a Granger causality of time series data. The
Granger causality is calculated from state-space parameters estimated using subspace identification
method, then permutation test is applied to classify zero and non-zero causality. We explored
the parameters that affect the performance of permutation test. It was shown that increasing the
number of permutations used in the test improves the performance. Moreover, we found that using
all possible permutations can give slightly better performance than randomly sample permutations
but only when there is no p-value correction applied. This result encourages us to use a MonteCarlo permutation test as the correction methods can improve the performance significantly. The
order of state-space models in estimation also affects the performance of permutation test. We
observed through simulations that underestimation of the model order yields worse performance than
overestimation. In the comparison between permutation test and GMM method, the performance
of permutation test was found to be generally higher when ground truth models have sparse GC
matrices. With large number of data, however, GMM method could give a very close performance
to permutation test. When ground truth models have denser GC matrices, the performance of
permutation test become worse when the length of time series data is longer and can be worse than
GMM method. For the computational cost, permutation test suffers heavily from the excessive use
of subspace identification algorithm.

