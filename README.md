# C++ implementation of the Weighted Histogram Analysis Method (WHAM)

## Validation of Code 
### Free Energy 
Validation of the Unbinned WHAM (UWHAM) using different optimization methods, specifically in this case for adaptive and LBFGS. Using the method of bootstrapping, error was also obtained for estimated free energy.
![wham](/test/validate.png)
### KL Divergence
Their KL divergences is also evaluated and shown below for the 8 biased simulations performed. (1 is unbiased) 
![kl](/test/KL.png)

## References 
1. Shirts, Michael R, and John D Chodera. “Statistically optimal analysis of samples from multiple equilibrium states.” The Journal of chemical physics vol. 129,12 (2008): 124105. doi:10.1063/1.2978177
2. AJ Patel, P Varilly, D Chandler, and S Garde, "Quantifying Density Fluctuations in Volumes of All Shapes and Sizes using Indirect Umbrella Sampling", Journal of Statistical Physics, 145, 265 (2011).
3. AJ Patel, P Varilly, and D Chandler, "Fluctuations of Water Near Extended Hydrophobic and Hydrophilic Surfaces", Journal of Physical Chemistry B 114, 1632 (2010).
