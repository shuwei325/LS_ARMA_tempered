# Two step estimation of locally stationary ARMA processes with tempered stable innovations

This is the accompanying GitHub repository to a paper by [Shu Wei Chou Chen](https://shuwei325.github.io/) and [Pedro A. Morettin](https://www.ime.usp.br/~pam/). 

[![licensebuttons
by](https://licensebuttons.net/l/by-nc/4.0//88x31.png)](https://creativecommons.org/licenses/by/4.0)


# Abstract

The class of locally stationary processes assumes a time-varying (tv) spectral representation and finite second moment. Different areas have observed phenomena with heavy tail distributions or infinite variance. Using the stable distribution as a heavy-tailed innovation is an attractive option. However, its estimation is difficult due to the absence of a closed expression for the density function and non-existence of second moment. In this paper, we propose the tvARMA model with tempered stable innovations, which have lighter tails than the stable distribution and have finite moments. A two-step method is proposed to estimate this parametric model. In the first step, we use the blocked Whittle estimation to estimate the time-varying structure of the process. In the second step, we recover residuals from the first step and use the maximum likelihood method to estimate the rest of the parameters related to the standardized classical tempered stable (stdCTS) innovations. We perform simulation studies to evaluate the consistency for the maximum likelihood estimation of independent stdCTS samples. Then, we execute simulations to study the two-step estimation of our model. Finally, an empirical application is illustrated.
	 
# Resources

## Code

* Code to run the application section from the paper can be found in this [link](https://github.com/shuwei325/LS_ARMA_tempered/tree/main/application).
* Code to run two simulations included in the paper can be found in this [link](https://github.com/shuwei325/LS_ARMA_tempered/tree/main/simulation).


## Data

* The dataset source of the application is [EMHIRES dataset](http://dx.doi.org/10.2790/831549).

