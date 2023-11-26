
# muchart <img src="logo.png" align="right" />

## Control Charts for Sample Mean of Gaussian, Inverse Gaussian, Gamma, and Beta Prime Distributions

### Description

<br>

Control Charts for Sample Mean of Gaussian, Inverse Gaussian, Gamma, and
Beta Prime Distributions. This package generates pseudo-random numbers
from the distributions described in the title in batches of samples. It
obtains maximum likelihood estimates that index the distributions of the
sample means of the respective distributions. It calculates control
limits and statistics for monitoring the performance of control charts,
such as ARL (Average Run Length), MRL (Median Run Length), and SDRL
(Standard Deviation Run Length).

It is possible to specify the true parameters if they are unknown or
simply not provided. In this case, the package obtains the estimates
using the maximum likelihood method. This is done using the following
functions:

1.  `stats_n()`: for the normal distribution;
2.  `stats_ig()`: for the inverse Gaussian distribution;
3.  `stats_ga()`: for the gamma distribution;
4.  `stats_bp()`: for the beta prime distribution.

For more details, please refer to the “Reference” tab of the package.

### Installation

<br>

You need install remotes package to install muchart package from github.

``` r
install.packages("remotes")
```

Then, you can install muchart package from github.

``` r
remotes::install_github("pedro-marinho/muchart")
```

### Authors

<br>

1.  [Pedro Rafael D. Marinho](https://prdm0.rbind.io/) (aut, cre)
2.  Luiz Medeiros de Araujo Lima Filho (aut)
3.  Marcelo Bourguignon Pereira (aut)
4.  Rodrigo Bernardo da Silva (aut)
