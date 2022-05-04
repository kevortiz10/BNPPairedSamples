# BNPPairedSamples

## Overview

The objective of the BNP.test function is to perform a Bayesian nonparametric hypothesis test for paired samples proposed by Pereira et al. (2020), to determine if there are significant differences between the marginal distributions before and after an intervention or treatment performed on the study population. These differences can be further studied through the shift function, which allows us to identify in which range of values the differences between the two distributions lie, if any. Finally, the contour function allows us to visualize the joint distribution of the paired data, providing us to have a complete picture and to identify the correlation that may exist between the marginal distributions.


## Installation

You can install the released version of BNPPairedSamples from github (https://github.com/kevortiz10/BNPPairedSamples) with:

``` r
devtools::install_github("kevortiz10/BNPPairedSamples")
```

## Example

This section presents a brief example of the use of the BNP.test,  plotshift.function and the contours.plot functions.


``` r
library(BNPPairedSamples)
library(datarium)
library(dplyr)
data("weightloss")

weights1 <- weightloss[37:48,] %>% pull(t1)
weights2 <-weightloss[37:48,] %>% pull(t3)

```

For this example we will use the $datarium$ library and weight loss data, which correspond to data taken by a researcher on the weight of a set of sedentary men who were subjected to 4 trials to identify the effect of diet and exercise, the weight was taken at 3 time points, at the beginning, at the end and in the middle of each trial. For details, see https://cran.r-project.org/web/packages/datarium/datarium.pdf. Then, the initial time (t1) and the final time (t3) of trial 4 corresponding to the combination between diet and exercise were taken to identify if there are differences in weight.

Vector x, corresponding to the weights taken at the initial time (t1), vector y, corresponding to the weights taken at the final time (t3), and n.mcm, which refers to the number of simulations to be performed, are added. 10000 or more simulations are recommended to ensure better convergence in the final result.

``` r
results <- BNP.test(x =weights1 , y =weights2 ,n.mcm = 10000)
```
![](README-files/hypothesis_weight_loss.PNG)<!-- -->

![](README-files/marginal_weigh_loss.png)<!-- -->

Finally the function returns a list of 3 elements, the first element corresponds to the parameters obtained in the Gibbs sampling, the second element corresponds to the posterior probability for the alternative hypothesis and finally the third element corresponds to the data entered in the function that will be useful for the other functions of the package.

Calling the variable where the results of the function were saved and adding the $ sign extracts the element to be displayed.

``` r
results[[1]][1:2]
  
results$posterior.probability.H1
  
results$data.init
```

![](README-files/sampling_parameters_weight_loss.PNG)<!-- -->

![](README-files/posterior_probability_weight_loss.PNG)<!-- -->

![](README-files/raw_data_weight_loss.PNG)<!-- -->

After the test identifies differences in the marginal distributions, the shift function is run.


``` r
plotshift.function(results_BNP=results)
```
![](README-files/shift_weight_loss.png)<!-- -->

Finally, we run the contour function for a more complete picture.

``` r
contours.plot(results_BNP=results)
```
![](README-files/contours_weight_loss.png)<!-- -->

## References

Pereira, L. A., Taylor-Rodríguez, D. & Gutiérrez, L. (2020), A Bayesian nonparametric testing procedure for paired samples. Biometrics 76(1), 1-14.
