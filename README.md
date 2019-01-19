# Instruction for developer #

* Before pushing make sure that works the following test can be run:
	1.	test/mixedeffect/test_NIGMixedEffect.R

# TODO #

## structure ##
* Add a file that automatically run all test files that needs to be run be fore pushing to bitbucket.

## Methods ##
* Add joint sampling of Gaussian components?
* Automatic switching between distributions based on TV distance
* Add support for weighted subsampling
* Add support for individual nu

## Interface ##
* Add reasonable starting values
* Add checks of input dimensions to all functions
* Add checks of initial values in estimation

## TO DO -- OA
- shape parameter is only 0.5 currently - future work
- all the params of gig can be estimated - future work

## NOTES -- OA
- make model specifications flexible, e.g. Normal, normal, Gaussian are all acceptable
- allow different locs in predict than those from the observed data
- continue with predict.derivative!!!
- setseed_ME - explanations to be added
- put examples
- simulate functions - remove?? see LDMod in simulate function!
- pSubsample2 seems un-used
- consider ~ -1 + ... in fixed input of ngme 
- add mean or median as an option to the predict.ngme function
- check calculation of coverages and width in predict.ngme
- plot for excursions
- put a note that if nu reaches the limit use normal
- newdata in predict
- fit.init is expected to be gaussian atm
- controls.init vs fit.init
- process nu never becomes 0.001 for normal nig normal model
- warning nu = 100 atm, not for 99.99
- get rid of the for loop at the end of predict.ngme
- ylabs in plot.ngme
- add forecasting and prediction  by fine interval by locs.pred
- rename columns and rows of fisher_est in ngme.fisher 