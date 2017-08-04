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

NOTES - OA
- continue with predict.derivative!!!
- setseed_ME ???
- put examples
- GHmixedInit can be deleted - JW approved
- simulate functions - remove?? see LDMod in simulate function!
- estimateLong doesnt recognise silent?
- remove estimateME

- remove nglda_est - TEST FIRST THE REPLACER
- remove estimate.wrapper - TEST FIRST THE REPLACER
- remove nglda_predict - TEST FIRST THE REPLACER
