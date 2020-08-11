### nucdataBaynet - R package

The package `nucdataBaynet`, short for *nuclear data in a Bayesian network*, facilitates the setup and handling of data structures that contain experimental measurements, statistical and systematic errors and associated uncertainties or Gaussian process specifications.
In particular, handlers enable on the basis of these data structures the construction of the covariance matrix for the systematic components and the sensitivity matrix to map from systematic components to measured observables.
These mathematical objects can be used in the Bayesian generalized least squares method to update the systematic components using the measurements.
Also an implementation of a modified Levenberg-Marquardt (LM) algorithm is implemented to locate the posterior maximum exactly in the case of a non-linear model and multivariate normal prior and likelihood.
This implementation of the LM algorithm also exploits the data structures mentioned above for efficient computations.

## Use of this package

This version (0.0.0) of the package has been developed along a prototype of a nuclear data evaluation pipeline.
Usage examples can be found in the repository *eval-fe56* containing the scripts of the pipeline.
More specifically, the tuning of experimental normalization uncertainties in step four and the adjustment of the hyperparameters imposed on energy-dependent TALYS model parameters in step six, and the adjustment of Talys parameters in step seven relied on the functionality of this package. 

The interface of this package will probably change in the future.
This version is nevertheless kept for documentation purposes and to guarantee the reproducibility of the prototype pipeline.
