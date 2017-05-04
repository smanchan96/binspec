A preprocessing/peak identification, classification, and feature ranking R library for mass spectra. The preprocessing algorithm is an extension of Yasui et al.'s peak-finding algorithm (https://www.biostat.wisc.edu/~kbroman/teaching/statgen/2004/refs/yasui.pdf). The classifiers available are radial kernel support vector machines (tuned with leave-one-out classification) and random forests (tuned with out-of-bag error), with the random forests providing a utility for ranking features based on mean decrease in Gini impurity.  We also provide a parallelized training function that allows for training classifiers on multiple sets of preprocessing parameters at once.

The package has uses the following libraries:
 * e1071 for support vector machines
 * randomForest for random forests
 * parallel for running tests with different parameter settings in parallel
 * roxygen2 and devtools for installing and documenting packages
which can be installed in the R console using the command:

```
install.packages(c("e1071", "randomForest", "parallel", "roxygen2", "devtools"))
```
In order to install the package in R on a unix-based system, run the following commands after installing roxygen2 and 
```
git clone https://github.com/smanchan96/binspec.git
cd binspec
./build.sh
```
If you are not on a unix-based system, you will not have access to the build script, so you will need to run build.R in the R console. 
Then you can get started in the R console by running
```
library(binspec)
```

Documentation is available [here](https://github.com/smanchan96/binspec/blob/master/binspec.pdf) and will automatically be available in your R console after installation, using `?function-name`.
