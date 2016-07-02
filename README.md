#binspec
A preprocessing/peak identification, classification, and feature ranking R library for mass spectra. The preprocessing algorithm is an extension of [Yasui et al.'s peak-finding algorithm] (https://www.biostat.wisc.edu/~kbroman/teaching/statgen/2004/refs/yasui.pdf). The classifiers available are radial kernel support vector machines (tuned with leave-one-out classification) and random forests (tuned with out-of-bag error), with the random forests providing a utility for ranking features based on mean decrease in Gini impurity.  We also provide a parallelized training function that allows for training classifiers on multiple sets of preprocessing parameters at once.  In order to install the package on your version of R, run the following commands.
```
git clone https://github.com/smanchan96/binspec.git
./build.sh
```
Additional dependencies may be necessary.
