# A Bayesian Approach to Estimating Abundance from Belt Transect Surveys

The rkscallop folder contains all code necessary to run the model found Duskey, Hart, Chang, and Sullivan 2022.  Briefly, the code can simulate density of units of interest over one- and two-dimensional coordinates, here imagined as depth and space, respectively.  The estimation model uses regression kriging in a Bayesian framework in order to estimate parameters of a generalized additive model (GAM) and a variogram with joint probability.

## Abstract

Populations of marine organisms are often characterized by pronounced spatial patterning.  This occurs at different scales, which may be driven by separate mechanisms.  Regression kriging is often used to generate fish population estimates which may parse variation among these scales.  A previous application of this and other frequentist methods to estimate the abundance of Atlantic sea scallops (Placopecten magellanicus) found that semi-parametric regression kriging performs best among model-based methods.  We expanded upon these results by adapting the model structure for a Bayesian framework.  This approach demonstrated the trade-off between abundance as a function of depth and broad spatial habitat features, and of fine-scale spatial aggregation.  We applied our model to the 2015 belt transect survey of Atlantic sea scallops in Georges Bank and the Mid-Atlantic Bight off the eastern coast of the United States.  The overall Bayesian predictions generally agreed with those of the frequentist method.  However, the former favored broad spatial trends coupled with very fine-scale aggregation, while the latter favored finer trends and broad aggregation as the dominant driver of scallop abundance.  Therefore, one can draw similar conclusions about abundance, yet disparate conclusions about the drivers.  Careful consideration of the benefits and limitations of each approach is warranted.

## Data

The actual data cannot be made available in a public database, though it is available upon request.  However, the code provided simulate data of similar character.  Many parameters may be manipulated to simulate data that vary widely in nature.  Generally, we used Gaussian functions to approximate peaks in numerical abundance of organisms in space and over and one-dimensional covariate that we refer to as depth.  We also used standard methods to generate autocorrelated residuals.  These were all multiplied to produce means with which we simulated Poisson counts.

## Requirements

The Bayesian model contained in the code is run in Stan via the R package cmdstanr.  See the following website for instructions on how to install cmdstanr and Stan:

https://mc-stan.org/cmdstanr/

## Usage

Follow these steps prior to first time usage to ensure the code runs properly on your machine:

1. Download the "rkscallop-main" zipped folder from Github in its entirety
2. Unzip and place "rkscallop-main" in your preferred directory
3. Navigate to the Code folder and open the file "depend.R" in R or RStudio
4. Confirm that all listed packages are installed on your machine; uncomment and run lines corresponding to packages you have not yet installed, then re-comment
5. Change the "mypath" variable to the file path which contains the rkscallop folder
6. Save and close "depend.R"

Upon all subsequent uses, we recommend running the .R files in the following order:

1. depend.R -- to set mypath; alternatively, change setwd() parameters in RUN.R to avoid this step
2. RUN.R -- to simulate data and run the model

All files not appearing here are sourced in one or more of the files listed above, and typically contain functions necessary to run the code in each script.

NOTE: The empty Data and Output folders contain a file called ".keep".  This is meant to contain and organize your own output.  Feel free to delete the ".keep" files after downloading.

## Contact us

If you have questions or concerns regarding this code, or would like help in re-formatting it for your own use, please do not hesitate to contant the corresponding author at:

elizabeth.duskey@slu.se

## License

Copyright 2022 Elizabeth Duskey

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
