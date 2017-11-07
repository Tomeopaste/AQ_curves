# AQ_curves
Fitting [multiple] photosynthetic light response curves in R [simultaneously].

In this repo you will find:
1) A script (file="AQ_curve_function.R") containing a function for fitting photosynthetic light response curves: fit_AQ_curve().

2) A script containing the function diagnostic_AQ_plot(). This takes the output from the fit_AQ_curve() along with the original curve data and produces a plot for each A-Q curve. Each plot has the fit curve plotted on top of the original A-Q points for visualizing how well the fit model matches the data. This has not yet (2017-11-01) had thorough testing, though appears to work for the curves I've thrown at it thus far.

3) A directory (AQ_curves/ExampleDataFiles) with example light curves used in testing out the above functions.

4) A directory (AQ_curves/OriginalPieces) with earlier pieces of the above functions before they were sufficiently functional.

5) An R-markdown doc (<FileNameGoesHere>) with an example analysis provided as an instructional guide for using the above functions. ***I need to make this (11/07/17)***
