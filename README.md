# AQ_curves
Fitting [multiple] photosynthetic light response curves in R [simultaneously].

In this repo you will find:
1) update on 10/06/2017, an additional script (file="AQ_curve_function.R") containing the model fitting as a function. I recommend using this version as it is much more user friendly.
2) two data files to test it out, one with simulated data (file="simulated_LRCs.RDS") and one with real world curves (file="athal_LRCs.RDS"). Find a more complete descriptions of these below.
3) an R script (file="Fitting_AQ_curves.R") walking through the fitting of multiple A-Q curves. This script is commented throughout, explaining as you go. 

# Which script to use?
If you just want to quickly fit some curves, use "AQ_curve_function.R" This is a rewrites the "Fitting_AQ_curves.R" script as a function, making it more user friendly.
If you want to play around with the fitting code, use "Fitting_AQ_curves.R" It is essentially the internals of the "AQ_curve_function.R" and is therefore more readily modifiable - and should break in a more sensical/understandable manner. 

# Provided data files
The "simulated_LRCs.RDS" file consists of 40 simulated curves created with the _Photo_ function (below) of the model of A-Q curves and altering each of the input parameters as:
_Parameter_ _defaultValue_  _alteredValues_
PhiCO2  0.9 (0.05 to 0.1 by 0.01
Asat  38.0  (10 to 40 by 2)
theta 0.6 (0.2 to 0.9 by 0.1)
Rd  1.0 (0.5 to 1.4 by 0.1)
PARi  NA  c(0, seq(10, 100, 10), seq(150, 400, 50), 500, 600, 800, 1000, 1250, 1500, 1800)
