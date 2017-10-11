# AQ_curves
Fitting [multiple] photosynthetic light response curves in R [simultaneously].

In this repo you will find:
1) update on 10/06/2017, a script (file="AQ_curve_function.R") containing the model fitting as a function. I recommend using this version as it is much more user friendly.
2) two data files to test it out, one with simulated data (file="simulated_LRCs.RDS") and one with real world curves (file="athal_LRCs.RDS"). Find a more complete descriptions of these below.
3) a script (file="Fitting_AQ_curves.R") walking through the fitting of multiple A-Q curves. This script is commented throughout with explainations of what is what. 

# Which script to use?
If you just want to quickly fit some curves, use "AQ_curve_function.R" This rewrites the "Fitting_AQ_curves.R" script as a function, making it more user friendly.
If you want to play around with the fitting code, use "Fitting_AQ_curves.R" It is essentially the internals of the "AQ_curve_function.R" and is therefore more readily modifiable - and should break in a more sensical/understandable manner. 

# Provided data files
(1)
The "simulated_LRCs.RDS" file consists of 40 simulated curves created with the _Photo_ function (below) of the model of A-Q curves and altering each of the input parameters as:
- PhiCO2: (range 0.05 to 0.1 by 0.01), default == 0.9
- Asat:   (range 10 to 40 by 2), default == 38.0
- theta:  (range 0.2 to 0.9 by 0.1), default == 0.6
- Rd:     (range 0.5 to 1.4 by 0.1), default == 1.0
- PARi:  Photo solved with the parameters above along the PARi sequence of c(0, seq(10, 100, 10), seq(150, 400, 50), 500, 600, 800, 1000, 1250, 1500, 1800)

using the function: NB- you can copy/paste this function and simulate your own curves.
Photo = function(PhiCO2, PARi, Asat, theta, Rd){ 
      # Function to simulate photosynthetic values from input parameters and a
      #     range PARi values
      ((
            PhiCO2 * PARi + Asat - 
                  sqrt((PhiCO2 * PARi + Asat)^2 - 4 *
                             PhiCO2 * theta * PARi * Asat)
      ) / (2*theta) - Rd)
}

(2)
The "athal_LRCs.RDS" file has eleven A-Q curves measured on Arabidopsis thaliana ecotypes. These were used as preliminary data to figure out how high to set PARi when measuring A-Ci curves. These curves only have the PARi and Photo columns from the original LiCor file plus a column with a designator for which ecotype was measured, a w/in ecotype replicate id, a designator for whether it has (relative to other ecotypes) high or low leaf mass per area (irrelavent here - I should have deleted this, sorry!), and a column with ecotype and replicate pasted together to produce a unique identifier for each of the curves --> this is what the AQ_curve_function iterates on and you need something that can serve this function in your data file(s).
