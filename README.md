# AQ_curves

Fitting [multiple] photosynthetic light response curves in R [simultaneously].

In this repo you will find:

1) A script (file="AQ_curve_function.R") containing a function for fitting photosynthetic light response curves: `fit_AQ_curve()`.   

2) A script containing the function `diagnostic_AQ_plot()`. This takes the output from the fit_AQ_curve() along with the original curve data and produces a plot for each A-Q curve. Each plot has the fit curve plotted on top of the original A-Q points for visualizing how well the fit model matches the data. This has not yet (2017-11-11) had thorough testing, though does work for the curves I've thrown at it.   

3) An R-markdown doc (Introduction_to_AQ_curve_fitting.Rmd) that introduces the `fit_AQ_curve()` and `diagnostic_AQ_plot()` functions with an example analysis provided as an instructional guide for using the above functions.   

4) A directory (AQ_curves/ExampleDataFiles) with example light curves used in testing out the above functions.   

5) A directory (AQ_curves/OriginalPieces) with earlier pieces of the above functions before they were sufficiently functional.   

## Documentation

`fit_AQ_curve` {AQ_curves}                                                                        R Documentation

### Description

Fits a non-rectangular hyperbola model to photosynthetic light response gas exchange data and solves for the commonly used parameters: light-saturated photosynthetic rate (A_sat), quantum efficiency (Phi), mitochondrial respiration in the light (R_d), curvature/convexivity of light saturation (Theta), light compensation point (LCP), and the irradiance required to saturate photosynthesis.
    
### Usage    
    
`fit_AQ_curve(df, group_id, Photo, PARi, fit_type = "onls")`
    
### Arguments    

Argument      | Description
------------- | ----------------------------------------------------------------
`df`          | A data frame with your gas exchange data that contains a column with unique identifiers for each curve (e.g., "curve_1", "curve_2", etc), a column of net photosynthetic values, and a column of irradiance/Q/PAR/PPFD values.
`group_id`    | The name of the column containing curve identification values.
`Photo`       | The name of the column containing net photosynthetic values.
`PARi`        | The name of the column containing irradiance values.
`fit_type`    | The type of regression model to use for fitting. Defaults to orthogonal difference with `onls::onls()`. Can optionally switch to `base::nls()` by explicityly setting to "nls".

### Details   

The model fit is the standard non-rectangular hyperbola (e.g., see Lobo et al., 2013, _Photosynthetica_, v51, [doi.org/10.1007/s11099-013-0045-y][link to Lobo, et al., 2013]) of the form:   
A_N = ((Phi x Q + A_sat - sqrt[(Phi x Q + A_sat)^2 - 4 x Theta x Phi x Q x A_sat])/(2 x Theta)] - R_d  
where, A_N (µmol CO_2 m^-2 s^-1) is the net photosynthetic rate, Phi (mol CO_2 mol^-1 photons) is quantum efficiency, Q (µmol photons m^-2 s^-1) is irradiance in the range of photosynthetically active radiation, A_sat (µmol CO_2 m^-2 s^-1) is the light-saturated photosynthetic rate, Theta (unitless, range 0-1) is a convexivity parameter relating the curvature of the response, and R_d (µmol CO_2 m^-2 s^-1) is the respiration rate in the day/light.
NB: GitHub needs support for LaTeX equation rendering in markdown files.


### Value   

Returns an object of class data.frame that contains nine columns:  

1. The group_id's  
2. Light saturated net photosynthesis (Asat)  
3. Quantum yield (Phi)  
4. Mitochondrial respiration in the light (Rd)   
5. The curvature/convexivity factor of the curve (theta)  
6. The residual sum-of-squares (resid_SSs) from the model fit  
7. The light compensation point (LCP)  
8. PARi at 75% saturation of photosynthesis (Q_sat_75) - reliable  
9. PARi at 85% saturation (Q_sat_85) - much less reliable  


### Examples   
Default fitting with orthogonal distance regression, i.e., `onls::onls()`:   
`fit_AQ_curve(df = soy_AQ_curves, group_id = "curve_ids", Photo = "Photo", PARi = "PARi")`   

Fitting with standard least-squares non-linear regression, i.e., `base::bls()`   
`fit_AQ_curve(df = soy_AQ_curves, group_id = "curve_ids", Photo = "Photo", PARi = "PARi", fit_type = "nls")`   

[link to Lobo, et al., 2013]: http://doi.org/10.1007/s11099-013-0045-y
