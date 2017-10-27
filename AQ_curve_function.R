################################################################################
######## Function to fit multiple photosynthetic light response (A-Q)   ########
########    curves. Written by Nick Tomeo, October 2017. The model      ########
########    code was originally modified from a script created by JM    ########
########    Heberling.                                                  ########
########    Contact me @Tomeopaste or TomeoNJ@gmail.com with questions. ########
################################################################################

# The function below will fit a photosynthetic light response curve to a non-
#       rectangular hyperbola model. Supply the function with 
#       1) a dataframe containing 
#             i) a column of PAR values, 
#             ii) a column of net photosynthesis values,
#             iii) and a column with a group/curve identifier values
#       2) the names of the above three columns 
#
# And the function will return a data frome that contains:
#     1) The group_id's
#     2) Light saturated net photosynthesis (Asat)
#     3) Quantum yield (Phi)
#     4) Mitochondrial respiration in the light (Rd)
#     5) The curvature/convexivity factor of the curve (theta)
#     6) The model residual sum-of-squares (resid_SSs)
#     7) The light compensation point (LCP)
#     8) PARi at 75% saturation of photosynthesis (Q_sat_75) - reliable
#     9) PARi at 85% saturation (Q_sat_85) - much less reliable
#
#
# Run the below function code and then you can fit your curves like so:
#
#  fit_AQ_curve(my_df_of_AQ_curves, 
#                  Photo = "Photo", # Enter the name of column of A_n vals
#                  PARi = "PARi", # Same, but for irradiance values
#                  group_id = "group_id") # Same, but for the group_id column

fit_AQ_curve <- function(df, group_id, Photo, PARi){
      AQ_curve_fits <- data.frame(ID = character(),
                                  Asat = numeric(),
                                  Phi = numeric(),
                                  Rd = numeric(),
                                  theta = numeric(),
                                  resid_SSs = numeric(),
                                  LCP = numeric(),
                                  Q_sat_75 = numeric(),
                                  Q_sat_85 = numeric(),  
                                  stringsAsFactors = FALSE
      )
      for(i in seq_along(unique(df[[group_id]]))){
            tryCatch({
            AQ_curve_fits[i, 1] <- unique(df[[group_id]])[i]
            # Subset by group_ID iteratively:
            single_curve1 <- df[df[[group_id]] == unique(df[[group_id]])[i],]
            single_curve1$assim <- single_curve1[[Photo]]
            single_curve1$PAR <- single_curve1[[PARi]]
            # Reorder to assure PAR goes from low to high, allowing for the 
                  # preliminary phi.as.slope estimate two lines below.
            single_curve <- single_curve1[order(single_curve1$PAR),]
            phi.as.slope <- with(single_curve, # use as starting value
                                    as.numeric(coef(lm(
                                          assim[1:5] ~ PAR[1:5]))[2]))
            # Fit the curve:
            temp.fit <- with(single_curve, # use the subset of a single curve
                             nls(assim ~ ((Phi * PAR + Asat - 
                                                 sqrt((Phi * PAR + Asat)^2 - 
                                                            4 * Phi * theta * 
                                                            Asat * PAR ))
                             )/(2*theta) - Rd,
                             start=list(
                                   Asat = (max(assim)),
                                   Phi = phi.as.slope,
                                   Rd = -min(assim),
                                   theta = 0.5),
                             control = list(maxiter = 50),
                             algorithm = "port")
            )
            AQ_curve_fits[i, 2] <- as.numeric(coef(temp.fit)[1]) # asat 
            AQ_curve_fits[i, 3] <- as.numeric(coef(temp.fit)[2]) # Phi
            AQ_curve_fits[i, 4] <- as.numeric(coef(temp.fit)[3]) # Rd
            AQ_curve_fits[i, 5] <- as.numeric(coef(temp.fit)[4]) # theta
            AQ_curve_fits[i, 6] <- sum(resid(temp.fit)^2)
            AQ_curve_fits[i, 7] <- (as.numeric(coef(temp.fit)[3]) *(
                  as.numeric(coef(temp.fit)[3]) * as.numeric(coef(temp.fit)[4]) - 
                                                                   as.numeric(coef(temp.fit)[1]))
            ) / (as.numeric(coef(temp.fit)[2]) * (
                  as.numeric(coef(temp.fit)[3]) - as.numeric(coef(temp.fit)[1])
            ))
            AQ_curve_fits[i, 8] <- (
                  (as.numeric(coef(temp.fit)[1]) * 0.75 + 
                         (as.numeric(coef(temp.fit)[3]))) * (
                               as.numeric(coef(temp.fit)[1]) * 0.75 *
                                     as.numeric(coef(temp.fit)[4]) +
                                     as.numeric(coef(temp.fit)[3]) *
                                     as.numeric(coef(temp.fit)[4]) -
                                     as.numeric(coef(temp.fit)[1])
                         )) / (
                               as.numeric(coef(temp.fit)[2])* (
                                     as.numeric(coef(temp.fit)[1]) * 0.75 +
                                           as.numeric(coef(temp.fit)[3]) -
                                           as.numeric(coef(temp.fit)[1])
                               ))
            
            AQ_curve_fits[i, 9] <- (
                  (as.numeric(coef(temp.fit)[1]) * 0.85 + 
                         (as.numeric(coef(temp.fit)[3]))) * (
                               as.numeric(coef(temp.fit)[1]) * 0.85 *
                                     as.numeric(coef(temp.fit)[4]) +
                                     as.numeric(coef(temp.fit)[3]) *
                                     as.numeric(coef(temp.fit)[4]) -
                                     as.numeric(coef(temp.fit)[1])
                         )) / (
                               as.numeric(coef(temp.fit)[2])* (
                                     as.numeric(coef(temp.fit)[1]) * 0.85 +
                                           as.numeric(coef(temp.fit)[3]) -
                                           as.numeric(coef(temp.fit)[1])
                               ))
      }, error = function(E){cat("Error: ", conditionMessage(E), "\n")})
      }
      return(AQ_curve_fits)
}


### Details ###
# The non-rectangular hyperbola is the standard, widely-used, A-Q model. As fit
      # here this is essentially the same as Equation 6 in Lobo et al. 2013 (
      # Photosynthetica, v51, doi.org/10.1007/s11099-013-0045-y),
      # Except, with the caveat that I am unsure how/why they modified the 
      # equation to solve for light saturation values. Because I cannot follow 
      # their logic, and it is not explicitly explained (or cited), I have 
      # chosen a slightly different tactic. Namely, the non-rectangular 
      # hyperbola model rearrannged for PARi solves to:
      #       PARi = ( (A_n + Rd) * (A_n * theta + Rd * theta - Asat)/
      #                     (Phi * (A_n + Rd - Asat)) )
      # To solve for PARi at photosynthetic saturation ("Q_sat_75" and
      # "Q_sat_85") I have entered entered Asat for A_n and reduced it 
      # proportionally, i.e., Q_sat_75 uses 75% of the fit Asat (Asat * 0.75)
      # and Q_sat_85 uses 85% of the fit Asat (Asat * 0.85). The values I have 
      # gotten out for Q_sat_75 make sense and appear resonable. The values for
      # Q_sat_85 are less so: they're often absurdly high, e.g., above 
      # maximum solar PAR at earth's surface. If you are fitting AQ curves to
      # figure out what light levels to use for measuring A_n under other 
      # circumstances (A-Ci curves) and want to measure at 'saturating' yet not
      # photoinhibitory PARi, my recommendation is to use Q_sat_75 and then add 
      # 20% or so to it. For example, if Q_sat_75 is on average 1000 µmol m^-2 
      # s^-1; adding 20% you should measure at 1200 µmol m^-2 s^-1.
