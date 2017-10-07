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
#                  Photo = "Photo", 
#                  PARi = "PARi", 
#                  group_id = "group_id")

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
            AQ_curve_fits[i, 1] <- unique(df[[group_id]])[i]
            # Subset by group_ID iteratively:
            single_curve <- df[df[[group_id]] == unique(df[[group_id]])[i],]
            single_curve$assim <- single_curve[[Photo]]
            single_curve$PAR <- single_curve[[PARi]]
            phi.as.slope <- with(single_curve,
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
                  as.numeric(coef(temp.fit)[3]) * as.numeric(coef(temp.fit)[4] - 
                                                                   as.numeric(coef(temp.fit)[1]))
            )) / (as.numeric(coef(temp.fit)[2]) * (
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
      }
      return(AQ_curve_fits)
}
