################################################################################
######## Script to iterate through photosynthetic light response (A-Q)  ########
########    curve fitting. Written by Nick Tomeo, September 2017, and   ########
########    heavily credited to, modified from, a script created by     ########
########    JM Heberling and widely distributed.                        ########
########    Contact me if you have questions or run into trouble:       ########
########    @Tomeopaste ... TomeoNJ@gmail.com                           ########
################################################################################

# To begin with, I am inlcuding two 'datasets' here:
#     - The first consists of 40 simulated curves created with a function
#           (below) of the model of A-Q curves and altering each of the input
#           parameters as:
#            - PhiCO2:   0.9   (range 0.05 to 0.1 by 0.01)
#            - Asat:     38.0  (range 10 to 40 by 2)
#            - theta:    0.6   (range 0.2 to 0.9 by 0.1)
#            - Rd:       1.0   (range 0.5 to 1.4 by 0.1)
#            - PARi: defined as PARi_short = c
#                 (0, seq(10, 100, 10), seq(150, 400, 50), 500, 600, 800, 1000, 
#                        1250, 1500, 1800)

Photo = function(PhiCO2, PARi, Asat, theta, Rd){ 
      # Function to simulate photosynthetic values from input parameters and a
      #     range PARi values
      ((
            PhiCO2 * PARi + Asat - 
                  sqrt((PhiCO2 * PARi + Asat)^2 - 4 *
                             PhiCO2 * theta * PARi * Asat)
      ) / (2*theta) - Rd)
}
# 
#     - The second has eleven A-Q curves measured on Arabidopsis thaliana
#           ecotypes. These were used as preliminary data to figure out how
#           high to set PARi when measuring A-Ci curves.
#
# Either of these 'datasets' can be used to test/trial the code below. The 
#     simulated curves of course produce near perfect fits.

sim_curves = readRDS("simulated_LRCs.RDS") # The simulated curves
athal.lrc = readRDS("athal_LRCs.RDS") # Arabidopsis curves

# Set up an empty data_frame to dump the fit-values for each curve into:
fits.sim <- data.frame(ID = character(), #      group_id
                   Asat = numeric(),
                   PhiCO2 = numeric(),
                   Rd = numeric(),
                   theta = numeric(),
                   resid_SSs = numeric(), #     resudual sum-of-squares for the
                   #    curve, maybe not the best fit statistic, but 
                   #    better than none...
                   stringsAsFactors = FALSE # (Oh, R...!)
)
fits.athal <- data.frame(ID = character(),
                       Asat = numeric(),
                       PhiCO2 = numeric(),
                       Rd = numeric(),
                       theta = numeric(),
                       resid_SSs = numeric(),
                       LCP = numeric(), # add light compensation point
                       stringsAsFactors = FALSE
)

# Iterate through the simulated curves:
for(i in seq_along(unique(sim_curves$curveID))){ # hit each group_id iteratively
      fits.sim[i, 1] <- unique(sim_curves$curveID)[i] # Column with group_ID's
      # Subset by group_ID:
      single_curve <- subset(
            sim_curves[sim_curves$curveID == unique(sim_curves$curveID)[i],]
            )
      # Initial estimate of PhiCO2 based on slope of A~PAR with the 1st 5 
      #     points, we'll use this below as the starting point for fitting
      phico2.as.slope <- with(single_curve,
                              as.numeric(coef(lm(Photo[1:5] ~ PARi[1:5]))[2]))
      # Fit the curve:
      temp.fit <- with(single_curve, # use the single curve subset
                       nls(Photo ~ ((PhiCO2 * PARi + Asat - # The model
                                           sqrt((PhiCO2 * PARi + Asat)^2 - 
                                                      4 * PhiCO2 * 
                                                      theta * Asat * PARi ))
                                    )/(2 * theta) - Rd, 
                           start=list( # Specify starting values parameter fits
                                 Asat = (max(Photo)), # maximum A in data
                                 PhiCO2 = phico2.as.slope, # Slope from above
                                 Rd = -min(Photo), # (+) lowest A value
                                 theta = 0.5), # A midrange value ... 
                           #  open to suggestions for a better theta option
                           control = list(maxiter = 50), # Up the max iterations
                           algorithm = "port") # Changing the algorithm. In my
                       #      experience the default fails more often.
      )
      fits.sim[i, 2] <- as.numeric(coef(temp.fit)[1]) # Pull out the model 
      fits.sim[i, 3] <- as.numeric(coef(temp.fit)[2]) #     coefficients and      
      fits.sim[i, 4] <- as.numeric(coef(temp.fit)[3]) #     dump them in 
      fits.sim[i, 5] <- as.numeric(coef(temp.fit)[4]) #     the data.frame
      fits.sim[i, 6] <- sum(resid(temp.fit)^2) #      And resid_SS's
}
fit.sim # NB: the first column here, the "ID" specifies which parameter was 
#     varied during the simulation and it's modified value, e.g., Asat12 means
#     that Asat had a value of 12 µmol m^-2 s^-1. Default values in the 
#     simulations are above, and are Asat==38, PhiCO2==0.9, theta==0.6, Rd==1. 


# Doing the same for the Arabidopsis curves as an example of using "real data"
head(athal.lrc)
# These curves only have the PARi and Photo columns from the original LiCor file
#     plus a column with a designator for which ecotype was measured, a w/in 
#     ecotype replicate id, a designator for whether it has (relative to other
#     ecotypes) high or low leaf mass per area, and a column with ecotype and 
#     replicate pasted together to produce a unique identifier for each of the
#     curves --> this is what we'll iterate on and you need something that can
#     function this way in your data file(s).
for(i in seq_along(unique(athal.lrc$ids))){
      fits.athal[i, 1] <- unique(athal.lrc$ids)[i] # Column with group_ID's
      # Subset by group_ID iteratively:
      single_curve <- subset(
            athal.lrc[athal.lrc$ids == unique(athal.lrc$ids)[i],]
            )
      # Initial estimate of PhiCO2 based on slope of A~PAR for the 1st 5 points:
      phico2.as.slope <- with(single_curve,
                              as.numeric(coef(lm(Photo[1:5] ~ PARi[1:5]))[2]))
      # Fit the curve:
      temp.fit <- with(single_curve, # use the subset of a single curve
                       nls(Photo ~ ((PhiCO2 * PARi + Asat - 
                                           sqrt((PhiCO2 * PARi + Asat)^2 - 
                                                      4 * PhiCO2 * theta * 
                                                      Asat * PARi ))
                                    )/(2*theta) - Rd,
                           start=list(
                                 Asat = (max(Photo)),
                                 PhiCO2 = phico2.as.slope,
                                 Rd = -min(Photo),
                                 theta = 0.5),
                           control = list(maxiter = 50),
                           algorithm = "port")
      )
      fits.athal[i, 2] <- as.numeric(coef(temp.fit)[1]) 
      fits.athal[i, 3] <- as.numeric(coef(temp.fit)[2])
      fits.athal[i, 4] <- as.numeric(coef(temp.fit)[3])
      fits.athal[i, 5] <- as.numeric(coef(temp.fit)[4])
      fits.athal[i, 6] <- sum(resid(temp.fit)^2)
      fits.athal[i, 7] <- ( # This is just the non-rectangular model 
            # solved for Q, and Photo is dropped since at the light
            # compensation point, net assimilation equals zero
            ((as.numeric(coef(temp.fit)[3])) * 2 * 
                   as.numeric(coef(temp.fit)[4]) + as.numeric(coef(temp.fit)[2])
             + as.numeric(coef(temp.fit)[1]))^2 +
                  (as.numeric(coef(temp.fit)[2]) + 
                         as.numeric(coef(temp.fit)[1]))^2
            ) / (
                  4 * as.numeric(coef(temp.fit)[4]) * 
                        as.numeric(coef(temp.fit)[1]))
}
fits.athal
# Notice the resid_SSs are substantially higher with real data...
#     These fits mostly make sense though - I should (and may eventually) 
#     compare them with fits obtained elsewhere using the old Excel sheet 
#     methods - but generally these parameter values are about what I'd expect.
#####
#####
# So that's it. Provide a data.frame with Photo and PARi values plus a curve 
#       identifier column. Create a data.frame to put the paramater estimates 
#       in, and point the first four lines of the for-loop toward your 
#       data.frame. 
