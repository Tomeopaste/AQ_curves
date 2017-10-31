################################################################################
######## Function to plot fits from AQ_curve_function.                  ########
########    Contact me @Tomeopaste or TomeoNJ@gmail.com with questions. ########
################################################################################

# This plotting function will create plots for diagnosing how well the fits
#  from fit_AQ_curve match your data. No ideal fit metrics exist for a non-
#  linear curve fit like this, and none that I am aware of come close to the
#  utility of just looking at the fits relative to your data. This function 
#  is meant to give you a quick and easy way to assess how well the
#  parameters of the model fit your data.
#
#
# Will probably break if any of the fits produced NA's. I need to test this and
#   troubleshoot a solution if that is indeed the case.
#
#
# Run this scritpt. Then to use the function, provide as input:
#   1) the original dataframe with the light response curves ("curve_data"), 
#   2) the dataframe of curve_fits output by fit_AQ_curve() ("fit_data"), 
#   3) the column headers for Photo, PARi, and group_id from the curve_data
#         dataframe just as with the fit_AQ_curve() function
# AND optionally:
#   4) save_to_pdf = TRUE which will put all of the plots produced into a 
#         pdf file at the location you specify with
#   5) save_path = "/Users/you/your_file_path" and the name you specify with
#   6) file_name = "my_file_name"
#   Including 4, 5, & 6 yields /Users/you/your_file_path/my_file_name.pdf with
#         all of the plots. 
#
#
# When you are running the function to output plots you will have something
#   like this:
# diagnostic_AQ_plot(curve_data = "my_AQ_curves_df", 
#                    fit_data = "my_fit_AQ_curves_output_df",
#                    Photo = "Photo", PARi = "PARi", group_id = "group_ids")
# OR:
#
# diagnostic_AQ_plot(curve_data = "my_AQ_curves_df", 
#                    fit_data = "my_fit_AQ_curves_output_df",
#                    Photo = "Photo", PARi = "PARi", group_id = "group_ids",
#                    save_to_pdf = TRUE,
#                    save_path = "/Users/me/R/LightCurves/diagnosticPlots/",
#                    file_name = "ProjectSunShade_AQ_plots")

diagnostic_AQ_plot <- function(curve_data, fit_data, Photo, PARi, group_id,
                               save_to_pdf = FALSE, save_path, file_name){
      if(save_to_pdf == FALSE){ 
            par(mar = c(3, 3, 1, 1), oma = c(1, 1, 1, 1))
            for(i in seq_along(1:length(unique(curve_data[[group_id]])))){
                  single_curve <- 
                        curve_data[curve_data[[group_id]] == 
                                         unique(curve_data[[group_id]])[i],]
                  plot(
                        single_curve[[Photo]] ~ single_curve[[PARi]] ,
                        xlim = c(-2, max(curve_data[[PARi]])), 
                        ylim = c(min(curve_data[[Photo]]) - 2,
                                 max(curve_data[[Photo]]) + 2),
                        pch = 3,
                        cex = 2,
                        xlab = "",
                        ylab = "",
                        main = paste("Data from curve ",
                                     as.character(
                                           unique(single_curve[[group_id]])))
                  )
                  mtext(expression("Photo (µmol "*CO[2]*" "*m^-2*" "*s^-1*")"),
                        line = 2.4, side = 2)
                  mtext(expression("PARi (µmol photons "*m^-2*" "*s^-1*")"),
                        line = 2.4, side = 1)
                  par(new = TRUE)
                  curve(((
                        fit_data$Phi[i] * PARi + fit_data$Asat[i] - 
                              sqrt((fit_data$Phi[i] * PARi + fit_data$Asat[i])^2 - 4 *
                                         fit_data$Phi[i] * fit_data$theta[i] * PARi *
                                         fit_data$Asat[i])
                  ) / (2*fit_data$theta[i]) - fit_data$Rd[i]),
                  from = 0, to = 1600, 
                  xname = "PARi",
                  xlab = "", ylab = "", 
                  xlim = c(-2, max(curve_data[[PARi]])), 
                  ylim = c(min(curve_data[[Photo]]) - 2,
                           max(curve_data[[Photo]]) + 2),
                  axes = FALSE,
                  col = "red",
                  lwd = 2
                  )
            }} else{
             if(dir.exists(save_path)){
      pdf(paste0(save_path, file_name, ".pdf"))
      par(mar = c(3, 3, 1, 1), oma = c(1, 1, 1, 1))
      for(i in seq_along(1:length(unique(curve_data[[group_id]])))){
            single_curve <- 
                  curve_data[curve_data[[group_id]] == 
                                   unique(curve_data[[group_id]])[i],]
            plot(
                  single_curve[[Photo]] ~ single_curve[[PARi]] ,
                  xlim = c(-2, max(curve_data[[PARi]])), 
                  ylim = c(min(curve_data[[Photo]]) - 2,
                           max(curve_data[[Photo]]) + 2),
                  pch = 3,
                  cex = 2,
                  xlab = "",
                  ylab = "",
                  main = paste("Data from curve ",
                               as.character(
                                     unique(single_curve[[group_id]])))
            )
            mtext(expression("Photo (µmol "*CO[2]*" "*m^-2*" "*s^-1*")"),
                  line = 2.4, side = 2)
            mtext(expression("PARi (µmol photons "*m^-2*" "*s^-1*")"),
                  line = 2.4, side = 1)
            par(new = TRUE)
            curve(((
                  fit_data$Phi[i] * PARi + fit_data$Asat[i] - 
                        sqrt((fit_data$Phi[i] * PARi + fit_data$Asat[i])^2 - 4 *
                                   fit_data$Phi[i] * fit_data$theta[i] * PARi *
                                   fit_data$Asat[i])
            ) / (2*fit_data$theta[i]) - fit_data$Rd[i]),
            from = 0, to = 1600, 
            xname = "PARi",
            xlab = "", ylab = "", 
            xlim = c(-2, max(curve_data[[PARi]])), 
            ylim = c(min(curve_data[[Photo]]) - 2,
                     max(curve_data[[Photo]]) + 2),
            axes = FALSE,
            col = "red",
            lwd = 2
            )
            }
      dev.off()
             } else {
            return(
                  "Warning: the file path provided to save_path does not exist"
            )}
            
}
}
