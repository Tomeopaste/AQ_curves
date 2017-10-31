################################################################################
######## Function to plot fits from AQ_curve_function.                  ########
########    Contact me @Tomeopaste or TomeoNJ@gmail.com with questions. ########
################################################################################

# provided at this point mostly without any comment. will add to this later.
# outputs a plot for each fit curve. will probably break if any of the fits 
#     produced NA's.
# Input the original curve dataframe, the curve_fits dataframe, and the column 
#     headers for Photo, PARi, and group_id from the curve dataframe


diagnostic_AQ_plot <- function(curve_data, fit_data, Photo, PARi, group_id){
      par(mar=c(3,3,1,1), oma=c(1,1,1,1))
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
            mtext(expression("Photo (µmol photons "*m^-2*" "*s^-1*")"),
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
}
