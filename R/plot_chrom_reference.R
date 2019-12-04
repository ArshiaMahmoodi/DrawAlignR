
#' A function for plotting the ms2 chromatogram of any peptide. It's usage is to plot the reference chromatogram,
#' which is the run that all other chromatograms are aligned to.
#'
#' @param chrom rds chromatogram file which contains a list of peptides with their retention times and intensities.
#'
#' @param peptide the precursor id of the specific peptide you would like to plot.
#'
#' @param Run_Id The run from which the peptide came from
#'
#' @param RT The retention time of the peptide
#'
#' @param Left_width Left width of the retention time
#'
#' @param Right_width Right width of the retention time
#'
#' @param mz The mass to charge value of the precursor
#'
#' @param sequence The sequence of the precursor
#
#' @return None. Produces plot output.
#'

plot_chrom_reference <- function(chrom, precursor, Run_ID, RT, Left_width, Right_width, mz, sequence) {
  
  precursor_string <- toString(precursor)
  
  plot(chrom$'Peptide'[[1]]$time, chrom$'Peptide'[[1]]$X73440, type='l',
       main= paste("Ms2 Chromatogram of Precursor",precursor_string,"Run", Run_ID, "\n",  sequence, ",", "M/Z= ", mz),
       xlab="Retention time (s)", ylab="Intensity",
       font.main=2, font.lab=4, col.main='blue',
       xlim= c(2350,2550))
  
  lines(chrom$'Peptide'[[2]]$time, chrom$'Peptide'[[2]]$X73523, col='red')
  lines(chrom$'Peptide'[[3]]$time, chrom$'Peptide'[[3]]$X73532, col='blue')
  lines(chrom$'Peptide'[[4]]$time, chrom$'Peptide'[[4]]$X73580, col='yellow')
  lines(chrom$'Peptide'[[5]]$time, chrom$'Peptide'[[5]]$X73635, col='purple')
  lines(chrom$'Peptide'[[6]]$time, chrom$'Peptide'[[6]]$X73664, col='green')
  
  abline(v=RT, lty= 2, col = 'red', lwd = 2)
  abline(v=Left_width , lty= 2)
  abline(v=Right_width, lty = 2)
  
}
