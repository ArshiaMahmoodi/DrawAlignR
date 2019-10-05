#' A function for the visualization of an aligned chromatogram to a reference chromatogram
#'
#' @param chrom rds chromatogram file which contains a list of peptides with their retention times and intensities.
#'
#' @param peptide the precursor id of the specific peptide you would like to plot.
#'
#' @param shift The specific retention time shift of the alignment, obtained upstream using DialignR
#'
#' @param Run_Id The run from which the peptide came from
#'
#' @param Reference_Id The reference run/precursor to which this peptide was aligned to.
#'
#' @param Reference_RT The retention time of the precursor in the reference run.
#'
#' @param RT The retention time of the peptide
#'
#' @param Left_width Left width of the retention time
#
#' @param Right_width Right width of the retention time
#'
#' @param mz The mass to charge value of the precursor
#'
#' @param The sequence of the precursor
#'
#' @return None. Produces plot output.
#'
#' @export

plot_aligned <- function(chrom, shift, precursor, Run_ID, Reference_ID, RT, Reference_RT,
                                 Left_width, Right_width, mz, sequence) {

  precursor_string <- toString(precursor)

  plot(chrom$'960'[[1]]$time, chrom$'960'[[1]]$X73440, type='l', col='gray',
       main= paste("Aligned Ms2 Chromatogram of Precursor",precursor_string, "Run", Run_ID, "\n",  sequence, ",", "M/Z= ", mz),
       xlab="Retention time (s)", ylab="Intensity",
       font.main=2, font.lab=4, col.main='blue',
       xlim= c(2350,2550),
       sub= paste("Shift =", shift, ", Error =", RT - shift - Reference_RT), col.sub='purple', font.sub=2)

  lines(chrom$'960'[[2]]$time, chrom$'960'[[2]]$X73523, col='gray')
  lines(chrom$'960'[[3]]$time, chrom$'960'[[3]]$X73532, col='gray')
  lines(chrom$'960'[[4]]$time, chrom$'960'[[4]]$X73580, col='gray')
  lines(chrom$'960'[[5]]$time, chrom$'960'[[5]]$X73635, col='gray')
  lines(chrom$'960'[[6]]$time, chrom$'960'[[6]]$X73664, col='gray')

  abline(v=RT, lty= 2, col = 'black', lwd = 2)


  lines(chrom$'960'[[1]]$time - shift, chrom$'960'[[1]]$X73440, col='black')
  lines(chrom$'960'[[2]]$time - shift, chrom$'960'[[2]]$X73523, col='red')
  lines(chrom$'960'[[3]]$time - shift, chrom$'960'[[3]]$X73532, col='blue')
  lines(chrom$'960'[[4]]$time - shift, chrom$'960'[[4]]$X73580, col='yellow')
  lines(chrom$'960'[[5]]$time - shift, chrom$'960'[[5]]$X73635, col='purple')
  lines(chrom$'960'[[6]]$time - shift, chrom$'960'[[6]]$X73664, col='green')

  abline(v=RT - shift, lty= 2, col = 'black', lwd = 2)

  abline(v=Reference_RT, col = 'green', lwd = 2)

  arrows(RT, 3000, RT-shift, 3000, length = 0.20, angle = 30, code = 2, col = 'purple', lwd = 3)

}
