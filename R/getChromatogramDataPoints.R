#// **********************************************************************************************
#//                         getChromatogramDataPoints.R
#// **********************************************************************************************
#//
#//
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing


#' Extract Chromatographic data (Retention Time and Intensity) from mzML or sqMass chromatograms.
#' This function can be used to extract chromatogram retention and intensity values for a given
#' list of fragment ids.
#'
#' @param filename A character vector of the absolute path and filename of the chromatogram file. (Must be .mzML or sqMass format)
#' @param frag_ids A list of a vector containing fragment ids
#' @return A list of fragment ids containing 2 vectors for Retention time and Intensity
#'
#' @author Justin Sing \url{https://github.com/singjc}
#'
getChromatogramDataPoints_ <- function( filename, frag_ids ){
  ## Get File Extension Type
  fileType <- gsub( '.*\\.', '', filename)
  ## Extract Chromatogram Data
  if ( tolower(fileType)=='sqmass' ){
    # Read in an sqMass Chromatogram ------------------------------------------
    cat('Reading in chromatogram of ', blue$bold$underline('sqmass type.\n', sep=''))
    chrom <- getChromatogramsbyIndice_( filename, frag_ids )
    return(chrom)
  } else if ( tolower(fileType)=='mzml' ){
    # Read in an mzML chromatogram --------------------------------------------
    cat('Reading in chromatogram of ', blue$bold$underline('mzML type.\n', sep=''))
    # Load Required mzR library for reading data from mzML chromatograms
    library(mzR)
    # Create an mzR object that stores all header information, and use ProteoWizard api to access data from MzML file
    mz_object <- openMSfile(filename, backend = "pwiz", verbose = T)
    # Get header information for chromtagograms
    chromHead <- chromatogramHeader(mz_object)
    # Extract all the indices of chromatograms that match the transition names of the ones found in the TargetPetides file
    chromatogramIndices <- chromHead$chromatogramIndex[ match(frag_ids[[1]], chromHead$chromatogramId)  ]
    # Check how many chromatogramIndices are present to extract
    if ( length(chromatogramIndices)==1 ){
      rawChrom <- list(chromatograms(mz_object, chromatogramIndices))
    } else if ( length(chromatogramIndices)>1 ) {
      rawChrom <- chromatograms(mz_object, chromatogramIndices)
    } else {
      cat( red$bold$underline('There was no Chromatogramphic data for the following fragment(s): ', base::paste(frag_ids[[1]], collapse = ', ')), sep='')
    }
    chrom <- list(); rawChrom_idx <- 1
    for ( fragment_id in frag_ids[[1]] ){
      chrom[[ fragment_id ]] <- list(RT=rawChrom[[rawChrom_idx]][,1],
                                     Int=rawChrom[[rawChrom_idx]][,2])
      rawChrom_idx <- rawChrom_idx + 1
    }
    rm(mz_object)
    return(chrom)
  } else {
    cat( red$bold$underline(fileType, ' FileType is not supported!!\n'), sep='')
  }
} ## End Function
