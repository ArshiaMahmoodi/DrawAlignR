#// **********************************************************************************************
#//                         getPepLibData.R
#// **********************************************************************************************
#//
#//
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing


#' Extract data from a PQP library file
#' This function can be used to extract information from the PQP library file
#'
#' @param libfile A character vector of the absolute path and filename of the library file. (Must be .pqp format)
#' @param peptide_id A character vector for extraction of a specific peptide. I.e. 'ANSSPTTNIDHLK'
#' @return A data.table containing spectral library information
#'
#' @author Justin Sing \url{https://github.com/singjc}
#'
getPepLibData_ <- function( libfile, peptide_id = ''){
  # Load Requried Libraries
  library(dplyr)
  library(dbplyr)
  # Connect to database
  lib_db <- DBI::dbConnect( RSQLite::SQLite(), libfile )
  # Add Query statement to extract a specific peptide
  if ( peptide_id != '') {
    peptide_query <- sprintf( "INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID AND PEPTIDE.UNMODIFIED_SEQUENCE=('%s')", peptide_id )
  } else {
    peptide_query <- 'INNER JOIN PEPTIDE ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = PEPTIDE.ID'
  }
  # Create Query
  lib_query =  sprintf( "
SELECT TRANSITION.ID AS TRANSITION_ID,
    TRANSITION.TRAML_ID,
    TRANSITION.PRODUCT_MZ,
    TRANSITION.CHARGE,
    TRANSITION.TYPE,
    TRANSITION.ORDINAL,
    TRANSITION.DETECTING,
    TRANSITION.IDENTIFYING,
    TRANSITION.QUANTIFYING,
    TRANSITION.LIBRARY_INTENSITY,
    TRANSITION.DECOY,
    TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID,
    PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID,
    PEPTIDE.MODIFIED_SEQUENCE,
    PEPTIDE.UNMODIFIED_SEQUENCE,
    PRECURSOR.LIBRARY_RT AS PRECURSOR_LIBRARY_RT,
    PRECURSOR.PRECURSOR_MZ,
    PRECURSOR.CHARGE AS PRECURSOR_CHARGE
FROM TRANSITION
INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION.ID=TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID AND TRANSITION.DECOY=0
INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID=PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID
%s
INNER JOIN PRECURSOR ON PRECURSOR.ID=PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID", peptide_query)
  # Query Databasse
  df_lib <- collect( tbl(lib_db, sql(lib_query)) )
  # Disconnect from database
  DBI::dbDisconnect(lib_db)
  return( df_lib )
}
