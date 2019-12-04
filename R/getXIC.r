#// **********************************************************************************************
#//                         getXIC.R
#// **********************************************************************************************
#//
#//
#// **********************************************************************************************
#// @Maintainer: Justin Sing
#// @Author: Justin Sing


#' Plot and Extracted Ion Chromatogram
#' This function can be used to plot an Extracted Ion Chromatogram
#'
#' @param graphic_obj A ggplot() graphics handle. Initialize with g <- ggplot() to create an empty ggplot handle
#' @param mod A character vector for specific peptide/modified peptide to extract I.e. 'ANSSPTTNIDHLK'/'ANS(UniMod:21)SPTTNIDHLK(UniMod:259)'. The MODIFIED_SEQUENCE column is used
#' @param df_lib A data.table containing spectral library information
#' @param chromatogram_file A character vector of the absolute path and filename of the chromatogram file. (Must be .mzML or sqMass format)
#' @param in_osw A character vector of the absolute path and filename of the OpenSwath Output file. (Must be .osw) @TODO maybe make this more robust for tsv files as well?
#' @param @TODO will add more decriptions for other params soon
#' @param
#' @param
#' @param
#'
#' @return A list containing graphic_obj = the graphic handle for the ggplot filled with data and max_Int = the maximun intensity
#'
#' @author Justin Sing \url{https://github.com/singjc}
#'
getXIC <- function( graphic_obj,
                    mod,
                    df_lib,
                    chromatogram_file,
                    in_osw=NULL,
                    transition_type='detecting',
                    intersecting_mz=NULL,
                    uni_mod_list=NULL,
                    max_Int=NULL,
                    smooth_chromatogram=TRUE,
                    doFacetZoom=FALSE,
                    FacetFcnCall=NULL,
                    top_trans_mod_list=NULL,
                    Isoform_Target_Charge=NULL,
                    RT_pkgrps=NULL,
                    plotIdentifying.Unique=NULL,
                    plotIdentifying.Shared=NULL,
                    plotIdentifying.Against=NULL,
                    show_n_transitions=NULL,
                    show_legend=T,
                    verbosity=0
){

  # filename <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/mzML_Chroms_Decomp/chludwig_K150309_013_SW_0_osw_chrom.mzML"

  # Helper Functions --------------------------------------------------------

  #' Check Dataframe, if empty return an object and stop function call
  #'
  #' @param df A data.frame/data.table/matrix object to check for number of rows
  #' @param return_item
  #' @param msg error message to return
  #' @return if data.frame has 0 rows, return true, otherwise return false
  #' @example
  #' checkDataframe( osw_df, graphic_obj, msg='There was no data found in OpenSwath Dataframe\n' )
  checkDataframe <- function( df, return_item, msg ){
    if ( dim(df)[1]==0 ){
      cat(bold(red(msg)))
      print(df)
      return( TRUE )
    } else {
      return( FALSE )
    }
  }

  #' Check numeric values if not NULL
  #'
  #' @param numeric_obj Check if an numeric object is not NULL a
  #' @return If numeric object is not null, round numeric object to 4 digits
  #' @example
  #' checkNumeric( osw_df_filtered$m_score[[1]] )
  checkNumeric <- function( numeric_obj ){
    if( !is.null( numeric_obj )){
      return( signif(numeric_obj, digits = 4) )
    } else {
      return( NULL )
    }
  }

  ## Filter library data for specific information on transition type selected
  if ( transition_type=='precursor' ){
    #cat( Verbose(threshold = verbosity), '--> Extracting Precursor Transition...\n')
    df_lib %>%
      dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
      dplyr::filter( TYPE=="" ) %>%
      dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge )-> df_lib_filtered
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found for precursor transition in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else if ( transition_type=='detecting_intersection' ){
    #cat( Verbose(threshold = verbosity), '--> Extracting Intersecting Detecting Transitions...\n')
    if ( !is.null( intersecting_mz ) ){
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( DETECTING==1 ) %>%
        dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
        dplyr::filter( PRODUCT_MZ %in% intersecting_mz ) -> df_lib_filtered
    } else {
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
        dplyr::filter( DETECTING==1 ) -> df_lib_filtered
    }
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found for common detecting transitions in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else if ( transition_type=='detecting_unique' ){
    #cat( Verbose(threshold = verbosity), '--> Extracting Unique Detecting Transitions...\n')
    df_lib %>%
      dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
      dplyr::filter( DETECTING==1 ) %>%
      dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
      dplyr::filter( !(PRODUCT_MZ %in% intersecting_mz) )-> df_lib_filtered
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found for unique detecting transitions in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else if ( transition_type=='detecting' ){
    #cat( Verbose(threshold = verbosity), '--> Extracting Unique Detecting Transitions...\n')
    df_lib %>%
      dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
      dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
      dplyr::filter( DETECTING==1 ) -> df_lib_filtered
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found detecting transitions in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else if ( transition_type=='identifying' ){
    #cat( Verbose(threshold = verbosity), '--> Extracting Identifying Transitions...\n')
    if ( !(is.null(top_trans_mod_list)) ){
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
        dplyr::filter( DETECTING==0 ) -> df_lib_filtered
    } else {
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
        dplyr::filter( DETECTING==0 ) -> df_lib_filtered

    }
    if ( checkDataframe( df_lib_filtered, graphic_obj, msg='There was no data found for identifying transitions in library\n' ) ){ return( list(graphic_obj=graphic_obj, max_Int=max_Int) ) }
  } else {
    ## OSW Information
    if ( !is.null( in_osw ) ){
      #cat( Verbose(threshold = verbosity), '--> Extracting OpenSwathResults Info...\n')
      df_lib %>%
        dplyr::filter( MODIFIED_SEQUENCE==mod ) %>%
        dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
        dplyr::filter( DETECTING==1 ) -> df_lib_filtered
      run_name <- gsub('_osw_chrom[.]sqMass$', '', basename(chromatogram_file))
      run <- gsub('_SW*|_SW_0|(*_-_SW[.]mzML[.]gz)', '', gsub('yanliu_I170114_\\d+_|chludwig_K150309_|lgillet_L\\d+_\\d+-Manchester_dirty_phospho_-_', '', run_name))
      mod_position <- getModificationPosition_(mod, character_index = T)+1
      mod_form_rename <- gsub('UniMod:4','Carbamidomethyl', gsub('UniMod:35','Oxidation', gsub('UniMod:259','Label:13C(6)15N(2)', gsub('UniMod:267','Label:13C(6)15N(4)', gsub('UniMod:21','Phospho', mod)))))
      ## Filter for precursor id that matches target charge state
      df_lib_filtered %>%
        dplyr::filter( PRECURSOR_CHARGE==Isoform_Target_Charge ) %>%
        select( PRECURSOR_ID ) %>%
        unique() %>% as.matrix() %>% as.numeric() -> target_charge_precursor
      ## @TODPO: Need to make this more robust for later
      if( mod==mod_form_rename ){
        osw_df <- getOSWData_( in_osw, run_name, precursor_id=target_charge_precursor, peptide_id='', mod_peptide_id='', mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=T )
        # # Original OSW Peptide Names
        # osw_pep_names <- gsub('UniMod:4','Carbamidomethyl', gsub('UniMod:35','Oxidation', gsub('UniMod:259','Label:13C(6)15N(2)', gsub('UniMod:267','Label:13C(6)15N(4)', gsub('UniMod:21','Phospho', osw_df$FullPeptideName)))))
        # # Keep only Rows that correspond to the correct Assay
        # osw_df %>% dplyr::filter( osw_pep_names == osw_df$ipf_FullPeptideName )
      } else {
        osw_df <- getOSWData_( in_osw, run_name, precursor_id=target_charge_precursor, peptide_id='', mod_peptide_id=c(mod,mod_form_rename), mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=T )
      }
      ## Check if openswath dataframe is empty
      if ( checkDataframe( osw_df, graphic_obj, msg='There was no data found in OpenSwath Dataframe\n' ) ){
        graphic_obj <- graphic_obj +
          ggtitle(  mod ) +
          labs(subtitle = paste('Run: ', run,
                                ' | Precursor: ', df_lib_filtered$PRECURSOR_ID,
                                ' | Peptide: ', df_lib_filtered$PEPTIDE_ID,
                                ' | Charge: ', df_lib_filtered$PRECURSOR_CHARGE, sep=''))
        return( list(graphic_obj=graphic_obj, max_Int=max_Int) )
      }
      if ( !is.null(Isoform_Target_Charge) ){
        osw_df %>%
          dplyr::filter( Charge==Isoform_Target_Charge ) -> osw_df
      }
      if ( !is.null( unlist(osw_df$ipf_pep) ) ){
        # Remove rows with NULL value in ipf_pep
        osw_df %>%
          dplyr::filter( !is.null(ipf_pep) ) %>%
          dplyr::filter( !is.nan(ipf_pep) ) -> osw_df

        osw_df %>%
          dplyr::filter( ipf_pep==min(ipf_pep) ) -> osw_df_filtered #### No longer filter by peak group
      } else {
        osw_df %>%
          dplyr::filter( peak_group_rank==1 ) -> osw_df_filtered #### No longer filter by peak group
      }
      m_score <- checkNumeric( osw_df_filtered$ms2_m_score[[1]] )
      ipf_pep <- checkNumeric( osw_df_filtered$ipf_pep[[1]] )
      ipf_m_score <- checkNumeric( osw_df_filtered$m_score[[1]] )
      ms2_pkgrp_rank <- checkNumeric( osw_df_filtered$peak_group_rank[[1]] )

      graphic_obj <- graphic_obj +
        geom_vline(xintercept = osw_df_filtered$RT, color='red', size = 1.5 ) +
        geom_vline(xintercept = osw_df_filtered$leftWidth, color='red', linetype='dotted', size=1 ) +
        geom_vline(xintercept = osw_df_filtered$rightWidth, color='red', linetype='dotted', size=1 ) + # default size 0.65
        ggtitle(  mod ) +
        labs(subtitle = paste(
          # 'Run: ', run,
          # ' | Precursor: ', df_lib_filtered$PRECURSOR_ID,
          # ' | Peptide: ', df_lib_filtered$PEPTIDE_ID,
          'Charge: ', osw_df_filtered$Charge,
          ' | m/z: ', osw_df_filtered$mz,
          ' | RT: ', osw_df_filtered$RT,
          # '\nlib_RT: ', round(osw_df_filtered$assay_iRT, digits = 4),
          ' | ms2-m-score: ', m_score,
          ' | q_value: ', ipf_m_score,
          ' | ipf_pep: ', ipf_pep,
          ' | ms2_pkgrp_rank: ', ms2_pkgrp_rank,
          sep=''))
      ### Plot Other peak rank groups
      if( !is.null(RT_pkgrps) ){
        osw_RT_pkgrps <- getOSWData_( in_osw, run_name, precursor_id=target_charge_precursor, peptide_id='', mod_peptide_id=c(mod,mod_form_rename), mod_residue_position='', peak_group_rank_filter=F, pep_list='', mscore_filter='', ipf_filter='', ms2_score=T, ipf_score=F )

        osw_RT_pkgrps %>%
          dplyr::filter( RT %in% RT_pkgrps ) %>%
          dplyr::filter( RT != osw_df_filtered$RT ) %>%
          dplyr::select( RT, leftWidth, rightWidth, peak_group_rank ) -> osw_RT_pkgrps_filtered
        if ( dim(osw_RT_pkgrps_filtered)[1]!=0 ){
          # Define unique set of colors to annotate different peak rank groups
          jBrewColors <- brewer.pal(n = dim(osw_RT_pkgrps_filtered)[1], name = "Dark2")
          # Y intcrements
          y_increment = 0
          for( RT_idx in seq(1, dim(osw_RT_pkgrps_filtered)[1],1) ){
            if ( any(osw_RT_pkgrps_filtered$RT[RT_idx] %in% c(5691.37, 5737.26)) ){ cat('Skipping: ', osw_RT_pkgrps_filtered$RT[RT_idx], '\n', sep=''); next }
            point_dataframe <- data.frame(RT=(osw_RT_pkgrps_filtered$RT[RT_idx]),
                                          y=((max(max_Int)/ 4 )-y_increment),
                                          # y=((1000)-y_increment),
                                          label=paste('Rank:',osw_RT_pkgrps_filtered$peak_group_rank[RT_idx],'\n',osw_RT_pkgrps_filtered$RT[RT_idx],sep=' '))
            graphic_obj <- graphic_obj +
              geom_vline(xintercept = osw_RT_pkgrps_filtered$RT[RT_idx], color=jBrewColors[RT_idx], alpha=0.65, size = 1.5 ) +
              geom_vline(xintercept = osw_RT_pkgrps_filtered$leftWidth[RT_idx], color=jBrewColors[RT_idx], linetype='dotted', alpha=0.85, size=1 ) +
              geom_vline(xintercept = osw_RT_pkgrps_filtered$rightWidth[RT_idx], color=jBrewColors[RT_idx], linetype='dotted', alpha=0.85, size=1 ) +
              # geom_rect( aes(xmin=osw_RT_pkgrps_filtered$leftWidth[RT_idx], xmax=osw_RT_pkgrps_filtered$rightWidth[RT_idx], ymin=0, ymax=Inf), col=jBrewColors[RT_idx], fill=jBrewColors[RT_idx], alpha=0.5) +
              geom_label(data=point_dataframe, aes(x=RT, y=y,label=label), alpha=0.7, fill=jBrewColors[RT_idx], size=3)
            y_increment = y_increment + 500
          }
        }
      }

      if ( doFacetZoom==TRUE ){
        if ( !is.null(FacetFcnCall) ){
          graphic_obj <- graphic_obj + FacetFcnCall
        } else {
          if ( max(max_Int) > 1000 ){
            graphic_obj <- graphic_obj +
              # facet_zoom(ylim = c(0, (1000) ))
              facet_zoom(ylim = c(0, (max(max_Int)/ 4 ) ))
            # facet_zoom(ylim = c(0, (max(max_Int)/ (max(max_Int)-mean(max_Int)) ) ))
          }
        }
      }
      if ( !is.null(top_trans_mod_list) ){
        graphic_obj <- graphic_obj +
          theme(plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5, size = 10),
                legend.text = element_text(size = 2),
                legend.key.size = unit(0.5, "cm")) +
          # theme( panel.background = element_rect( fill='lightblue', color='lightblue', size=0.5, linetype='solid'),
          #       panel.border = element_rect( fill=NA, color='black', size=0.5, linetype='solid'),
          #       panel.grid.major = element_line( color='white', size=0.5, linetype='solid'),
          #       panel.grid.minor = element_line( color='white', size=0.25, linetype='solid') )
          theme_bw()

        if ( show_legend==FALSE ){
          graphic_obj <- graphic_obj + theme(legend.position="none")
        }

      } else {
        graphic_obj <- graphic_obj +
          theme(plot.title = element_text(hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5, size = 10)) +
          # theme(panel.background = element_rect( fill='#c1e1ec', color='#c1e1ec', size=0.5, linetype='solid'),
          #       panel.border = element_rect( fill=NA, color='black', size=0.5, linetype='solid'),
          #       panel.grid.major = element_line( color='white', size=0.5, linetype='solid'),
          #       panel.grid.minor = element_line( color='white', size=0.25, linetype='solid') )
          # guides( Transition=FALSE ) +
          theme_bw()
        if ( show_legend==FALSE ){
          graphic_obj <- graphic_obj + theme(legend.position="none")
        }
      }

      return( list(graphic_obj=graphic_obj, max_Int=max_Int) )
    } else {
      cat(red( bold(underline(transition_type)), ' is not a supported argument for transition_type!!!\n'), sep='')
      return( list(graphic_obj=graphic_obj, max_Int=max_Int) )
    }
  }

  ## Get Transition IDs for chromatogram data extraction
  #cat( Verbose(threshold = verbosity), '----> Getting Transition IDs and Extracting Chromatogram Data...\n')
  frag_ids <- list()
  if( length(as.character( df_lib_filtered$TRANSITION_ID ))>1 ){
    frag_ids[[1]] <- as.character( df_lib_filtered$TRANSITION_ID )
  } else {
    frag_ids[[1]] <- list(as.character( df_lib_filtered$TRANSITION_ID ))
  }

  #***************************************#
  #***    Extract Chromatogram Data    ***#
  #***************************************#

  chrom <- getChromatogramDataPoints_( chromatogram_file, frag_ids  )

  if (smooth_chromatogram==TRUE){
    #cat( Verbose(threshold = verbosity), '----> Smoothing Chromatogram Data...\n')
  }
  ## Smooth Intensity values to make Chromatogram look nice
  for (i in seq(1:length(chrom))){
    names(chrom[[i]]) <- c('RT','Int')
    if (smooth_chromatogram==TRUE){
      chrom[[i]]$Int <- signal::sgolayfilt(chrom[[i]]$Int, p = 4, n = 9)
    }
  }

  df_plot <- bind_rows(mclapply(chrom, data.frame), .id='Transition')

  if ( transition_type=='precursor' ){
    df_lib_filtered %>%
      dplyr::filter( TRANSITION_ID %in% unique(df_plot$Transition) ) %>%
      select( TRANSITION_ID, PRECURSOR_CHARGE ) -> transition_info
    transition_ids <- paste( paste('0 - ',transition_info$TRANSITION_ID,sep=''), paste(transition_info$PRECURSOR_CHARGE,'+',sep=''), sep='_')
    tmp <- sapply(seq(1,length(df_plot$Transition)), function(i){ transition_ids[grepl(paste('0 - ',df_plot$Transition[i],'_*',sep=''), transition_ids)] } )
  } else {
    # if ( transition_type=='identifying' ){
    #   if( !(is.null(top_trans_mod_list)) ){
    #     cat( '----> Extracting Top Transitions with low PEPs..\n')
    #     df_plot %>%
    #       dplyr::filter( Transition %in% (top_trans_mod_list[[mod]]$transition_id[top_trans_mod_list[[mod]]$transition_pep<1])[1:10] ) -> df_plot
    #   }
    # }
    df_lib_filtered %>%
      dplyr::filter(  TRANSITION_ID %in% unique(df_plot$Transition) ) %>%
      select( TRANSITION_ID, CHARGE, TYPE, ORDINAL, PRODUCT_MZ ) -> transition_info
    transition_ids <- paste( transition_info$TRANSITION_ID, paste(transition_info$CHARGE,'+',sep=''), paste(transition_info$TYPE, transition_info$ORDINAL,sep=''), transition_info$PRODUCT_MZ, sep='_')
    tmp <- sapply(seq(1,length(df_plot$Transition)), function(i){ transition_ids[grepl(paste('^',df_plot$Transition[i],'_*',sep=''), transition_ids)] } )
  }
  df_plot$TRANSITION_ID <- df_plot$Transition
  df_plot$Transition <- tmp

  ## Check Max Intensity
  if ( transition_type!='precursor' & dim(df_plot)[1]>0 ){
    # if ( max(df_plot$Int) > max_Int ){ max_Int <- max(df_plot$Int) }
    max_Int <- c(max_Int, max(df_plot$Int[is.finite(df_plot$Int)]))
  }

  ## Plotting Action
  if ( transition_type=='precursor' ){
    graphic_obj <- graphic_obj +
      geom_line( data=df_plot, aes(RT, Int, group=Transition), show.legend = T, alpha=0.65, linetype='solid', col='black' )
  } else if ( transition_type=='detecting_intersection' ){
    graphic_obj <- graphic_obj +
      geom_line(data=df_plot, aes(RT, Int, group=Transition), show.legend = F, alpha=0.95, col='gray')
  } else if ( transition_type=='detecting_unique' ){
    graphic_obj <- graphic_obj +
      geom_line(data=df_plot, aes(RT, Int, group=Transition), alpha=0.75, col='black', linetype='dotted', show.legend = F)
    # theme( legend.position = 'bottom' )
  } else if ( transition_type=='detecting' ){
    graphic_obj <- graphic_obj +
      geom_line(data=df_plot, aes(RT, Int, group=Transition), alpha=0.75, col='gray', linetype='solid', show.legend = F)
  } else{

    df_plot <- merge( df_plot, select( df_lib_filtered[df_lib_filtered$TRANSITION_ID %in% unique(df_plot$TRANSITION_ID),], c(TRANSITION_ID, TRAML_ID) ), by='TRANSITION_ID' )

    if ( !(is.null(top_trans_mod_list)) ){

      #### Unique Transitions
      if ( plotIdentifying.Unique==TRUE ){
        df_plot %>%
          dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')),
                               gsub('.*\\{|\\}.*','',TRAML_ID)) &
                           (!grepl(glob2rx('*\\|*'),
                                   gsub('.*\\{|\\}.*','',TRAML_ID))) ) %>%
          dplyr::filter( TRANSITION_ID %in% (top_trans_mod_list[[mod]]$transition_id[top_trans_mod_list[[mod]]$transition_pep<1]) ) -> tmp_plot

        # df_plot %>%
        #   dplyr::filter( grepl(glob2rx("*_\\{ESTAEPDSLS(Phospho)R(Label:13C(6)15N(4))\\}_*"), TRAML_ID) ) -> tmp_plot
        #
        #
        # df_plot %>%
        #   dplyr::filter( grepl(glob2rx("*_\\{ESTAEPDS(Phospho)LSR(Label:13C(6)15N(4))\\}_*"), TRAML_ID) ) -> tmp_plot

        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }

        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int

        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot



        if ( dim(tmp_plot)[1] > 0){
          graphic_obj <- graphic_obj +
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition), linetype='solid', alpha=0.5, size=1.5, show.legend = show_legend) #+
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
          #

          # graphic_obj + geom_line(data=tmp_plot, aes(RT, Int, col=Transition), linetype='solid', alpha=0.5, size=1.5, show.legend = show_legend) + facet_zoom(ylim = c(0, 5000))

        }
      }
      #### Shared Transitions
      if ( plotIdentifying.Shared==TRUE ){
        # if ( transition_type=='identifying' ){
        #   if( !(is.null(top_trans_mod_list)) ){
        #     df_plot <- df_plot_org
        #     df_plot %>%
        #       dplyr::filter( TRANSITION_ID %in% (top_trans_mod_list[[mod]]$transition_id[top_trans_mod_list[[mod]]$transition_pep<1])[1:10] ) -> df_plot
        #   }
        # }
        df_plot %>%
          dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')),
                               gsub('.*\\{|\\}.*','',TRAML_ID)) &
                           grepl(glob2rx('*\\|*'),
                                 gsub('.*\\{|\\}.*','',TRAML_ID)) ) %>%
          dplyr::filter( TRANSITION_ID %in% (top_trans_mod_list[[mod]]$transition_id[top_trans_mod_list[[mod]]$transition_pep<0.6]) ) -> tmp_plot

        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }

        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int

        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot

        if ( dim(tmp_plot)[1] > 0){
          graphic_obj <- graphic_obj +
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition), linetype='F1', alpha=0.5, show.legend = show_legend) #+
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
      ##### Transitions against
      if ( plotIdentifying.Against==TRUE ){
        df_plot %>%
          dplyr::filter( !grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')),
                                gsub('.*\\{|\\}.*','',TRAML_ID)) ) %>%
          dplyr::filter( TRANSITION_ID %in% (top_trans_mod_list[[mod]]$transition_id[ top_trans_mod_list[[mod]]$transition_pep<1 ]) ) -> tmp_plot

        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }

        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int

        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot

        if ( dim(tmp_plot)[1] > 0){
          graphic_obj <- graphic_obj +
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition), linetype='dashed', show.legend = show_legend) #+
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
    } else {
      #### Unique Identifying Transitions
      if ( plotIdentifying.Unique==TRUE ){
        df_plot %>%
          dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')),
                               gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID))) &
                           (!grepl(glob2rx('*\\|*'),
                                   gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID)))) ) -> tmp_plot


        mods_present <- names(uni_mod_list)
        alternate_mod <- mods_present[ !( mods_present %in% mod)]
        alternate_mod_index <- unlist(lapply(alternate_mod, getModificationPosition_))
        current_mod_index <- getModificationPosition_(mod)

        identification_traml_ids <- str_split( gsub('.*\\{|\\}.*', '', df_plot$TRAML_ID), '\\|' )
        current_mod_transitions <- unlist(lapply( identification_traml_ids, function( traml_id ){
          mods_for_ident_transition <- unique(unlist(lapply(traml_id, getModificationPosition_)))
          if ( any( current_mod_index == mods_for_ident_transition ) ){
            if ( any( (mods_for_ident_transition %in% alternate_mod_index) ) ){
              return( FALSE )
            } else {
              return( TRUE )
            }
          } else {
            return( FALSE )
          }
        } ) )


        df_plot %>%
          dplyr::filter( current_mod_transitions ) -> df_plot
        tmp_plot <- df_plot

        # len_unmod <- nchar( unique(df_lib$UNMODIFIED_SEQUENCE) )
        # uni_mod_list <- getSiteDeterminingIonInformation_( mod, len_unmod )
        # ion_ordinal_type <- gsub( '*_\\d+\\.\\d+', '', gsub('.*\\d+_\\d[\\+]_*', '', df_plot$Transition))

        # df_plot %>%
        #   dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')),
        #                        gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID))) &
        #                    ( gsub('[[:digit:]]', '', ion_ordinal_type) =='y' & as.numeric(gsub('[^[:digit:]]', '', ion_ordinal_type)) < as.numeric(uni_mod_list$y[4]) ) |
        #                    ( gsub('[[:digit:]]', '', ion_ordinal_type) =='b' & as.numeric(gsub('[^[:digit:]]', '', ion_ordinal_type)) < as.numeric(uni_mod_list$b[3]) )
        #                  ) -> tmp_plot

        # df_plot %>%
        #   dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')),
        #                        gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID))) &
        #                    ( gsub('[[:digit:]]', '', ion_ordinal_type) =='y' & as.numeric(gsub('[^[:digit:]]', '', ion_ordinal_type)) %in% c(2,3) ) |
        #                    ( gsub('[[:digit:]]', '', ion_ordinal_type) =='b' & as.numeric(gsub('[^[:digit:]]', '', ion_ordinal_type)) %in% c(8,9) )
        #   ) -> tmp_plot


        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }

        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int

        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot

        if ( dim(tmp_plot)[1] > 0){

          graphic_obj <- graphic_obj +
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition), linetype='solid', alpha=0.5, size=1.5, show.legend = show_legend) #+
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
      ##### Shared Identifying Transitions
      if ( plotIdentifying.Shared==TRUE ){
        df_plot %>%
          dplyr::filter( grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')),
                               gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID))) &
                           grepl(glob2rx('*\\|*'),
                                 gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',df_plot$TRAML_ID))) ) -> tmp_plot

        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }

        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int

        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot

        if ( dim(tmp_plot)[1] > 0){

          graphic_obj <- graphic_obj +
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition), linetype='F1', show.legend = show_legend) #+
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }

      ##### Transitions Against
      if ( plotIdentifying.Against==TRUE ){
        df_plot %>%
          dplyr::filter( !grepl(glob2rx( paste('*', gsub('\\(UniMod:259\\)|\\(UniMod:267\\)|\\(Label.*)','', gsub('UniMod:35', 'Oxidation', gsub('UniMod:4', 'Carbamidomethyl', gsub('UniMod:21','Phospho',mod)))) ,'*',sep='')),
                                gsub('\\(Label:13C\\(6\\)15N\\(4\\)\\)', '', gsub('.*\\{|\\}.*','',TRAML_ID)) ))  -> tmp_plot

        # Number of Transitions to display
        if ( !is.null(show_n_transitions) ){
          if ( show_n_transitions==-1 ){
            show_n_transitions_val <- length(unique(tmp_plot$TRAML_ID))
          } else {
            show_n_transitions_val <- show_n_transitions
          }
        } else {
          show_n_transitions_val <- 6
        }

        tmp_plot %>%
          group_by(TRANSITION_ID) %>%
          top_n(n=1, wt=Int) %>%
          arrange(desc(Int)) %>%
          select(TRANSITION_ID) %>%
          head(n=show_n_transitions_val) -> Ordered_top_Int

        tmp_plot %>%
          dplyr::filter( TRANSITION_ID %in% as.matrix(Ordered_top_Int) ) -> tmp_plot

        if ( dim(tmp_plot)[1] > 0){
          graphic_obj <- graphic_obj +
            geom_line(data=tmp_plot, aes(RT, Int, col=Transition), linetype='dashed', show.legend = show_legend) #+
          # theme(legend.text = element_text(size = 2),
          #       legend.key.size = unit(0.5, "cm"))
        }
      }
    }
  }

  return( list(graphic_obj=graphic_obj, max_Int=max_Int) )
}
