#' LTMPDataTools.R is a utility library for filtering, processing and grouping LTMP data
#'
#' @description 
#' Data processing functions to deal with the LTMP (Long Term Monitoring 
#' Program) and MMP (marine Monitoring Program) Data from AIMS (Australian Institute
#' of Marine Sciences).
#'
#'
#' @section Current features:
#' \itemize{
#'   \item data processing for splitting data into recovery trajectories from 
#'     site-level using LTMP disturbance classifications.
#'   \item Filtering of recovery trajectories based on upper-bound for initial 
#'     percent hard coral cover, lower-bound for final percent hard coral cover, and lower 
#'     on number of observations.
#'   \item Enables the assignment of trjectory groups base on external file or k-means clustering 
#' }
#' @section Warning:
#' All functions process trajectories to the site level, however, functions require transect level
#' data to enumerate the standard error of the site observations and disturbance significance testing. 
#' @section Important disclaimer:
#' This code assumes *.Rdata files of the form provied by Juan Ortiz and Angus Thompson, modifications \emph{may} be
#' required if data tables differ from this. We have done our best to make codes as generic as possible, but occasionally,
#' assumptions may have been made along the way (particularily with respect to table column names).
#' @author
#' \itemize{
#'      \item David J. Warne[1,2,3] (\email{david.warne@qut.edu.au})
#'      \item Grace E. M. Heron[1,3] (\email{g.heron@qut.edu.au})
#' }
#' \enumerate{
#'  \item School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#'  \item Centre for Data Science, Queensland University of Technology
#'  \item ARC Centre of Excellence for Mathematical and Statistical Frontiers (ACEMS)
#' }
#'
#' @note Functions are specifically focused on the filtering of AIMS coral cover 
#' data at a recovery trajectory level (sequences between two disturbance periods
#' for a given site). 
#' @docType data
#' @name AA_Library_Summary
NULL

# ROXYGEN_STOP


# ROXYGEN_START

#-------------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------------

#' Extract recovery trajectories by transect
#' 
#' @description  Extraction of recovery trajectories using transect level data.
#' Uses differences in successive transect observations at each site over
#' time to extract to detect start and end of site level recovery trajectories. 
#' Assigns a recovery id to each trajectory at a site level.
#'
#' @param disturbance        data frame of disturbance records
#' @param samples             data frame of visit time/place records
#' @param transect.cover     data frame of percent cover at transect level, can be at
#'                           major group level or benthic code (ker code) level
#'                          by major group ("GROUP_CODE") or ker code ("BENTHOS_CODE")
#' @param code               code base trajectory tests on
#' @param alpha              confidence level for paired t-test
#'
#' @note currently only supports GROUP_CODE grouping, future versions should 
#' implement BENTHOS_CODE grouping (but will probably need to specify groups to 
#' aggregate)
#' @note consider alternative approach to disturbance detection, e.g., a mixture
#' of t-test and disturbance records 
#' @family Data extraction
extract_recovery_trajectories_transect <- function(disturbance,samples,
                                                   transect.cover,
                                                   code = "HC",alpha = 0.95) {

    # to avoid factor warnings
    combined <- sort(union(union(levels(disturbance$REEF),levels(samples$REEF)),
                           levels(transect.cover$REEF)))
    disturbance <- disturbance %>% 
        mutate(REEF=factor(REEF,levels=combined))
    samples <- samples %>% 
        mutate(REEF=factor(REEF,levels=combined))
    transect.cover <- transect.cover %>% 
        mutate(REEF=factor(REEF,levels=combined))
   
    # extend disturbance factor levels to distinguish between recorded distubances
    # and our own entries
    combined <- sort(union(levels(disturbance$DISTURBANCE), c('U','T','S')))
    disturbance <- disturbance %>% 
        mutate(DISTURBANCE=factor(DISTURBANCE,levels=combined))


    # join disturbances and samples to transects
    samples.disturbance.transect.cover <- right_join(
                      right_join(disturbance%>%select(-c("P_CODE")),samples,by = c("REEF","VISIT_NO")),
                      select(transect.cover, P_CODE,REEF,DEPTH,VISIT_NO,SITE_NO),
                      by = c("REEF","DEPTH","VISIT_NO")) %>%
                      mutate(RP_ID = 0,DIFF = 0, Pval = 0)
    print(samples.disturbance.transect.cover)
    #return(samples.disturbance.transect.cover)
    recovery.trajectories.site <- data.frame()
    # for every REEF/SITE/DEPTH triple by visit
    unique_sites <- unique(select(samples.disturbance.transect.cover,
                                   REEF,DEPTH,SITE_NO))
    rp_id <- 1
    for (i in 1:nrow(unique_sites)){
        visit.sequence.site <- samples.disturbance.transect.cover %>% 
                                  filter(REEF == unique_sites$REEF[i],
                                         DEPTH == unique_sites$DEPTH[i],
                                         SITE_NO == unique_sites$SITE_NO[i]) %>%
                                  unique() %>%
                                  arrange(VISIT_NO)
        # step through sequence applying one-sided paired t-test to subsequent pairs
        num_dist <- 1
        visit.seq <- visit.sequence.site$VISIT_NO
        for (j in 1:length(visit.seq)) {
            # first visit is always a new recovery trajectory
            if (j == 1) { 
                rp_id <- rp_id + 1
                num_dist <- 1
                visit.sequence.site$DISTURBANCE[j] <- 'S' 
            } else { 
                # Get transect cover for level code provided for time t and t-1
                transect.cover.pair <- transect.cover %>% 
                          filter(REEF == visit.sequence.site$REEF[1],
                                 DEPTH == visit.sequence.site$DEPTH[1],
                                 SITE_NO == visit.sequence.site$SITE_NO[1],
                                 VISIT_NO == visit.seq[j] |
                                 VISIT_NO == visit.seq[j-1],
                                 GROUP_CODE == code) %>%
                          select(TRANSECT_NO,VISIT_NO,COVER)
                transect.cover.t1 <- transect.cover.pair %>%
                          filter(VISIT_NO == visit.seq[j-1])
                transect.cover.t2 <- transect.cover.pair %>%
                          filter(VISIT_NO == visit.seq[j])

                #compute drop/growth in cover of more than (drop > 0 and growth < 0)
                diff <- mean(transect.cover.t1$COVER) - mean(transect.cover.t2$COVER) 
                visit.sequence.site$DIFF[j] <- diff # store difference

                if (visit.sequence.site$DISTURBANCE[j] %in% 
                       c('d','s','b','c','m','u','f')) { # recorded disturbance
                    rp_id <- rp_id + 1
                    num_dist <- 1
                } else if (diff >= 5.0) { # unrecorded absolute drop >= 5%
                    rp_id <- rp_id + 1
                    num_dist <- num_dist + 1
                    # use U for unknown disturbance that was not recorded 
                    visit.sequence.site$DISTURBANCE[j] <- 'U' 
                } else { # otherwise perform paired t-test to check for small but 
                         # statistically significant drops

                    # catch cases of missing transect observation, ensure we keep
                    # only the pairs that match to ensure the paired t-test is
                    # validly applied
                    if (length(transect.cover.t1$VISIT_NO) < 
                        length(transect.cover.t2$VISIT_NO)) {
                        transect.cover.t2 <- transect.cover.t2 %>%
                               filter(TRANSECT_NO %in% transect.cover.t1$TRANSECT_NO)
                    } else if (length(transect.cover.t1$VISIT_NO) > 
                               length(transect.cover.t2$VISIT_NO)) {
                        transect.cover.t1 <- transect.cover.t1 %>%
                               filter(TRANSECT_NO %in% transect.cover.t2$TRANSECT_NO)
                    }
                    transect.cover.t1 <- transect.cover.t1 %>% arrange(TRANSECT_NO)
                    transect.cover.t2 <- transect.cover.t2 %>% arrange(TRANSECT_NO)
                
                    # perform paired t-test if there are enough transects
                    if (length(transect.cover.t1$TRANSECT_NO) > 1 &&
                        length(transect.cover.t2$TRANSECT_NO) > 1) {
                        # run paired t-test
                        res <- t.test(x = transect.cover.t1$COVER, y = transect.cover.t2$COVER,
                                      alternative = "greater", paired = TRUE, conf.level = alpha)
                        # any statistically significant pair update the RP_ID
                        if (res$p.value < 1.0 - alpha) {
                            rp_id <- rp_id + 1
                            num_dist <- num_dist + 1
                            # use T for small unknown disturbance as per t-test 
                            visit.sequence.site$DISTURBANCE[j] <- 'T'
                        }
                        # store p-value and diff estimate regardless of test result 
                        visit.sequence.site$Pval[j] <- res$p.value
                        visit.sequence.site$DIFF[j] <- res$estimate
                    }
                } 
            }
            visit.sequence.site$RP_ID[j] <- rp_id
            visit.sequence.site$NUM_DIST[j] <- num_dist
        }    
        recovery.trajectories.site <- bind_rows(recovery.trajectories.site,
                                          visit.sequence.site)
    }
    return(recovery.trajectories.site)
}

#' Extract time series by site and transect
#' 
#' @description  Extraction of site level time-series with unique id given for each.
#' @details No splitting of time series, simply combines disturbances, samples 
#' and site data together assigning each \code{REEF/SITE/DEPTH} series a unique ID.
#'
#' @param disturbance        data frame of disturbance records
#' @param samples             data frame of visit time/place records
#' @param transect.cover     data frame of percent cover at transect level, can be at
#'                           major group level or benthic code (ker code) level
#'                          by major group ("GROUP_CODE") or ker code ("BENTHOS_CODE")
#'
#' @family Data extraction
extract_time_series_transect <- function(disturbance,samples,
                                                   transect.cover) {

    # to avoid factor warnings
    combined <- sort(union(union(levels(disturbance$REEF),levels(samples$REEF)),
                           levels(transect.cover$REEF)))
    disturbance <- disturbance %>% 
        mutate(REEF=factor(REEF,levels=combined))
    samples <- samples %>% 
        mutate(REEF=factor(REEF,levels=combined))
    transect.cover <- transect.cover %>% 
        mutate(REEF=factor(REEF,levels=combined))
    
    # join disturbances and samples to transects
    samples.disturbance.transect.cover <- right_join(
                      right_join(disturbance%>%select(-c("P_CODE")),samples,by = c("REEF","VISIT_NO")),
                      select(transect.cover, P_CODE,REEF,DEPTH,VISIT_NO,SITE_NO),
                      by = c("REEF","DEPTH","VISIT_NO")) %>%
                      mutate(RP_ID = 0,DIFF =0,Pval=0)

    print(samples.disturbance.transect.cover)
    recovery.trajectories.site <- data.frame()
    # for every REEF/P_CODE/SITE/DEPTH 4-tuple by visit
    unique_sites <- unique(select(samples.disturbance.transect.cover,
                                   REEF,DEPTH,SITE_NO))
    for (i in 1:nrow(unique_sites)){
        visit.sequence.site <- samples.disturbance.transect.cover %>% 
                                  filter(REEF == unique_sites$REEF[i],
                                         DEPTH == unique_sites$DEPTH[i],
                                         SITE_NO == unique_sites$SITE_NO[i]) %>%
                                  unique() %>%
                                  mutate(RP_ID = i)
        recovery.trajectories.site <- bind_rows(recovery.trajectories.site,
                                          visit.sequence.site)
    }
    return(recovery.trajectories.site)
}

#' Filter recovery trajectories transect
#' 
#' @description  Filter recovery periods from the Austrialian Institute
#'  for Marine Science (AIMS) Long Term Monitoring Program (LTMP) data and Marin 
#'  Monitoring Program (MMP) data.
#'
#' @details  Applies filter conditions to the recovery trajectory data set that
#'  is derived using the \code{extract_recovery_trajectories_transect} function. filter conditions are specified
#'  by a generic user function.
#'
#' @param recovery.trajectories data frame as produced by 
#'                            \code{extract_recovery_trajectories_transect()}
#' @param transect.cover      transect level percent cover data grouped by group codes
#'                             or benthic codes
#' @param derive_func         user prescribed function to append derived quantities
#'                            to trajectories data. This function receives all cover 
#'                            cover data for a single recovery trajectory with RP_ID 
#'                            already appended
#' @param filter_func         user prescribed function to allow a trajectory through
#'                            the filter. This function receives the output from
#'                            \code{derive_func()} as input
#'
#' @return a data frame of trajectories with coral cover data that satisfied the 
#'         filter conditions
#'
#' @examples
#' # Example user functions to filter based on obs > 5 and initial cover < 5%
#' min_obs <- 5
#' max_init <- 5
#' # A derive function that appends the number of observations to each trajectorie
#' append_num_obs <- function(traj) {
#'    num_obs <- length(unique(traj$VISIT_NO))
#'    visits <- traj %>% select(VISIT_NO,Date) %>% unique() %>%  arrange(VISIT_NO)
#'    diff_in_days = as.numeric(visits$Date[num_obs]-visits$Date[1],units = "days")
#'    return(mutate(traj,NUM_OBS = num_obs,DURATION = diff_in_days))
#'}
#' # filter function that thresholds the initial hard coral cover and number of 
#' # observantions
#' trajectory_filter_condition_no_end <- function(traj) {
#'    # Short-cut if not enough observations
#'    if (traj$NUM_OBS[1] < min_obs) {
#'        return(FALSE)
#'    } else {
#'        # condition total coral cover < thresh_cover 
#'        first_visit <- min(traj$VISIT_NO)
#'        init_state <- traj %>% filter(VISIT_NO == first_visit,GROUP_CODE == 'HC')
#'        # site level cover is derived from the mean of transect level cover
#'        tf <- (mean(init_state$COVER) < max_init)
#'        return(tf)
#'    }
#'}
#' @note currently only supports GROUP_CODE grouping, but should implement filter
#' based on BENTHOS_CODE
#' @family Data filtering
filter_recovery_trajectories_transect <- function(recovery.trajectories,
                                    transect.cover,
                                    derive_func = function(x){return(x)},
                                    filter_func = function(x){return(FALSE)}) {

    # create new data set containing only recovery trajectories with the 
    # required conditions
    filt.rec.traj <- data.frame() 
    # to avoid factor warnings/errors
    combined <- sort(union(levels(recovery.trajectories$REEF),
                           levels(transect.cover$REEF)))
    recovery.trajectories <- recovery.trajectories %>% 
                           mutate(REEF=factor(REEF,levels=combined))
    transect.cover <- transect.cover %>% 
                      mutate(REEF=factor(REEF,levels=combined))

    for (rp_id in unique(recovery.trajectories$RP_ID)) {
        
        # extract the trajectory
        visit.sequence <- recovery.trajectories %>%
                          filter(RP_ID == rp_id) %>%
                          arrange(VISIT_NO)
        
        # And coral cover data associated with this trajectory
        cover.dat <- transect.cover %>%
                     filter(REEF == visit.sequence$REEF[1],
                            DEPTH == visit.sequence$DEPTH[1],
                            SITE_NO == visit.sequence$SITE_NO[1],
                            VISIT_NO %in% visit.sequence$VISIT_NO) %>%
                     mutate(RP_ID = rp_id)
        
        ## join to get the dates for each visit to be available for the user
        cover.dat <- right_join(select(visit.sequence,VISIT_NO,Date,DIFF),
                                         cover.dat, by = "VISIT_NO")

        cover.dat <- derive_func(cover.dat)
            
        if (filter_func(cover.dat) == TRUE) {
             filt.rec.traj <- bind_rows(filt.rec.traj,cover.dat)
        }

    }
    filt.rec.traj <- right_join(select(recovery.trajectories,
                                       RP_ID,VISIT_NO,DISTURBANCE),
                                filt.rec.traj, by = c("RP_ID","VISIT_NO") )
    return(filt.rec.traj)
}

#' Reformat recovery trajectories 
#'
#' @description  Reformat a single recovery trajectory as produced by the function
#' \code{filter_recovery_trajectoies_transect} for easier analysis.
#' @details  Rhe transect data is combined but group code to compute the site level
#' data and the standard error included as columns, so that each visit is a single row. 
#' Also the variable T is included as the number of days since the start of recovery.
#'
#' @param traj a single recovery trajectory (\code{RP_ID/SITE_NO} is unique)
#' @param codes string for code field to use
#' @param annotate custom names
#'
#' @return a data table such that result\code{[i,]} is the data for the ith visit in this
#' recovery sequence.
#' @note include option to create new groups from sets of BENTHIC CODES to
#' separate out Achroporites.
#' @family Data reformatting
reformat_recovery_trajectories <- function(traj, codes = 'GROUP_CODE', annotate = "") {

    # extract the constant data and unique visits (date is unique for each visit
    # so the number of rows is the same as the number of visits
    output.seq <- traj %>%
                  select(P_CODE,REEF,SITE_NO,DEPTH,RP_ID,Date,VISIT_NO,
                         DURATION,DISTURBANCE,DIFF) %>%
                  unique() %>%
                  arrange(VISIT_NO) %>%
                  mutate('T' =  as.double(Date - Date[1]))
    groups <- levels(traj[[codes]])
    # append a column for each group and populate
    for (g in groups){
        res <- traj %>% 
                filter(!!rlang::sym(codes) == g) %>%
                group_by(VISIT_NO) %>%
                summarise(mu = mean(COVER),
                          sd = sd(COVER), 
                          se = sd(COVER)/sqrt(length(COVER))) %>%
                arrange(VISIT_NO)
        if(nrow(res) == 0) {
            print(traj)
        } else {
        output.seq <- mutate(output.seq,!!paste(g,annotate,sep="") := res$mu, 
                             !!paste(g,"_sd",annotate,sep="") := res$sd,
                             !!paste(g,"_se",annotate,sep="") := res$se)
        }
    }
    return(output.seq)
}

reformat_recovery_trajectories2 <- function(traj, codes = 'GROUP_CODE', annotate = "") {
  
  # extract the constant data and unique visits (date is unique for each visit
  # so the number of rows is the same as the number of visits
  output.seq <- traj %>%
    select(P_CODE,REEF,SITE_NO,DEPTH,RP_ID,Date,VISIT_NO,
           DURATION,DISTURBANCE,DIFF) %>%
    unique() %>%
    arrange(VISIT_NO) %>%
    mutate('T' =  as.double(Date - Date[1]))
  groups <- levels(traj[[codes]])
  # append a column for each group and populate
  for (g in groups){
    res <- traj %>% 
      filter(!!rlang::sym(codes) == g) %>%
      group_by(VISIT_NO,TRANSECT_NO) %>%
      summarise(mu = mean(COVER),
                sd = sd(COVER), 
                se = sd(COVER)/sqrt(length(COVER))) %>%
      arrange(VISIT_NO)
    if(nrow(res) == 0) {
      print(traj)
    } else {
      output.seq <- mutate(output.seq,!!paste(g,annotate,sep="") := res$mu, 
                           !!paste(g,"_sd",annotate,sep="") := res$sd,
                           !!paste(g,"_se",annotate,sep="") := res$se)
    }
  }
  return(output.seq)
}

#' Create reef group ID
#'
#' @description  Create a reef group ID using kmeans clustering, k = \code{55} for \code{136} unique reefs
#'           REEF_GROUP of zero for NAN lat or lon 
#'           option to export: save(reef_groups, file = \code{"reef_groups.RData"})
#'           future work will include other clustering methods 
#'
#' @param samples data frame (or any df) with REEF, LATITUDE, LONGITUDE 
#' @family Recovery trajectory grouping
create_reef_group_ID <- function(samples){
    ## Initialize all groups as zero (to avoid the reefs with NAN lat/lon)
    samples$REEF_GROUP <- 0
    ## Simple kmeans clustering for reef groups (this can be improved)
    stepone <- cbind(scale(samples$LATITUDE), scale(samples$LONGITUDE))
    steptwo <- na.omit(stepone)
    my_clusters <- kmeans(steptwo, 55)
    
    samples[!is.na(samples$LATITUDE), "REEF_GROUP"] <- my_clusters$cluster
    
    ## Unique row for each ID
    return(as_tibble(select(samples[!duplicated(samples$REEF), ], REEF, REEF_GROUP)))
}

#' Get Group IDs From File
#'
#' @description  Update the Group IDs with clustering in original data from pipeline
#' @section Warning: 
#' This function has some specific elements related to processing the file \code{"groups to murray 4 july 2017.csv"} 
#' to work with the cover data we had been provided with. We suspect AIMS has a more automatic processof doing this and that the \code{manual_reef_key} parameter should not be required. 
#'
#' @details Update the CombinedGroups ID column in the pipeline file with the REEF_GROUP ID from the clustering
#'          done within the \code{create_reef_group_ID()} function
#' 
#' @param manual_reef_key dataframe from the "murraygroups.csv" manual key file that maps different name formats
#' @param data_lonlat any dataframe that has REEF, Longitude and Latitude information
#' @param original_file dataframe of with reef name in different format to our cover data (e.g., the file Angus provided for grouping data from his analysis pipeline. 
#' @family Recovery trajectory grouping
GetGroupIDsFromFile <- function(manual_reef_key, data_lonlat, original_file){
    
    ## Using lon/lat to determine clusters 
    my_REEF_GROUPS <- create_reef_group_ID(data_lonlat) 
    
    ## Add on Angus' REEF naming convention (AngusReef column)
    originalgrouping <- left_join(my_REEF_GROUPS, manual_reef_key, by = "REEF")
    
    ## Now add the clustered REEF_GROUP to Angus' original file (to be put back into the original pipeline)
    left_join(original_file, originalgrouping) %>% 
        mutate(CombinedGroups = REEF_GROUP, REEF = AngusReef) %>% ## Renamed and cleaned 
        select(REEF, CombinedGroups) %>% 
        return()
}
