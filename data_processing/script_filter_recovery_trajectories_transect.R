
##
# @file script_filter_recovery_trajectories.R
# @brief script to filter recovery periods from the Austrialian Institute
#  for Marine Science (AIMS) Long Term Monitoring Program (LTMP) data and Marin 
#  Monitoring Program (MMP) data.
#
# @details applies filter conditions to the recovery trajectory data set that
#  is derived using the extract_recoveries.R script. This script is intended to
#  be exporative, so it does not provide generic filter capability, just particular
#  features of interest.
#
# @author David J. Warne (david.warne@qut.edu.au)
#           ARC Centre of Excellence for Mathematical and Statistical Frontiers
#           School of Mathematical Sciences 
#           Science and Engineering Faculty
#           Queensland University of Technology
#
# @author Grace E. M. Heron (g.heron@qut.edu.au)
#           ARC Centre of Excellence for Mathematical and Statistical Frontiers
#           School of Mathematical Sciences 
#           Science and Engineering Faculty
#           Queensland University of Technology
#
#

# User functions for filtering
#-------------------------------------------------------------------------------

##
# @brief append derived data function
# @detail This function just adds the number of observations
append_num_obs <- function(traj) {
    num_obs <- length(unique(traj$VISIT_NO))
    visits <- traj %>% 
                select(VISIT_NO,Date) %>% 
                unique() %>% 
                arrange(VISIT_NO)
    diff_in_days = as.numeric(visits$Date[num_obs]-visits$Date[1],units = "days")
    return(mutate(traj,NUM_OBS = num_obs,DURATION = diff_in_days))
}


##
# @brief Trajectory filter function
# @detail User defined filter function implementing the condition
#   HC_0 < 5 % and HC_tn > 15% and n >= 5
#
# @param traj coral cover data subseddt that matches a rp_id/site pair
trajectory_filter_condition <- function(traj) {
    # Short-cut if not enough observations
    if (traj$NUM_OBS[1] < min_obs) {
        return(FALSE)
    } else {
        # condition total coral cover % = total HC % + total SC % < thresh_cover %
        first_visit <- min(traj$VISIT_NO)
        last_visit <- max(traj$VISIT_NO)
        init_state <- traj %>% 
                      filter(VISIT_NO == first_visit,GROUP_CODE == 'HC')
        final_state <- traj %>% 
                       filter(VISIT_NO == last_visit,GROUP_CODE == 'HC')
        # site level cover is derived from the mean of transect level cover
        return(mean(init_state$COVER) < max_init 
                && mean(final_state$COVER) > min_final)
    }
}

##
# @brief Trajectory filter function
# @detail User defined filter function implementing the condition
#   HC_0 < 5 % and n >= 5
#
# @param traj coral cover data subseddt that matches a rp_id/site pair
trajectory_filter_condition_no_end <- function(traj) {
    # Short-cut if not enough observations
    if (traj$NUM_OBS[1] < min_obs) {
        return(FALSE)
    } else {
        # condition total coral cover % = total HC % + total SC % < thresh_cover %
        first_visit <- min(traj$VISIT_NO)
        last_visit <- max(traj$VISIT_NO)
        init_state <- traj %>% 
                      filter(VISIT_NO == first_visit,GROUP_CODE == 'HC')
        final_state <- traj %>% 
                       filter(VISIT_NO == last_visit,GROUP_CODE == 'HC')
        # site level cover is derived from the mean of transect level cover
        tf <- (mean(init_state$COVER) < max_init 
               && mean(final_state$COVER) > mean(init_state$COVER))
        # check no complete zeros
        for (visit in traj$VISIT_NO) {
            state <- traj %>% 
                      filter(VISIT_NO == visit,GROUP_CODE == 'HC')
            tf <- tf && (mean(state$COVER) > 0.0)
        }
        return(tf)
    }
}


# process filter
#-------------------------------------------------------------------------------

# load recovery trajectory data, cover data, and spatial data
load(paste(PROC_DATA_DIR,REC_TRAJ_DAT_FILE,sep=""))
load(paste(DATA_DIR,GRP_LEV_TRANS_DAT_FILE,sep=""))


# perform filtering 
filt.rec.traj <- filter_recovery_trajectories_transect(recovery.trajectories,
                                     groups.transect,
                                     derive_func = append_num_obs,
                                     filter_func = trajectory_filter_condition_no_end)
# save the trajectory data
save(filt.rec.traj,file=paste(PROC_DATA_DIR,
                        sprintf(FILTER_OUT_DATA_FMT,max_init,min_final,min_obs)
                        ,".RData",sep=""))
