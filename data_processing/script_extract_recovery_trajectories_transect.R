##
# @file script_extract_recovery_trajectories.R
# @brief script to extract out recovery periods from the Austrialian Institute
#  for Marine Science (AIMS) Long Term Monitoring Program (LTMP) data and Marin 
#  Monitoring Program (MMP) data.
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



# import the samples and disturbance data
load(paste(DATA_DIR,DIST_DAT_FILE,sep=""))
load(paste(DATA_DIR,SAMPLE_DAT_FILE,sep=""))
load(paste(DATA_DIR,GRP_LEV_TRANS_DAT_FILE,sep=""))


# extract the recovery data based on LTMP disturbance records
recovery.trajectories <- extract_recovery_trajectories_transect(disturbance,
                                                                samples,
                                                                groups.transect)

# save resulting data table
save(recovery.trajectories, file=paste(PROC_DATA_DIR,REC_TRAJ_DAT_FILE,sep=""))
