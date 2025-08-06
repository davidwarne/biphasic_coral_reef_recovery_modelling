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
#
# import the samples and disturbance data
load(paste(DATA_DIR,DIST_DAT_FILE,sep=""))
load(paste(DATA_DIR,SAMPLE_DAT_FILE,sep=""))
load(paste(DATA_DIR,GRP_LEV_TRANS_DAT_FILE,sep=""))
load(paste(DATA_DIR,FAMILY_LEV_TRANS_DAT_FILE,sep=""))

# Extract the transect cover for the family Acroporidae and join as new group
# Therefore ACR = "Acroporidae hard coral" HC - ACR = "other hard corals"
gR <- unique(groups.transect$REEF)
fR <- unique(family.cover.transect$REEF)

temp.acr.hc <- family.cover.transect %>% filter(FAMILY_DESC == "ACROPORIDAE") %>% mutate(GROUP_CODE = "ACR") %>% select(-FAMILY_DESC) %>% bind_rows(groups.transect) 
temp.hc <- temp.acr.hc %>% filter(GROUP_CODE == "HC")
temp.acr <- temp.acr.hc %>% filter(GROUP_CODE == "ACR")
# IMPORTANT! Some inconsistencies in depth recoreds for LTMP Group vs family data
temp.acr$DEPTH[temp.acr$P_CODE != "IN"] <- 6
temp.hc$DEPTH[temp.hc$P_CODE != "IN"] <- 6
# Now the join works (*.x = HC *.y = ACR )
temp.oth <- left_join(temp.hc,temp.acr, by=c("REEF","SITE_NO","TRANSECT_NO","DEPTH","VISIT_NO")) %>% mutate(COVER = COVER.x - COVER.y,GROUP_CODE = "OTH")
# It seems that it is possible to have ACR > HC (from further investigation it seems some SC are classified as ACR)
# this seems to be an exception, so I am 
temp.oth$COVER.y[temp.oth$COVER < 0] <- 0
temp.oth$COVER[temp.oth$COVER < 0] <- 0
temp.acr <- temp.oth %>% 
    select(-c(COVER.x,GROUP_CODE.x,P_CODE.y,COVER,GROUP_CODE)) %>% 
    mutate(COVER = COVER.y, GROUP_CODE = GROUP_CODE.y,P_CODE = P_CODE.x) %>% 
    select(-c(COVER.y,GROUP_CODE.y,P_CODE.x)) 
temp.oth <- temp.oth %>% 
    select(-c(COVER.x,GROUP_CODE.x,P_CODE.y,COVER.y,GROUP_CODE.y)) %>% 
    mutate(P_CODE = P_CODE.x) %>% 
    select(-c(P_CODE.x)) 
acr.cover.transect <- groups.transect %>% 
    bind_rows(temp.acr) %>% 
    bind_rows(temp.oth)
acr.cover.transect$DEPTH[acr.cover.transect$P_CODE != "IN"] <- 6
acr.cover.transect$GROUP_CODE <- as.factor(acr.cover.transect$GROUP_CODE)

# extract the recovery data based on LTMP disturbance records
recovery.trajectories <- extract_recovery_trajectories_transect(disturbance,
                                                                samples,
                                                                acr.cover.transect)

# save group level cover data with extract ACR 
save(acr.cover.transect,file=paste(PROC_DATA_DIR,GRP_ACR_LEV_TRANS_DAT_FILE,sep=""))
# save resulting data table
save(recovery.trajectories, file=paste(PROC_DATA_DIR,REC_TRAJ_DAT_FILE,sep=""))
