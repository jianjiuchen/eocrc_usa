loc_run=T
library(tidyverse)
if(loc_run==T){
  dir_resu=dir_data=dir_code=""
  TargetGender='Male'
  VersionTag="v1"
}else{
  dir_resu=dir_data=dir_code=''
  TargetGender='gen_val'
  VersionTag='ver_val'
}
n_core=if_else(loc_run==T,1,50);core_id=if(loc_run==T) 1 else NTH_CORE
VersionTagNum=as.numeric(str_sub(VersionTag,2,-1))
seer_data=read_csv(str_c(dir_data,"crc_3version.csv")) %>% filter(version=='13') %>% select(-version) %>%
  mutate(race=factor(race,levels=c('White','Black')))
source(str_c(dir_code,"crc_mvm.R"))
load(str_c(dir_data,"p2_rf_list.RData"))#load risk factor data
rf_data=p2_rf_list
FixedPeriod=c(6,7,8,9,10);FixedAge=c(6,7,8,9)
#Period: 6 -> 1992-1996, 7 -> 1997-2001, 8 -> 2002-2006, 9 -> 2007-2011, 10 -> 2012-2016
#Age: 6 -> 30-34, 7 -> 35-39, 8 -> 40-44, 9 -> 45-49
gender=TargetGender#Male/Female
MinNumOfStudyFactor=1
dir_utils=str_c(dir_code,'utils.R')
source(dir_utils)
ref_age=get_ref_grp(VersionTag)[['age']]
ref_race=get_ref_grp(VersionTag)[['race']]
CorCoefThreshold=get_cutoff(VersionTag)
if(VersionTag %in% c('v2')){
  StudyFactorModelCombo=list(c("cs_m.nhis-amc_lag2","cs_m.nhis-amc_cum10y5ylag2","-"),
                             c("adj_avg_alc_m.brfss-amc_lag2","adj_avg_alc_m.brfss-amc_cum10y5ylag2",
                               'drinkhea.m2w2.brfss-amc_lag2','drinkhea.m2w2.brfss-amc_cum10y5ylag2',
                               'drinkhea.m2w1.brfss-amc_lag2','drinkhea.m2w1.brfss-amc_cum10y5ylag2',
                               'drinkbin.brfss-amc_lag2','drinkbin.brfss-amc_cum10y5ylag2',"-"),
                             c("overweight_up_m.nhis-amc_no_lag","overweight_up_m.nhis-amc_lag52",
                               "overweight_up_m.nhis-amc_lag2","overweight_up_m.nhis-amc_cum10y5ylag2",
                               "obese_m.nhis-amc_no_lag","obese_m.nhis-amc_lag52","obese_m.nhis-amc_lag2","obese_m.nhis-amc_cum10y5ylag2","-"),
                             c("hbp_m.nhis-amc_no_lag","hbp_m.nhis-amc_lag52","hbp_m.nhis-amc_lag2","hbp_m.nhis-amc_cum10y5ylag2","-"),
                             c("dia_m.nhis-amc_no_lag","dia_m.nhis-amc_lag52","dia_m.nhis-amc_lag2","dia_m.nhis-amc_cum10y5ylag2","-"),
                             c("thy_m.fill-amc_no_lag","thy_m.fill-amc_lag52","thy_m.fill-amc_cum10y5ylag2","-"),
                             c("calcium_m.fill-amc_cum10y5ylag2","-"),c("fibe_m.fill-amc_cum10y5ylag2","-"),c("wholegrain_m.nhanes-amc_cum10y5ylag2","-"),
                             c("dairy_m.nhanes-amc_cum10y5ylag2","-"),c("rednproc_m.nhanes-amc_cum10y5ylag2","-"))
}else if(VersionTag %in% c('v3')){
  StudyFactorModelCombo=list(c("cs_m.brfss-amc_lag2","cs_m.brfss-amc_cum10y5ylag2","-"),
                             c("adj_avg_alc_m.brfss-amc_lag2","adj_avg_alc_m.brfss-amc_cum10y5ylag2",
                               'drinkhea.m2w2.brfss-amc_lag2','drinkhea.m2w2.brfss-amc_cum10y5ylag2',
                               'drinkhea.m2w1.brfss-amc_lag2','drinkhea.m2w1.brfss-amc_cum10y5ylag2',
                               'drinkbin.brfss-amc_lag2','drinkbin.brfss-amc_cum10y5ylag2',"-"),
                             c("overweight_up_m-amc_no_lag","overweight_up_m-amc_lag52","overweight_up_m.fill-amc_cum10y5ylag2",
                               "obese_m-amc_no_lag","obese_m-amc_lag52","obese_m.fill-amc_cum10y5ylag2","-"),
                             c("hbp_m.nhanes.fill-amc_no_lag","hbp_m.nhanes.fill-amc_lag52","hbp_m.nhanes.fill-amc_cum10y5ylag2","-"),
                             c("dia_m.nhanes.fill-amc_no_lag","dia_m.nhanes.fill-amc_lag52","dia_m.nhanes.fill-amc_cum10y5ylag2","-"),
                             c("thy_m.fill-amc_no_lag","thy_m.fill-amc_lag52","thy_m.fill-amc_cum10y5ylag2","-"),
                             c("calcium_m.fill-amc_cum10y5ylag2","-"),c("fibe_m.fill-amc_cum10y5ylag2","-"),c("wholegrain_m.nhanes-amc_cum10y5ylag2","-"),
                             c("dairy_m.nhanes-amc_cum10y5ylag2","-"),c("rednproc_m.nhanes-amc_cum10y5ylag2","-"))
}else if(VersionTag %in% c('v1')){
  StudyFactorModelCombo=list(c("cs_m.brfss-amc_lag2","cs_m.brfss-amc_cum10y5ylag2","-"),
                             c("adj_avg_alc_m.brfss-amc_lag2","adj_avg_alc_m.brfss-amc_cum10y5ylag2",
                               'drinkhea.m2w2.brfss-amc_lag2','drinkhea.m2w2.brfss-amc_cum10y5ylag2',
                               'drinkhea.m2w1.brfss-amc_lag2','drinkhea.m2w1.brfss-amc_cum10y5ylag2',
                               'drinkbin.brfss-amc_lag2','drinkbin.brfss-amc_cum10y5ylag2',"-"),
                             c("overweight_up_m.brfss.fill-amc_no_lag","overweight_up_m.brfss.fill-amc_lag52","overweight_up_m.brfss.fill-amc_cum10y5ylag2",
                               "obese_m.brfss.fill-amc_no_lag","obese_m.brfss.fill-amc_lag52","obese_m.brfss.fill-amc_cum10y5ylag2","-"),
                             c("hbp_m_no_preg.brfss-amc_no_lag","hbp_m_no_preg.brfss-amc_lag52","hbp_m_no_preg.brfss-amc_lag2",
                               "hbp_m_no_preg.brfss-amc_cum10y5ylag2","-"),
                             c("dia_m.brfss.fill-amc_no_lag","dia_m.brfss.fill-amc_lag52","dia_m.brfss.fill-amc_cum10y5ylag2","-"),
                             c("thy_m.fill-amc_no_lag","thy_m.fill-amc_lag52","thy_m.fill-amc_cum10y5ylag2","-"),
                             c("calcium_m.fill-amc_cum10y5ylag2","-"),c("fibe_m.fill-amc_cum10y5ylag2","-"),c("wholegrain_m.nhanes-amc_cum10y5ylag2","-"),
                             c("dairy_m.nhanes-amc_cum10y5ylag2","-"),c("rednproc_m.nhanes-amc_cum10y5ylag2","-"))
}

tm_before=Sys.time()

OutCombo=crc_mvm(rf_data=rf_data,seer_data=seer_data,gender=gender,
                 StudyFactorModelCombo=StudyFactorModelCombo,
                 FixedPeriod=FixedPeriod,FixedAge=FixedAge,ref_age=ref_age,ref_race=ref_race,
                 CorCoefThreshold=CorCoefThreshold,MinNumOfStudyFactor=MinNumOfStudyFactor,
                 n_core=n_core,core_id=core_id,
                 dir_utils=dir_utils)


if(n_core>1){
  eval(parse(text=str_c("OutCombo_",VersionTag,"_",gender,"_core",core_id,"=OutCombo")))
  eval(parse(text=str_c("save(OutCombo_",VersionTag,"_",gender,"_core",core_id,",file='",dir_resu,"OutComb_",VersionTag,"_",gender,"_core",core_id,".RData')")))
  
  #get the file name for the output of all other parallel cores
  #if all are ready, combine & save
  exist_filename=str_c('OutComb_',VersionTag,'_',gender,'_core') %>% str_c(.,1:n_core) %>% str_c(.,'.RData')
  if(all(exist_filename %in% list.files(dir_resu))){
    for(file_i in exist_filename) load(str_c(dir_resu,file_i))
    exist_OutCombo=str_c('OutCombo_',VersionTag,'_',gender,'_core') %>% str_c(.,1:n_core)
    OutCombo=NULL
    for (OutCombo_i in exist_OutCombo)  eval(parse(text=str_c('OutCombo=bind_rows(OutCombo,',OutCombo_i,')')))
    
    eval(parse(text=str_c("OutCombo_",VersionTag,"_",gender,"=OutCombo")))
    eval(parse(text=str_c("save(OutCombo_",VersionTag,"_",gender,",file='",dir_resu,"OutComb_",VersionTag,"_",gender,".RData')")))
  }
}else{
  eval(parse(text=str_c("OutCombo_",VersionTag,"_",gender,"=OutCombo")))
  eval(parse(text=str_c("save(OutCombo_",VersionTag,"_",gender,",file='",dir_resu,"OutComb_",VersionTag,"_",gender,".RData')")))
}

print(paste0(VersionTag," for ", gender," took ",round(as.numeric(Sys.time()-tm_before,units="hours"),2), " hours."))



