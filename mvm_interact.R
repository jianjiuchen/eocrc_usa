loc_run=T
library(tidyverse)
if(loc_run==T){
  dir_resu=dir_data=dir_code=""
  TargetGender='Male'
  VersionTag="v1b"
}else{
  dir_resu=dir_data=dir_code=''
  TargetGender='gen_val'
  VersionTag='ver_val'
}
n_core=if_else(loc_run==T,1,50);core_id=if(loc_run==T) 1 else NTH_CORE
crc_mvm_tst_interact_version="v1"
gender=TargetGender
seer_data=read_csv(str_c(dir_data,"crc_3version.csv")) %>% filter(version=='13') %>% select(-version) %>%
  mutate(race=factor(race,levels=c('White','Black')))
load(str_c(dir_data,"p2_rf_list.RData"))
rf_data=p2_rf_list
source(str_c(dir_code,'crc_mvm_other.R'))
source(str_c(dir_code,'utils.R'))
dir_utils=str_c(dir_code,'utils.R')
crc_mvm_tst_interact=get(str_c("crc_mvm_tst_interact_",crc_mvm_tst_interact_version))

tm_before=Sys.time()
AllModels=find_nest_mod(proc_mvm_res(VersionTag,TargetGender,Var2Drop=get_Var2Drop(VersionTag,TargetGender),dir_data=dir_data,dir_utils=dir_utils))
OutCombo=crc_mvm_tst_interact(AllModels,rf_data=rf_data,seer_data=seer_data,
                              n_core=n_core,core_id=core_id,dir_utils=dir_utils)


if(n_core>1){
  eval(parse(text=str_c("OutCombo_",VersionTag,"_",gender,"_interact_core",core_id,"=OutCombo")))
  eval(parse(text=str_c("save(OutCombo_",VersionTag,"_",gender,"_interact_core",core_id,
                        ",file='",dir_resu,"OutComb_",VersionTag,"_",gender,"_interact_core",core_id,".RData')")))
  
  #get the file name for the output of all other parallel cores
  #if all are ready, combine & save
  exist_filename=str_c('OutComb_',VersionTag,'_',gender,'_interact_core') %>% str_c(.,1:n_core) %>% str_c(.,'.RData')
  if(all(exist_filename %in% list.files(dir_resu))){
    print(str_c('===================================FINAL FILE HERE!========================================='))
    for(file_i in exist_filename) load(str_c(dir_resu,file_i))
    exist_OutCombo=str_c('OutCombo_',VersionTag,'_',gender,'_interact_core') %>% str_c(.,1:n_core)
    OutCombo=NULL
    for (OutCombo_i in exist_OutCombo)  eval(parse(text=str_c('OutCombo=bind_rows(OutCombo,',OutCombo_i,')')))
    
    eval(parse(text=str_c("OutCombo_",VersionTag,"_",gender,"_interact=OutCombo")))
    eval(parse(text=str_c("save(OutCombo_",VersionTag,"_",gender,"_interact",
                          ",file='",dir_resu,"OutComb_",VersionTag,"_",gender,"_interact.RData')")))
  }
}else{
  eval(parse(text=str_c("OutCombo_",VersionTag,"_",gender,"_interact=OutCombo")))
  eval(parse(text=str_c("save(OutCombo_",VersionTag,"_",gender,"_interact",
                        ",file='",dir_resu,"OutComb_",VersionTag,"_",gender,"_interact.RData')")))
}

print(paste0(VersionTag," for ", gender," took ",round(as.numeric(Sys.time()-tm_before,units="hours"),2), " hours."))
