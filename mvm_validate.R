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
##########
#prepare
##########
source(str_c(dir_code,'crc_mvm_other.R'))
source(str_c(dir_code,'utils.R'))
dir_utils=str_c(dir_code,'utils.R')


##########
#get AllModelsTmp
##########
print(Sys.time())
AllModelsTmp=find_nest_mod(proc_mvm_interact_res(VersionTag,TargetGender,dir_data=dir_data,dir_utils=dir_utils))
parsim_mod_ls_attr=attributes(AllModelsTmp)$parsim_mod_ls
gender_attr=attributes(AllModelsTmp)$gender
RankBy_attr=attributes(AllModelsTmp)$RankBy
cutoff_attr=attributes(AllModelsTmp)$cutoff
RfSource_attr=attributes(AllModelsTmp)$RfSource
vTag_attr=attributes(AllModelsTmp)$vTag
ref_grp_attr=attributes(AllModelsTmp)$ref_grp
#
clean_dum_var_nam=function(cov_nam_raw){
  cov_nam_raw %>% str_replace(.,'rBlack1','rBlack') %>% str_replace(.,'rWhite1','rWhite') %>% 
    str_replace(.,'a61','a6') %>% str_replace(.,'a71','a7') %>% str_replace(.,'a81','a8') %>% str_replace(.,'a91','a9')
}
AllModelsTmp=AllModelsTmp %>% mutate(cov_nam_raw=clean_dum_var_nam(cov_nam_raw))
print(Sys.time())

##########
#based on cov_nam_raw get rf data ready (apply lag, amc, etc)
##########
ParsimModels=AllModelsTmp %>% filter(parsim_mod)
all_study_factor_model_nameS=unique(ParsimModels$cov_nam_raw)
all_study_factor_model_nameS=all_study_factor_model_nameS[str_detect(all_study_factor_model_nameS,pattern='\\:',negate=T)]
AllRF=as.data.frame(str_split(all_study_factor_model_nameS,'_amc'))[1,] %>% unlist %>% {x=.;x=as.character(x);names(x)=NULL;x} %>% unique
load(str_c(dir_data,"p2_rf_list.RData"))
col_nam_PerAgeGenRac=c("period","age_grp","gender","race")
col_nam_PerAgeRac=c("period","age_grp","race")
std_fun=function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)
AllRfData=p2_rf_list %>% select(str_c(col_nam_PerAgeGenRac),str_c(AllRF)) %>% filter(gender==TargetGender) %>% select(-gender)
all_study_factor_nameS=as.data.frame(str_split(all_study_factor_model_nameS,"_amc"))[1,] %>% unlist %>% {x=.;x=as.character(x);names(x)=NULL;x}#may have duplicates
#all_study_factor_nameS doesn't contain duplicates, but two elements in it may have the same study_factor or the same model
all_model_nameS=as.data.frame(str_split(all_study_factor_model_nameS,"_amc"))[2,] %>% unlist %>% {x=.;x=as.character(x);names(x)=NULL;x}#may have duplicates
all_model_nameS=str_c('amc',all_model_nameS)
#note: unlist rtrn char on Mac and fctc on cluster; so coerce to char
#the above two var (all_study_factor_nameS & all_model_nameS) are meant to contain duplicates, eg,
#adj_avg_alc_m.brfss, adj_avg_alc_m.brfss, fibe_m.fill
#amc_lag2, amc_cum10y5ylag2, amc_cum10y5ylag2
TmpData=as.list(1:length(all_study_factor_model_nameS))
for (i in 1:length(all_study_factor_model_nameS)) {
  Fn_tmp=get(str_c("fun_",all_model_nameS[i]))
  #TmpData[[i]]=Fn_tmp(AllRfData,all_study_factor_nameS[i],PeriodForOutput=PeriodForOutput,AgeForOutput=AgeForOutput) %>%
  #  rename(all_study_factor_model_nameS[i]=all_study_factor_nameS[i])
  eval(parse(text=str_c("TmpData[[i]]=Fn_tmp(AllRfData,all_study_factor_nameS[i],PeriodForOutput=c(6,7,8,9,10),AgeForOutput=c(6,7,8,9)) %>%
                          rename(",all_study_factor_model_nameS[i],"=",all_study_factor_nameS[i],")")))
  #Fn_tmp requires col_nam_PerAgeRac & std_fun in the outer env; I did it.
  if(i==1){
    DataAllStudyFactorModel=TmpData[[i]]
  } else{
    DataAllStudyFactorModel=suppressMessages(inner_join(DataAllStudyFactorModel,TmpData[[i]]))
  }
}
rm(TmpData)
DataAllStudyFactorModel=DataAllStudyFactorModel %>% drop_na
print(Sys.time())
print('Finished getting exposure data!')
##########
#get crc data ready
#combine exposure and crc data & dummy code
##########
seer_data=suppressMessages(read_csv(str_c(dir_data,"crc_3version.csv"))) %>% 
  mutate(version=if_else(version=="9_13_18","13_18",version),race=factor(race,levels=c('White','Black')))
seer_data=seer_data %>% filter(version=="13",gender==TargetGender) %>% select(-gender) %>% as_tibble
DataForReg=suppressMessages(inner_join(seer_data,DataAllStudyFactorModel))
#dummy coding for age, period, and race
DataForReg=dummy_code_AgePerRac(DataForReg,ref_age=get_ref_grp(VersionTag)[['age']],ref_race=get_ref_grp(VersionTag)[['race']])
cv_fold=5
period_ls=as.list(numeric(cv_fold))
for (i in 1:length(period_ls)) {
  period_ls[[i]]=combn(6:10,ifelse(cv_fold==10,2,1))[,i]
}
PopSet=TestSet=TrainSet=row4period=as.list(numeric(length(period_ls)))
for (i in 1:length(period_ls)) {
  TrainSet[[i]]=DataForReg %>% filter(!(period %in% period_ls[[i]]))
  TestSet[[i]]=DataForReg %>% filter(period %in% period_ls[[i]]) %>% .$count
  PopSet[[i]]=DataForReg %>% filter(period %in% period_ls[[i]]) %>% .$pop
  row4period[[i]]=DataForReg %>% mutate(id=period %in% period_ls[[i]]) %>% .$id#the rows corresponding to the test set
}
AllDemoTerm=str_c(attributes(DataForReg)$AllAgeTerm,attributes(DataForReg)$AllRacTerm,sep='+')#for later use; save time
##########
#loop thr all parsim models, give a out-of-sample fit measure
#get a df: 1col being model_rank & another col being out-of-sample fit; 1row for each model_rank
##########
cal_out_logl=function(period_ls,DataForReg,covars,TrainSet,PopSet,TestSet,row4period){
  #ARGUMENTS:
  #period_ls: [[1]] [1] 6 7 [[2]] [1] 6 8 ...; len=10
  #DataForReg: contain crc and rf data (afte lag & amc, etc), and after dummy coding for age, period, race
  #covars: example: 'a7+a8+a9+rBlack+drinkbin.brfss_amc_cum10y5ylag2+hbp_m.nhanes.fill_amc_lag52+overweight_up_m.nhis_amc_lag2+
  #cont:             rBlack:drinkbin.brfss_amc_cum10y5ylag2+rBlack:overweight_up_m.nhis_amc_lag2'
  #TrainSet: a list of 10 (or 5) training sets
  #PopSet: a list of 10 (or 5) vectors
  #TestSet: a list of 10 (or 5) vectors
  #row4period: a list of 10 (or 5) boolean vectors; each indicating the rows in the 40-row df that correspond to the two chosen periods
  #NOTE/VALUE:
  #calculate out-of-sample model performance given a model (covars determined by model_rank)
  
  #fit the fullset, ie, DataForReg
  suppressWarnings(eval(parse(text=str_c('RegResFull=MASS::glm.nb(count~',covars,'+offset(log(pop)),data=DataForReg)'))))
  #get design matrix for the 10 test sets; store in a list
  DesiMatrFull=model.matrix(RegResFull)
  DesiMatr=as.list(numeric(length(period_ls)))
  for (i in 1:length(period_ls)){
    DesiMatr[[i]]=DesiMatrFull[row4period[[i]],]
  }
  #main
  out_logl_vec=numeric(length(period_ls))
  RegRes=as.list(numeric(length(period_ls)))
  for(i in 1:length(period_ls)){
    #RegRes[[i]]=MASS::glm.nb(formula,data=TrainSet[[i]])
    suppressWarnings(eval(parse(text=str_c('RegRes[[i]]=MASS::glm.nb(count~',covars,'+offset(log(pop)),data=TrainSet[[i]])'))))
    out_logl_vec[i]=sum(dnbinom(x=TestSet[[i]],
                                mu=exp(DesiMatr[[i]] %*% matrix(coef(RegRes[[i]]),ncol=1))*PopSet[[i]],
                                size=RegRes[[i]]$theta,log=T))
  }
  sum(out_logl_vec)
}


df_out_logl=ParsimModels %>% select(model_rank) %>% distinct %>% mutate(out_logl=NA_real_)
for (mod_rank_i in df_out_logl$model_rank) {
  covars=str_c(AllDemoTerm,
               AllModelsTmp %>% filter(model_rank==mod_rank_i) %>% .$cov_nam_raw %>% str_c(.,collapse='+'),
               sep='+')
  df_out_logl[df_out_logl$model_rank==mod_rank_i,"out_logl"]=cal_out_logl(period_ls,DataForReg,covars,TrainSet,PopSet,TestSet,row4period)
}
AllModels=left_join(AllModelsTmp,df_out_logl)

##########
#I need to visualize goodness-of-fit and cross-validation performance for some models, especially top ones (by BIC or CV)
#I prepare such data here
#build the function below for this end
##########
get_visual_GOFnCV=function(model_rank_i){
  #VALUES
  #list(Data=df,FullModPred=FullModPred,CVPred=CVPred); these are needed for plotting
  #Data has period, age_grp, race, count, pop, rf in the top model
  #FullModPred contains prediction of the model trained using full data
  #CVPred contains prediction of the model trained using part of data; CV_i col in CVPred tells the fold
  #
  #ARGUMENT
  #model_rank_i is a specific value in model_rank col of df_out_logl or of AllModels
  #
  #NOTES
  #this fun uses many objs in Glob Env, eg, AllModelsTmp
  cov_nam_raw_tmp=AllModelsTmp %>% filter(model_rank==model_rank_i) %>% .$cov_nam_raw
  cov_nam_raw_main_eff_tmp=cov_nam_raw_tmp[str_detect(cov_nam_raw_tmp,patter='\\:',negate=T)]
  df=DataForReg %>% select(c('period','age_grp','race','count','pop'),str_c(cov_nam_raw_main_eff_tmp))
  covars_tmp=str_c(AllDemoTerm,
                   cov_nam_raw_tmp %>% str_c(.,collapse='+'),
                   sep='+')
  suppressWarnings(eval(parse(text=str_c('RegResFull=MASS::glm.nb(count~',covars_tmp,'+offset(log(pop)),data=DataForReg)'))))
  DesiMatrFull=model.matrix(RegResFull)
  FullModPred=(exp(DesiMatrFull %*% matrix(coef(RegResFull),ncol=1))*df$pop) %>% as.vector %>% 
    tibble(FullMod_mu=.,FullMod_size=RegResFull$theta) %>% bind_cols(df[,c('period','age_grp','race')],.)
  DesiMatr=as.list(numeric(length(period_ls)))
  for (i in 1:length(period_ls)){
    DesiMatr[[i]]=DesiMatrFull[row4period[[i]],]
  }
  RegRes=as.list(numeric(length(period_ls)))
  CVPred=NULL
  for(i in 1:length(period_ls)){
    #RegRes[[i]]=MASS::glm.nb(formula,data=TrainSet[[i]])
    suppressWarnings(eval(parse(text=str_c('RegRes[[i]]=MASS::glm.nb(count~',covars_tmp,'+offset(log(pop)),data=TrainSet[[i]])'))))
    bind_cols(as_tibble(df[row4period[[i]],c('period','age_grp','race')]),
              CV_mu=(exp(DesiMatr[[i]] %*% matrix(coef(RegRes[[i]]),ncol=1))*PopSet[[i]]) %>% as.vector) %>% 
      mutate(CV_size=RegRes[[i]]$theta,CV_i=i) %>%
      bind_rows(CVPred,.) -> CVPred
  }
  CVPred=CVPred %>% arrange(period,age_grp,race)
  list(Data=df,FullModPred=FullModPred,CVPred=CVPred)
}


#use get_visual_GOFnCV
top_model_rank_i_BIC=AllModels %>% filter(parsim_mod) %>% .$model_rank %>% unique %>% .[1:10]
top_model_rank_i_CV=AllModels %>% arrange(desc(out_logl)) %>% filter(parsim_mod) %>% .$model_rank %>% unique %>% .[1:10]
visual_GOFnCV_BIC=as.list(1:length(top_model_rank_i_BIC))
visual_GOFnCV_CV=as.list(1:length(top_model_rank_i_CV))
for (i in top_model_rank_i_BIC)  visual_GOFnCV_BIC[[which(top_model_rank_i_BIC==i)]]=get_visual_GOFnCV(i)
for (i in top_model_rank_i_CV)  visual_GOFnCV_CV[[which(top_model_rank_i_CV==i)]]=get_visual_GOFnCV(i)


##########
#output
##########
attr(AllModels,'parsim_mod_ls')=parsim_mod_ls_attr
attr(AllModels,"gender")=gender_attr
attr(AllModels,'RankBy')=RankBy_attr
attr(AllModels,'cutoff')=cutoff_attr
attr(AllModels,'RfSource')=RfSource_attr
attr(AllModels,'vTag')=vTag_attr
attr(AllModels,'ref_grp')=ref_grp_attr
attr(AllModels,'visual_GOFnCV_BIC')=visual_GOFnCV_BIC
attr(AllModels,'visual_GOFnCV_CV')=visual_GOFnCV_CV

save(AllModels,file=str_c(dir_resu,'OutComb_',VersionTag,'_',TargetGender,'_interact_',cv_fold,'FoldCV.RData'))
print(str_c('Finished ',TargetGender,' and ',VersionTag))






