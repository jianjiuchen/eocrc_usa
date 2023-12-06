##################
#data processing
##################
proc_mvm_res=function(vTag,gender,RankBy='BIC',Var2Drop=NULL,AMC_only=T,
                      dir_data,
                      dir_utils){
  #
  #ARGUMENTS
  #vTag & gender (Male/Female): needed to import data
  #RankBy='AIC'/'BIC'
  #AMC_only means only keeping AMC models, ie, dropping all quintile or non-amc cont models
  #
  #NOTE
  #can be run on cluster and Mac
  
  #########
  #prepare
  #########
  require(tidyverse)
  if(str_sub(vTag,-1,-1) %in% letters) vTagFile=str_sub(vTag,1,-2) else vTagFile=vTag # if contain letter, we process results with interaction terms in models
  load(str_c(dir_data,'OutComb_',vTagFile,'_',gender,'.RData'))
  eval(parse(text=str_c('data=OutCombo_',vTagFile,'_',gender)))
  rm(list=ls(pattern='^OutCombo'))
  source(dir_utils)
  gender=if(gender=='Male') 'Men' else if(gender=='Female') 'Women'
  cutoff=get_cutoff(vTag)
  cutoff_dbl=as.numeric(cutoff)
  RfSource=get_RfSource(vTag)
  ref_grp=get_ref_grp(vTag)
  
  ########
  #subsetting cutoff, Var2Drop, & rank by GOF, drop high VIF, etc
  ########
  data=data %>% filter(max_abs_cor_coef<=cutoff_dbl)
  n1=length(unique(data$model))
  #
  data=data %>% {
    #drop models that contain certain study factor, eg, overw or obes
    x=.
    if(is.null(Var2Drop)){
      x
    }else{
      x %>% group_by(model) %>% mutate(AnyVar2Drop=any(Var2Drop %in% covariate_name)) %>% ungroup %>% filter(!AnyVar2Drop) %>% select(-AnyVar2Drop)
    }
  } %>% {
    #dropping all quintile or non-amc cont models
    x=.
    if(AMC_only){
      x %>% mutate(amcID=str_detect(covariate_name,'_amc_')) %>% group_by(model) %>% mutate(amcID_all=all(amcID)) %>% ungroup %>% 
        filter(amcID_all) %>% dplyr::select(-amcID,-amcID_all)
    }else x
  } %>% {
    #rank by aic or bic
    x=.
    if(RankBy=='AIC') x %>% arrange(AIC,BIC,covariate_name) else if(RankBy=='BIC') x %>% arrange(BIC,AIC,covariate_name)
  } %>% select(model,everything()) %>% 
    #the line below deals with this: some models are too saturated, and some IRRs are NAs; (actually, this doesn't happen until cutoff>=0.7)
    group_by(model) %>% mutate(WithNA=if(any(is.na(IRR))) T else F) %>% ungroup %>% filter(WithNA==F) %>% select(-WithNA)
  n2=length(unique(data$model))
  
  ########
  #create variables
  ########
  AllModels=data %>% mutate(model_rank=as.numeric(factor(model,levels=unique(model)))) %>%
    mutate(cov_nam_with_linebreak=recode_cov_nam_2_fct(covariate_name,inclu_line_break=T),
           cov_nam_no_linebreak=recode_cov_nam_2_fct(covariate_name,inclu_line_break=F)) %>%
    group_by(model_rank) %>% mutate(eq=str_c('EOCRC ~ ',str_c(cov_nam_no_linebreak,collapse=' + '))) %>% ungroup %>%
    mutate(IRRstr=turn_irr_to_str(IRR),IRR_low_str=turn_irr_to_str(IRR_low),IRR_up_str=turn_irr_to_str(IRR_up)) %>% {
      #rank by aic or bic
      x=.
      if(RankBy=='AIC'){
        x %>% mutate(P=if_else(`p-value`<0.05,"P"," "),B=if_else(delta_BIC>0,"B"," "),A=if_else(delta_AIC>0,"*"," "))
      }else if(RankBy=='BIC'){
        x %>% mutate(P=if_else(`p-value`<0.05,"P"," "),A=if_else(delta_AIC>0,"A"," "),
                     B=case_when(delta_BIC>2&delta_BIC<=6 ~ '*',delta_BIC>6&delta_BIC<=10 ~ '**',delta_BIC>10 ~ '***',T ~ ''))
      }
    } %>%
    select(model_rank,cov_nam_with_linebreak,IRR,IRR_low,IRR_up,IRRstr,IRR_low_str,IRR_up_str,
           pvalue=`p-value`,dBIC=delta_BIC,dAIC=delta_AIC,P,B,A,BIC,AIC,cov_nam_raw=covariate_name,cov_nam_no_linebreak,eq)
  
  ########
  #attr
  ########
  attr(AllModels,"gender")=gender
  attr(AllModels,'RankBy')=RankBy
  attr(AllModels,'cutoff')=cutoff
  attr(AllModels,'RfSource')=RfSource
  attr(AllModels,'vTag')=vTag
  attr(AllModels,'ref_grp')=ref_grp
  attr(AllModels,'Flow of model counts')=str_c(n1,' models initially;\n',n1-n2,' dropped due to Var2Drop or AMC_only or IRR=NA;\n',n2,' models remained.')
  
  print(str_c('For ',vTag,' ',gender,' RfSource_',str_c(RfSource,collapse='&'),' cutoff',cutoff))
  print(str_c('The reference group used were ',ref_grp[['age']],' for age & ',ref_grp[['race']],' for race'))
  print(str_c('Dropped ',str_c(Var2Drop,collapse=' & ')))
  print(attr(AllModels,'Flow of model counts'))
  
  AllModels
}

get_mod4plot_est=function(df,n_mod2show=10){
  #VALUE
  #a df holding the top models; ensure each model has the same set of cov_nam, some including NAs
  #need to plot est, eg, IRRs
  #
  #ARGUMENT
  #df is returned by proc_mvm_res and then find_nest_mod
  #df can also be returned by proc_mvm_interact_res and find_nest_mod
  #df can also be a product of mvm_validate.R (this script only adds a col to the df above)
  #n_mod2show: number of models to show; top 5 or 10? default to 10
  #
  #NOTE
  #here, I excluded models nested in higher ranking models
  #here the output has no attr
  
  #get models that are (1) parsimonious AND (2) top
  #join with NAtibble
  Mod2Show=df %>% filter(parsim_mod)
  RankOfLastMod2Show=unique(Mod2Show$model_rank)[n_mod2show]
  RankOfMod2Show=unique(Mod2Show$model_rank)[1:n_mod2show]
  Mod2Show=Mod2Show %>% filter(model_rank %in% RankOfMod2Show) #CONDITION 'model_rank<=RankOfLastMod2Show' is wrong when deal with df by mvm_validate.R
  
  #change model_rank from eg, 1,1,3,3,4,4 to 1,1,2,2,3,3 OR from 3,3,1,1,4,4 to 1,1,2,2,3,3
  Mod2Show=suppressMessages(left_join(Mod2Show,tibble(model_rank=unique(Mod2Show$model_rank)) %>% mutate(model_rank_new=1:n()),
                                      by='model_rank'))  %>%
    select(-model_rank) %>% select(model_rank=model_rank_new,everything())
  
  n_of_cov=length(unique(Mod2Show$cov_nam_with_linebreak))
  NAtibble=tibble(model_rank=rep(unique(Mod2Show$model_rank),each=n_of_cov),
                  cov_nam_with_linebreak=factor(rep(unique(Mod2Show$cov_nam_with_linebreak),n_mod2show),levels=levels(Mod2Show$cov_nam_with_linebreak))
                  )
  
  Mod2Show=suppressMessages(left_join(NAtibble,Mod2Show,by=c('model_rank','cov_nam_with_linebreak')))
  attr(Mod2Show,'RankOfLastMod2Show')=RankOfLastMod2Show
  Mod2Show
}

proc_mvm_interact_res=function(vTag,gender,RankBy='BIC',
                               dir_data,
                               dir_utils){
  #
  #ARGUMENTS
  #vTag & gender (Male/Female) : needed to import data
  #RankBy='AIC'/'BIC'
  #
  #NOTE
  ##can be run on cluster and Mac
  
  #########
  #prepare
  #########
  require(tidyverse)
  load(str_c(dir_data,'OutComb_',vTag,'_',gender,'_interact.RData'))
  eval(parse(text=str_c('data=OutCombo_',vTag,'_',gender,'_interact')))
  rm(list=ls(pattern='^OutCombo'))
  source(dir_utils)
  gender=if(gender=='Male') 'Men' else if(gender=='Female') 'Women'
  cutoff=get_cutoff(vTag)
  RfSource=get_RfSource(vTag)
  Var2Drop=str_c(get_Var2Drop(vTag),collapse=' & ')
  
  ########
  #misc processing
  ########
  data=data %>%
    group_by(model) %>% mutate(WithNA=if(any(is.na(IRR))) T else F) %>% ungroup %>% filter(WithNA==F) %>% select(-WithNA)
  
  n_of_mod=data$model %>% unique %>% length
  
  ########
  #create variables
  ########
  #NOTE: #dBIC can't be calculated in some instances (when a numeric interact with a factor, can't drop main eff of the numeric)
  AllModels=data %>% mutate(cov_nam_no_linebreak=recode_interact_cov_nam(covariate_name,inclu_line_break=F),
                            cov_nam_with_linebreak=recode_interact_cov_nam(covariate_name,inclu_line_break=T)) %>% {
    x=.
    if(RankBy=='AIC'){
      x %>% mutate(model_rank_x_main_eff=rank_mod(AIC))
    }else if(RankBy=='BIC'){
      x %>% mutate(model_rank_x_main_eff=rank_mod(BIC))
    }
  } %>% arrange(model_rank_x_main_eff) %>% mutate(IRRstr=turn_irr_to_str(IRR),IRR_low_str=turn_irr_to_str(IRR_low),IRR_up_str=turn_irr_to_str(IRR_up)) %>%
    select(model_rank=model_rank_x_main_eff,cov_nam_with_linebreak,IRR,IRR_low,IRR_up,IRRstr,IRR_low_str,IRR_up_str,
           pvalue=`p-value`,dBIC=delta_BIC,dAIC=delta_AIC,BIC,AIC,cov_nam_raw=covariate_name,cov_nam_no_linebreak,main_eff)
  
  
  ########
  #attr
  ########
  attr(AllModels,"gender")=gender
  attr(AllModels,'RankBy')=RankBy
  attr(AllModels,'cutoff')=cutoff
  attr(AllModels,'RfSource')=RfSource
  attr(AllModels,'vTag')=vTag
  attr(AllModels,'Number of models')=n_of_mod
  attr(AllModels,'GOF_max')=if(RankBy=='AIC') max(AllModels$AIC) else max(AllModels$BIC)
  attr(AllModels,'GOF_min')=if(RankBy=='AIC') min(AllModels$AIC) else min(AllModels$BIC)
  
  print(str_c('For ',vTag,' ',gender,' RfSource_',str_c(RfSource,collapse='&'),' cutoff',cutoff))
  print(str_c('Dropped ',Var2Drop))
  print(str_c(attr(AllModels,'Number of models'), ' models (allowing for interaction) have been tested.'))
  
  AllModels
}


##################
#test interaction for a given result set, defined by gender, cutoff, seer_version, RfSource, RankBy
##################
crc_mvm_tst_interact_v1=function(AllModels,rf_data,seer_data,max_dof=10,calc_dGOF=F,
                                 n_core,core_id,
                                 dir_utils,
                                 age_grp_bound=c(5,9),FixedPeriod=c(6,7,8,9,10),FixedAge=c(6,7,8,9),RtrnDemoEst=F,
                                 IncluPeriodEffect=F,IncluAgePerInter=F,IncluRacPerInter=F,IncluAgeRacInter=F,P_THRESHOLD=0.05
){
  #AIM
  #test interaction for all big models (not nested in a higher ranking model) in a given result set, defined by gender, cutoff, RfSource, RankBy
  
  #NOTE
  #this function is adapted from crc_mvm.R
  #only test two-way interaction
  #a model with interaction can include exposure-exposure (and/or exposure-demo_variable) interaction
  #handle amc models only
  
  #ARGUMENTS
  #AllModels returned by proc_mvm_res and then find_nest_mod, an attr of AllModels is parsim_mod_ls
  #max_dof: maximum # of degree-of-freedom for a model; dof=# of covariates (ignore intercept)
  #calc_dGOF default to F b/c 1) save time for calling Drop1StudyFactorThenCalcGOF AND 2) dGOF can't be computed for numeric main_eff if this main_eff
  #n_core: the number of cores used to run; n_core=if_else(loc_run==T,1,50)
  #core_id: the nth core used to run; core_id=if_else(loc_run,1,NTH_CORE); NTH_CORE assigned by bash script
  #continued: interact with a factor
  #a few arg use default to simplify code: age_grp_bound=c(5,9),FixedPeriod=c(6,7,8,9,10),FixedAge=c(6,7,8,9),RtrnDemoEst=F (setting for MultiReg function)
  #a few arg kept here but not used: IncluPeriodEffect=F,IncluAgePerInter=F,IncluRacPerInter=F,IncluAgeRacInter=F,P_THRESHOLD=0.05
  
  ##############
  #misc prepare
  ##############
  require(tidyverse)
  source(dir_utils)
  gender=attributes(AllModels)$gender
  mod_4_tst=attributes(AllModels)$parsim_mod_ls
  vTag=attributes(AllModels)$vTag
  ref_age=attributes(AllModels)$ref_grp[['age']]; ref_race=attributes(AllModels)$ref_grp['race']
  AllModels=NULL#relief space
  if(tolower(gender) %in% c("women","female")){
    TargetGender="Female"
  }else if(tolower(gender) %in% c("men","male")){
    TargetGender="Male"
  }else stop("Wrong gender value!")
  
  col_nam_PerAgeGenRac=c("period","age_grp","gender","race")
  col_nam_PerAgeRac=c("period","age_grp","race")
  std_fun=function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)
  
  ###########
  #get seer data for one gender
  ###########
  seer_data_use=seer_data %>% dplyr::filter(gender==TargetGender) %>% dplyr::select(-gender) %>% as_tibble
  
  ###########
  #get all risk factor data for one gender
  ###########
  all_study_factor_model_nameS=unique(unlist(mod_4_tst))
  #above, study factor nam and model nam are connected with _
  AllRF=str_split(all_study_factor_model_nameS,'_amc') %>% as.data.frame %>% .[1,] %>% unlist %>% {x=.;x=as.character(x);names(x)=NULL;x} %>% unique
  #note: unlist rtrn char on Mac and fctc on cluster; so coerce to char
  AllRfData=rf_data %>% dplyr::select(str_c(col_nam_PerAgeGenRac),str_c(AllRF)) %>%
    dplyr::filter(gender==TargetGender) %>% dplyr::select(-gender) %>% as_tibble
  ########################
  #the regression models
  ########################
  fun_no_lag=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    out=data.agg %>% dplyr::filter(age_grp %in% AgeForOutput,period %in% PeriodForOutput) %>%
      dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    #build a function to get quintiles
    get_quintile=function(x){
      out=x
      non.na.id=!is.na(x)
      #consider quintiles or tertiles
      breaks_for_cut=quantile(x[non.na.id],probs=(0:5)/5) %>% {x=.;x[1]=x[1]-0.1;x}#deal with low boundary problem
      if(length(unique(breaks_for_cut))<6){#some quintile cut points are the same
        out[non.na.id]=cut(x[non.na.id],breaks=quantile(x[non.na.id],probs=(0:3)/3) %>% {x=.;x[1]=x[1]-0.1;x})
        return(factor(out,levels=paste0(1:3)))
      }else{#the six quintile cut points are all different
        out[non.na.id]=cut(x[non.na.id],breaks=breaks_for_cut)#assign fctr to numeric, out becomes numeric with 1:5, including NAs
        return(factor(out,levels=paste0(1:5)))
      }
    }
    out[,TargetRF]=get_quintile(unlist(out[,TargetRF]))
    out
  }
  fun_lag=function(data.agg,TargetRF,UseWithOthFun=F,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    out=data.agg %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    out[,TargetRF]=NA_real_
    
    period_miss=base::setdiff(2:10,unique(out$period))
    if(length(period_miss)>0){
      for (period_miss_i in period_miss) {
        out=bind_rows(out,
                      out %>% dplyr::filter(period==out$period[1]) %>% mutate(period=period_miss_i))
      }
      out=out %>% distinct
    }
    for (i in 1:nrow(out)) {
      #fix period, age, race group, 
      tmp_demo=out %>% dplyr::select(str_c(col_nam_PerAgeRac)) %>% slice(i)
      #find period-2 and age-2 but same race
      data_ear=data.agg %>% dplyr::filter(period==tmp_demo$period-2,age_grp==tmp_demo$age_grp-2,race==tmp_demo$race) %>% dplyr::select(str_c(TargetRF))
      if(nrow(data_ear)==1) out[i,TargetRF]=data_ear#sometimes, there are NAs
    }
    if(UseWithOthFun==F){
      out %>% dplyr::filter(age_grp %in% AgeForOutput,period %in% PeriodForOutput) %>% fun_no_lag(.,TargetRF,PeriodForOutput,AgeForOutput)
    }else{
      out# out will be given to amc_lag2; age_grp and period bounds will be considered in that fun;
    }
  }
  fun_lag5=function(data.agg,TargetRF,UseWithOthFun=F,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    out=data.agg %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    out[,TargetRF]=NA_real_
    
    period_miss=base::setdiff(2:10,unique(out$period))
    if(length(period_miss)>0){
      for (period_miss_i in period_miss) {
        out=bind_rows(out,
                      out %>% dplyr::filter(period==out$period[1]) %>% mutate(period=period_miss_i))
      }
      out=out %>% distinct
    }
    
    for (i in 1:nrow(out)) {
      #fix period, age, race group, 
      tmp_demo=out %>% dplyr::select(str_c(col_nam_PerAgeRac)) %>% slice(i)
      #find period-1 and age-1 but same race
      data_ear=data.agg %>% dplyr::filter(period==tmp_demo$period-1,age_grp==tmp_demo$age_grp-1,race==tmp_demo$race) %>% dplyr::select(str_c(TargetRF))
      if(nrow(data_ear)==1) out[i,TargetRF]=data_ear#sometimes, there are NAs
    }
    
    if(UseWithOthFun==F){
      out %>% dplyr::filter(age_grp %in% AgeForOutput,period %in% PeriodForOutput) %>% fun_no_lag(.,TargetRF,PeriodForOutput,AgeForOutput)
    }else{
      out# out will be given to amc_lag52; age_grp and period bounds will be considered in that fun;
    }
  }
  fun_amc_no_lag=function(data.agg,TargetRF,UseWithOthFun=F,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    if(UseWithOthFun==F){
      out=data.agg %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput) %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    }else{
      out=data.agg %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
      # output will be given to amc_lag2 or amc_cum10y5ylag2; age_grp and period bounds will be considered in that fun
    }
    
    #group by age_grp and race
    #get mean
    age_race_mean=out %>% dplyr::select(-period) %>% group_by(age_grp,race) %>% summarise_all(mean,na.rm=T)
    #match the mean to out
    tmp=left_join(x=out,y=age_race_mean,by=c("age_grp"="age_grp","race"="race"))
    #minus mean for each col
    tmp.x=tmp %>% dplyr::select(contains(".x"))#don't contain demographic varaibles
    tmp.y=tmp %>% dplyr::select(contains(".y"))
    tmp=tmp.x-tmp.y#overwrite!
    if(UseWithOthFun==F) tmp=apply(tmp,2,std_fun)#standardize and overwrite!
    colnames(tmp)=TargetRF#b/c I changed colname earliers
    
    bind_cols(out[,col_nam_PerAgeRac],as_tibble(tmp))
  }
  fun_amc_lag2=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){ # do amc before lag
    out=fun_lag(fun_amc_no_lag(data.agg,TargetRF,UseWithOthFun=T),TargetRF,UseWithOthFun=T)
    out=out %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput)
    eval(parse(text=str_c('out$',TargetRF,'=std_fun(out$',TargetRF,')')))
    out
  }
  fun_amc_lag52=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){ # do amc before lag
    out=fun_lag5(fun_amc_no_lag(data.agg,TargetRF,UseWithOthFun=T),TargetRF,UseWithOthFun=T)
    out=out %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput)
    eval(parse(text=str_c('out$',TargetRF,'=std_fun(out$',TargetRF,')')))
    out
  }
  fun_amc_cum10y5ylag2=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    tmp=data.agg
    #fill in period 4 (deal with missingness; linear interpolation)
    #identify if needed; if so, fill in; do this for all rf
    tmp_tmp=tmp %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    tmp_tmp=tmp_tmp[!is.na(unlist(tmp_tmp[,TargetRF])),]#overwrite and keep rows not NA
    period_vec=unique(tmp_tmp$period)#periods with data for rf i
    if(all( (c(3,4,5) %in% period_vec)==c(T,F,T) )){
      tmp_tmp_pad_4=tmp_tmp %>% dplyr::filter(period %in% c(3,5)) %>% dplyr::select(-period) %>% 
        group_by(age_grp,race) %>% summarize_all(mean,na.rm=F) %>% dplyr::mutate(period=4)#na.rm doesn't matter b/c no NA
      tmp_tmp=bind_rows(tmp_tmp,tmp_tmp_pad_4)#overwrite!
    }
    tmp=tmp_tmp#overwrite! tmp can't be overwriten while for-loop was run
    out=tmp=fun_amc_no_lag(tmp,TargetRF,UseWithOthFun=T)
    out[,TargetRF]=NA_real_
    
    period_miss=base::setdiff(2:10,unique(out$period))
    if(length(period_miss)>0){
      for (period_miss_i in period_miss) {
        out=bind_rows(out,
                      out %>% dplyr::filter(period==out$period[1]) %>% mutate(period=period_miss_i))
      }
      out=out %>% distinct
    }
    
    for (i in 1:nrow(out)) {
      #fix period, age, race group, 
      tmp_demo=out %>% dplyr::select(str_c(col_nam_PerAgeRac)) %>% slice(i)
      # 5 yr lag, 10 y cum
      data_cum=tmp %>% dplyr::filter(
        (period==tmp_demo$period-1 & age_grp==tmp_demo$age_grp-1) | 
          (period==tmp_demo$period-2 & age_grp==tmp_demo$age_grp-2), race==tmp_demo$race) %>% dplyr::select(str_c(TargetRF))
      if(nrow(data_cum)==2)  out[i,TargetRF]=colSums(data_cum,na.rm=F)#sometimes, there are NAs
    }
    out=out %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput)
    eval(parse(text=str_c('out$',TargetRF,'=std_fun(out$',TargetRF,')')))
    out
  }
  
  ###########################
  #create data for all studyfactor-models
  ###########################
  PeriodForOutput=FixedPeriod
  AgeForOutput=FixedAge
  
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
    eval(parse(text=str_c("TmpData[[i]]=Fn_tmp(AllRfData,all_study_factor_nameS[i],PeriodForOutput=PeriodForOutput,AgeForOutput=AgeForOutput) %>%
                          rename(",all_study_factor_model_nameS[i],"=",all_study_factor_nameS[i],")")))
    if(i==1){
      DataAllStudyFactorModel=TmpData[[i]]
    } else{
      DataAllStudyFactorModel=suppressMessages(inner_join(DataAllStudyFactorModel,TmpData[[i]]))
    }
  }
  rm(TmpData)
  DataAllStudyFactorModel=DataAllStudyFactorModel %>% drop_na
  
  ######################
  #get data for regression by merging studyfactor-model data and seer data and dummy code demo var
  ######################
  DataForReg1=suppressMessages(inner_join(seer_data_use,DataAllStudyFactorModel))
  #dummy coding for age, period, and race
  DataForReg2=dummy_code_AgePerRac(DataForReg1,ref_age=ref_age,ref_race=ref_race)
  
  ################################################
  #build a function that take in rf_model_names and return non-interaction and interaction models in a list
  #################################################
  add_mod_with_InterTerm=function(study_factor_model_nameS,DataForReg,max_dof){
    #eg, study_factor_model_nameS=c("adj_avg_alc_m.brfss_amc_cum10y5ylag2","overweight_up_m.brfss_amc_no_lag",'thy_m_amc_no_lag') 
    #eg, retrun list(c("adj_avg...","overweight_up...",'thy...',"adj_avg...:overweight_up..."),
    #                c("adj_avg...","overweight_up...",'thy...',"thy...:overweight_up..."),
    #                c('adj_avg...','overweight_up...','thy...','adj_avg...:rBlack'),...)
    
    #ARUGUMENTS:
    #need attributes of DataForReg
    #max_dof: maximum # of degree-of-freedom for a model; dof=# of covariates + 1 (intercept)
    #
    #NOTE:
    #only consider two-way interaction, 
    #considered exposure-exposure interaction
    #considered exposure-age interaction, eg, adj_avg...:a7+adj_avg...:a8+adj_avg...:a9
    #continued: 'adj_avg...:a7+adj_avg...:a8+adj_avg...:a9' is one element in this vector, ie, out[[i]]
    #considered exposure-race interaction, eg, adj_avg...:rBlack
    
    NofMainEff=length(study_factor_model_nameS)#number of main effects
    #get all exposure-exposure interaction
    InterTerm1=if(NofMainEff>=2) apply(combn(study_factor_model_nameS,2),2,FUN=(function(x) str_c(x[1],':',x[2]))) else NULL
    #get all exposure-age interaction
    AllAgeTerm=str_split(attributes(DataForReg)$AllAgeTerm,'\\+')[[1]]#turn 1 element list to a vector
    InterTerm2=character(NofMainEff)
    for (i in 1:length(InterTerm2))  InterTerm2[i]=str_c(study_factor_model_nameS[i],AllAgeTerm,sep=':') %>% str_c(.,collapse='+')
    #get all exposure-race interaction
    InterTerm3=str_c(study_factor_model_nameS,attributes(DataForReg)$AllRacTerm,sep=':')
    
    #exhaust all interaction combos
    exhaust_combo_InterTerm=function(InterTerm,NofMainEff,max_dof){
      #
      #VALUE
      #out is initially a list (length=length(InterTerm)+1), 1st element is NA_character_, nth element is a vector (each containing n-1 covariate for interaction)
      #cont': last element is a vector (length=1; it contains length(InterTerm) covariates for interaction)
      #cont': 'adj_avg...:a7+adj_avg...:a8+adj_avg...:a9' is considered 1 covariate
      #we return unlist(out)
      #we cap length of out (as a list) to reduce burden on memory
      
      # minus 4 b/c of 3 age effects and 1 race effect
      max_n_of_inter_term=max(max_dof-NofMainEff-4,1)
      out_len=length(InterTerm)+1
      out_len=min(out_len,max_n_of_inter_term+1)
      out=as.list(rep(NA_character_,out_len))
      if(out_len>1){
        for (i in 2:out_len)  out[[i]]=apply(combn(InterTerm,i-1),2,FUN=(function(x) str_c(x,collapse='+')))
      }#if len==1 do nothing
      unlist(out)
    }
    
    #get combos
    df=expand.grid(exhaust_combo_InterTerm(InterTerm1,NofMainEff,max_dof),
                   exhaust_combo_InterTerm(InterTerm2,NofMainEff,max_dof),
                   exhaust_combo_InterTerm(InterTerm3,NofMainEff,max_dof))
    
    #get the list I want
    out=as.list(1:nrow(df))
    
    #concatenate main and interaction effects, getting a character vector
    #drop NAs in the character vector
    #drop models that are too complex; use NA to replace dropped models 
    
    for (i in 1:nrow(df)){
      out[[i]]=na.omit(c(study_factor_model_nameS,df[i,] %>% unlist %>% unname %>% as.character)) %>% {x=.;attributes(x)=NULL;x} %>%
        drop_complex_mod(.,max_dof=max_dof)
    }
    
    #drop NAs in the list
    out[!is.na(out)]
  }
  
  #############
  #build a multi var reg model for the 2-layer for-loop
  #############
  MultiReg=function(DataForReg,TargetGender,study_factor_model_nameS,calc_dGOF=F,
                    IncluPeriodEffect=F,IncluAgePerInter=F,IncluRacPerInter=F,IncluAgeRacInter=F,P_THRESHOLD=0.05,RtrnDemoEst=F){
    ###############################
    #create a function to calculate GOF, AIC&BIC when one study factor is dropped; this allows me to get delta_aic/bic later
    ###############################
    Drop1StudyFactorThenCalcGOF=function(data,study_factor_model_nameS,IncluPeriodEffect,FullRegCoefNum){
      #1) if ignore quintile or tertile models, AIC_vec/BIC_vec will initially be a list with len=len(study_factor_model_nameS), 
      #then will be a vec with len=len(study_factor_model_nameS); in such case, RepNum=1
      #2) if there is a quintile model, which happens to be dropped, RepNum=4, AIC_vec/BIC_vec will initially be a list with len=len(study_factor_model_nameS),
      #then will be a vec with longer length
      #3) if study_factor_model_nameS contain, eg, 'adj_avg...:a7+adj_avg...:a8+adj_avg...:a9', I break it into 3, and add to length(study_factor_model_nameS)
      
      study_factor_model_nameS=str_split(study_factor_model_nameS,'\\+') %>% unlist
      NumOfStudyFactor=length(study_factor_model_nameS)
      AIC_vec=BIC_vec=list()
      for (i in 1:NumOfStudyFactor) {
        if(IncluPeriodEffect==T){
          suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(data)$AllAgePerRacTerm,SigInterTerm,
                                                  "+",str_c(study_factor_model_nameS[-i],collapse="+"),"+offset(log(pop)),data=data)"))))
        }else{
          suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(data)$AllAgeTerm,"+",attributes(data)$AllRacTerm,
                                                  "+",str_c(study_factor_model_nameS[-i],collapse="+"),"+offset(log(pop)),data=data)"))))
        }
        irr_out=exp(coef(res_tmp))
        sam_siz=nrow(data)
        k_model=(length(irr_out)+1)#+1 free param for error
        RepNum=FullRegCoefNum-(k_model-1)#whether the study factor is amc or quintile or tertile; =1 for cont =3 or 4 otherwise
        BIC_vec[[i]]=if(RepNum==0) NA_real_ else rep(k_model*log(sam_siz)-summary(res_tmp)$twologlik,RepNum)
        AIC_vec[[i]]=if(RepNum==0) NA_real_ else rep(2*k_model-summary(res_tmp)$twologlik,RepNum)
      }
      list(AIC_vec=unlist(AIC_vec),BIC_vec=unlist(BIC_vec))
    }
    
    
    ########
    #run!!!
    ########
    suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg)$AllAgeTerm,"+",attributes(DataForReg)$AllRacTerm,
                                            "+",str_c(study_factor_model_nameS,collapse="+"),"+offset(log(pop)),data=DataForReg)"))))
    ########
    #output
    ########
    irr_out=exp(coef(res_tmp))
    irr_ci_out=suppressMessages(exp(confint(res_tmp)))
    irr_p_out=summary(res_tmp)$coefficients[,4]
    irr_p_out=if(length(irr_p_out)<length(irr_out)) NA_real_ else irr_p_out
    #when one pvalue is NA, the model is too saturated to be useful, so output NA for all pvalue
    
    if(calc_dGOF==F){
      GOF_vec=list(AIC_vec=NA_real_,BIC_vec=NA_real_)
    }else{
      GOF_vec=Drop1StudyFactorThenCalcGOF(data=DataForReg,study_factor_model_nameS=study_factor_model_nameS,IncluPeriodEffect=IncluPeriodEffect,
                                          FullRegCoefNum=length(irr_out))
    }
    
    sam_siz=nrow(DataForReg)
    k_model=(length(irr_out)+1)#+1 free param for error
    BIC=k_model*log(sam_siz)-summary(res_tmp)$twologlik
    AIC=2*k_model-summary(res_tmp)$twologlik
    
    periods_available=str_c(PeriodForOutput,collapse="+")
    num_periods_available=length(PeriodForOutput)
    
    out_demo_studyfactor=tibble(gender=TargetGender,covariate_name=names(irr_out),
                                IRR=irr_out,IRR_low=irr_ci_out[,1],IRR_up=irr_ci_out[,2],`p-value`=irr_p_out,BIC=BIC,AIC=AIC)#,
    
    if(RtrnDemoEst==F){
      out_demo_studyfactor %>% dplyr::filter(grepl("lag",covariate_name)) %>%#b/c all study factor name contain "lag", including exposure-demo interaction!!
        mutate(delta_AIC=GOF_vec$AIC_vec-AIC,delta_BIC=GOF_vec$BIC_vec-BIC) %>%
        dplyr::select(gender,covariate_name,IRR,IRR_low,IRR_up,`p-value`,contains("delta"),everything())
    }else{
      out_demo_studyfactor
    }
  }#end of MultiReg
  
  ##########
  #get all models ready for testing
  ##########
  #mod_ls_with_InterTerm is initially a list of lists
  #each element in up-level list correspond to a main effect model
  #each low-level list is a list of character vectors
  #mod_ls_with_InterTerm is then coerce to a list of character vectors
  mod_ls_with_InterTerm=as.list(1:length(mod_4_tst))
  for (i in 1:length(mod_4_tst)) {
    mod_ls_with_InterTerm[[i]]=add_mod_with_InterTerm(mod_4_tst[[i]],DataForReg=DataForReg2,max_dof=max_dof)
    #rhs returns a list, making mod_ls_with_InterTerm a list of lists
    #rhs may return a list with len=0, the following line handels that well
  }
  mod_ls_with_InterTerm=unlist(mod_ls_with_InterTerm,recursive=F,use.names=F)
  
  #############
  #loop thr all the combos for testing
  #allow parallel computing
  #############
  models=str_c('model',1:length(mod_ls_with_InterTerm))
  NumOfModel=length(models)
  
  if(n_core==1){
    ;
  }else{#n_core>=1
    NumOfModelEachCore=NumOfModel %/% n_core#the last core will run a bit more models
    if(core_id==n_core){#core_id is the last one, eg, 20 out of 1:20
      models=                              models[((core_id-1)*NumOfModelEachCore+1):length(models)]
      mod_ls_with_InterTerm=mod_ls_with_InterTerm[((core_id-1)*NumOfModelEachCore+1):length(mod_ls_with_InterTerm)]
      NumOfModel=length(models)
    }else{
      models=                              models[((core_id-1)*NumOfModelEachCore+1):(core_id*NumOfModelEachCore)]
      mod_ls_with_InterTerm=mod_ls_with_InterTerm[((core_id-1)*NumOfModelEachCore+1):(core_id*NumOfModelEachCore)]
      NumOfModel=length(models)
    }
    
  }
  
  OutCombo=NULL
  for (model_i in models) {
    i_on_this_core=which(models==model_i)
    if(i_on_this_core %% 100==0) {
      print(str_c("Test interaction: Finished ",i_on_this_core," models for ",vTag," and ",TargetGender,"; ",
                  NumOfModel," in total on this core; ",n_core," core(s)."))
    }
    study_factor_model_nameS=mod_ls_with_InterTerm[[i_on_this_core]]
    out=MultiReg(DataForReg=DataForReg2,TargetGender=TargetGender,study_factor_model_nameS=study_factor_model_nameS,RtrnDemoEst=RtrnDemoEst)
    out=out %>% mutate(model=model_i,main_eff=study_factor_model_nameS[!str_detect(study_factor_model_nameS,'\\:')] %>% str_c(.,collapse='\n'))
    OutCombo=bind_rows(OutCombo,out)
  }
  OutCombo
}

##################
#get predicted crc for one model
##################

crc_mvm4adj_crc=function(gender,
                         cov_nam_raw,
                         boot_n=1,
                         age_grp_bound=c(5,9),
                         IncluPeriodEffect=F,IncluAgePerInter=F,IncluRacPerInter=F,IncluAgeRacInter=F,P_THRESHOLD=0.05,
                         FixedPeriod=c(6,7,8,9,10),FixedAge=c(6,7,8,9),MinNumOfStudyFactor=1,
                         CorCoefThreshold=NULL,
                         dir_utils,dir
){
  ###############
  #misc prepare
  ###############
  require(tidyverse)
  if(tolower(gender) %in% c("women","female")){
    TargetGender="Female"
  }else if(tolower(gender) %in% c("men","male")){
    TargetGender="Male"
  }else stop("Wrong gender value!")
  
  col_nam_PerAgeGenRac=c("period","age_grp","gender","race")
  col_nam_PerAgeRac=c("period","age_grp","race")
  std_fun=function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T)
  source(dir_utils)
  
  ###########
  #get seer data for main and sensitivity analysis for one gender
  ###########
  load(str_c(dir,"p2_rf_list.Rdata"))
  rf_data=p2_rf_list %>% as_tibble
  seer_data=suppressMessages(read_csv(str_c(dir,"crc_3version.csv")) %>% mutate(version=if_else(version=="9_13_18","13_18",version)))
  seer_data_use=seer_data %>% dplyr::filter(version=='13',gender==TargetGender) %>% dplyr::select(-gender) %>% as_tibble
  
  ###########
  #get all risk factor data for one gender
  ###########
  all_study_factor_model_nameS=cov_nam_raw[str_detect(cov_nam_raw,'\\:',negate=T)]
  all_inter_term=cov_nam_raw[str_detect(cov_nam_raw,'\\:',negate=F)]
  all_study_factor_nameS=as.data.frame(str_split(all_study_factor_model_nameS,"_amc"))[1,] %>% unlist %>% {x=.;names(x)=NULL;x}
  all_model_nameS=as.data.frame(str_split(all_study_factor_model_nameS,"_amc"))[2,] %>% unlist %>% {x=.;names(x)=NULL;x}
  all_model_nameS=str_c('amc',all_model_nameS)
  AllRfData=rf_data %>% dplyr::select(all_of(col_nam_PerAgeGenRac),all_of(all_study_factor_nameS)) %>%
    dplyr::filter(gender==TargetGender) %>% dplyr::select(-gender) %>% as_tibble
  ########################
  #the regression models
  ########################
  fun_no_lag=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    out=data.agg %>% dplyr::filter(age_grp %in% AgeForOutput,period %in% PeriodForOutput) %>%
      dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    #build a function to get quintiles
    get_quintile=function(x){
      out=x
      non.na.id=!is.na(x)
      #consider quintiles or tertiles
      breaks_for_cut=quantile(x[non.na.id],probs=(0:5)/5) %>% {x=.;x[1]=x[1]-0.1;x}#deal with low boundary problem
      if(length(unique(breaks_for_cut))<6){#some quintile cut points are the same
        out[non.na.id]=cut(x[non.na.id],breaks=quantile(x[non.na.id],probs=(0:3)/3) %>% {x=.;x[1]=x[1]-0.1;x})
        return(factor(out,levels=paste0(1:3)))
      }else{#the six quintile cut points are all different
        out[non.na.id]=cut(x[non.na.id],breaks=breaks_for_cut)#assign fctr to numeric, out becomes numeric with 1:5, including NAs
        return(factor(out,levels=paste0(1:5)))
      }
    }
    out[,TargetRF]=get_quintile(unlist(out[,TargetRF]))
    out
  }
  fun_lag=function(data.agg,TargetRF,UseWithOthFun=F,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    out=data.agg %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    out[,TargetRF]=NA_real_
    
    period_miss=base::setdiff(2:10,unique(out$period))
    if(length(period_miss)>0){
      for (period_miss_i in period_miss) {
        out=bind_rows(out,
                      out %>% dplyr::filter(period==out$period[1]) %>% mutate(period=period_miss_i))
      }
      out=out %>% distinct
    }
    for (i in 1:nrow(out)) {
      #fix period, age, race group, 
      tmp_demo=out %>% dplyr::select(str_c(col_nam_PerAgeRac)) %>% slice(i)
      #find period-2 and age-2 but same race
      data_ear=data.agg %>% dplyr::filter(period==tmp_demo$period-2,age_grp==tmp_demo$age_grp-2,race==tmp_demo$race) %>% dplyr::select(str_c(TargetRF))
      if(nrow(data_ear)==1) out[i,TargetRF]=data_ear#sometimes, there are NAs
    }
    if(UseWithOthFun==F){
      out %>% dplyr::filter(age_grp %in% AgeForOutput,period %in% PeriodForOutput) %>% fun_no_lag(.,TargetRF,PeriodForOutput,AgeForOutput)
    }else{
      out# out will be given to amc_lag2; age_grp and period bounds will be considered in that fun;
    }
  }
  fun_lag5=function(data.agg,TargetRF,UseWithOthFun=F,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    out=data.agg %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    out[,TargetRF]=NA_real_
    
    period_miss=base::setdiff(2:10,unique(out$period))
    if(length(period_miss)>0){
      for (period_miss_i in period_miss) {
        out=bind_rows(out,
                      out %>% dplyr::filter(period==out$period[1]) %>% mutate(period=period_miss_i))
      }
      out=out %>% distinct
    }
    
    for (i in 1:nrow(out)) {
      #fix period, age, race group, 
      tmp_demo=out %>% dplyr::select(str_c(col_nam_PerAgeRac)) %>% slice(i)
      #find period-1 and age-1 but same race
      data_ear=data.agg %>% dplyr::filter(period==tmp_demo$period-1,age_grp==tmp_demo$age_grp-1,race==tmp_demo$race) %>% dplyr::select(str_c(TargetRF))
      if(nrow(data_ear)==1) out[i,TargetRF]=data_ear#sometimes, there are NAs
    }
    
    if(UseWithOthFun==F){
      out %>% dplyr::filter(age_grp %in% AgeForOutput,period %in% PeriodForOutput) %>% fun_no_lag(.,TargetRF,PeriodForOutput,AgeForOutput)
    }else{
      out# out will be given to amc_lag52; age_grp and period bounds will be considered in that fun;
    }
  }
  fun_amc_no_lag=function(data.agg,TargetRF,UseWithOthFun=F,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    if(UseWithOthFun==F){
      out=data.agg %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput) %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    }else{
      out=data.agg %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
      # output will be given to amc_lag2 or amc_cum10y5ylag2; age_grp and period bounds will be considered in that fun
    }
    
    #group by age_grp and race
    #get mean
    age_race_mean=out %>% dplyr::select(-period) %>% group_by(age_grp,race) %>% summarise_all(mean,na.rm=T)
    #match the mean to out
    tmp=left_join(x=out,y=age_race_mean,by=c("age_grp"="age_grp","race"="race"))
    #minus mean for each col
    tmp.x=tmp %>% dplyr::select(contains(".x"))#don't contain demographic varaibles
    tmp.y=tmp %>% dplyr::select(contains(".y"))
    tmp=tmp.x-tmp.y#overwrite!
    if(UseWithOthFun==F) tmp=apply(tmp,2,std_fun)#standardize and overwrite!
    colnames(tmp)=TargetRF#b/c I changed colname earliers
    
    bind_cols(out[,col_nam_PerAgeRac],as_tibble(tmp))
  }
  fun_amc_lag2=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){ # do amc before lag
    out=fun_lag(fun_amc_no_lag(data.agg,TargetRF,UseWithOthFun=T),TargetRF,UseWithOthFun=T)
    out=out %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput)
    eval(parse(text=str_c('out$',TargetRF,'=std_fun(out$',TargetRF,')')))
    out
  }
  fun_amc_lag52=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){ # do amc before lag
    out=fun_lag5(fun_amc_no_lag(data.agg,TargetRF,UseWithOthFun=T),TargetRF,UseWithOthFun=T)
    out=out %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput)
    eval(parse(text=str_c('out$',TargetRF,'=std_fun(out$',TargetRF,')')))
    out
  }
  fun_amc_cum10y5ylag2=function(data.agg,TargetRF,PeriodForOutput=3:10,AgeForOutput=age_grp_bound[1]:age_grp_bound[2]){
    tmp=data.agg
    #fill in period 4 (deal with missingness; linear interpolation)
    #identify if needed; if so, fill in; do this for all rf
    tmp_tmp=tmp %>% dplyr::select(str_c(col_nam_PerAgeRac),str_c(TargetRF))
    tmp_tmp=tmp_tmp[!is.na(unlist(tmp_tmp[,TargetRF])),]#overwrite and keep rows not NA
    period_vec=unique(tmp_tmp$period)#periods with data for rf i
    if(all( (c(3,4,5) %in% period_vec)==c(T,F,T) )){
      tmp_tmp_pad_4=tmp_tmp %>% dplyr::filter(period %in% c(3,5)) %>% dplyr::select(-period) %>% 
        group_by(age_grp,race) %>% summarize_all(mean,na.rm=F) %>% dplyr::mutate(period=4)#na.rm doesn't matter b/c no NA
      tmp_tmp=bind_rows(tmp_tmp,tmp_tmp_pad_4)#overwrite!
    }
    tmp=tmp_tmp#overwrite! tmp can't be overwriten while for-loop was run
    out=tmp=fun_amc_no_lag(tmp,TargetRF,UseWithOthFun=T)
    out[,TargetRF]=NA_real_
    
    period_miss=base::setdiff(2:10,unique(out$period))
    if(length(period_miss)>0){
      for (period_miss_i in period_miss) {
        out=bind_rows(out,
                      out %>% dplyr::filter(period==out$period[1]) %>% mutate(period=period_miss_i))
      }
      out=out %>% distinct
    }
    
    for (i in 1:nrow(out)) {
      #fix period, age, race group, 
      tmp_demo=out %>% dplyr::select(str_c(col_nam_PerAgeRac)) %>% slice(i)
      # 5 yr lag, 10 y cum
      data_cum=tmp %>% dplyr::filter(
        (period==tmp_demo$period-1 & age_grp==tmp_demo$age_grp-1) | 
          (period==tmp_demo$period-2 & age_grp==tmp_demo$age_grp-2), race==tmp_demo$race) %>% dplyr::select(str_c(TargetRF))
      if(nrow(data_cum)==2)  out[i,TargetRF]=colSums(data_cum,na.rm=F)#sometimes, there are NAs
    }
    out=out %>% dplyr::filter(age_grp %in% AgeForOutput, period %in% PeriodForOutput)
    eval(parse(text=str_c('out$',TargetRF,'=std_fun(out$',TargetRF,')')))
    out
  }
  ###########################
  #create data for all studyfactor-models
  ###########################
  PeriodForOutput=FixedPeriod
  AgeForOutput=FixedAge
  
  for (i in 1:length(all_study_factor_model_nameS)) {
    Fn_tmp=get(str_c("fun_",all_model_nameS[i]))
    #TmpDatai=Fn_tmp(AllRfData,all_study_factor_nameS[i],PeriodForOutput=PeriodForOutput,AgeForOutput=AgeForOutput) %>%
    #rename(all_study_factor_model_nameS[i]=all_study_factor_nameS[i])
    eval(parse(text=str_c("TmpData",i,"=Fn_tmp(AllRfData,all_study_factor_nameS[i],PeriodForOutput=PeriodForOutput,AgeForOutput=AgeForOutput) %>%
                          rename(",all_study_factor_model_nameS[i],"=",all_study_factor_nameS[i],")")))
    if(i==1){
      DataAllStudyFactorModel=TmpData1
    } else{
      #DataAllStudyFactorModel=inner_join(DataAllStudyFactorModel,TmpDatai)
      suppressMessages(eval(parse(text=str_c("DataAllStudyFactorModel=inner_join(DataAllStudyFactorModel,TmpData",i,")"))))
    }
  }
  DataAllStudyFactorModel=DataAllStudyFactorModel %>% drop_na
  
  
  #############
  #build a multi var reg model for the big for-loop
  #############
  MultiReg=function(rf_data_use,seer_data_use,TargetGender,study_factor_model_nameS,InterTerm=InterTerm,boot_n,
                    IncluPeriodEffect=F,IncluAgePerInter=F,IncluRacPerInter=F,IncluAgeRacInter=F,P_THRESHOLD=0.05,RtrnDemoEst=F){
    ##########################################################
    #create a function for dummy coding age, period, and race#
    ##########################################################
    dummy_code_AgePerRac=function(df){
      periods=sort(unique(df$period))[-1]#no dummy coding for reference level
      age_grps=sort(unique(df$age_grp))[-1]
      for (i_periods in periods) {
        #df$pi_periods=if_else(df$period==i_periods,'1','0') %>% factor
        eval(parse(text=str_c("df$p",i_periods,"=if_else(df$period==i_periods,'1','0') %>% factor")))
      }
      for (i_age_grps in age_grps) {
        #df$ai_age_grps=if_else(df$age_grp==i_age_grps,'1','0')
        eval(parse(text=str_c("df$a",i_age_grps,"=if_else(df$age_grp==i_age_grps,'1','0') %>% factor")))
      }
      df$rBlack=if_else(df$race=="Black",'1','0') %>% factor
      
      AllAgeTerm=str_c(str_c("a",age_grps),collapse="+")#works even if length(age_grps) is 1
      AllPerTerm=str_c(str_c("p",periods ),collapse="+")
      AllAgePerRacTerm=str_c(AllAgeTerm,"+",AllPerTerm,"+rBlack")
      attr(df,"AllAgePerRacTerm")=AllAgePerRacTerm
      attr(df,"AllAgeTerm")=AllAgeTerm
      attr(df,"AllPerTerm")=AllPerTerm
      attr(df,"AllRacTerm")="rBlack"
      df
    }
    ########################################################
    #create a function to get interaction terms for testing#
    ########################################################
    get_inter_term=function(df,IncluAgePerInter,IncluRacPerInter,IncluAgeRacInter){
      if(IncluAgePerInter==T & IncluRacPerInter==F & IncluAgeRacInter==F){
        InterTerm=str_c("(",attributes(df)$AllAgeTerm,")*(",attributes(df)$AllPerTerm,")")
      }else if(IncluAgePerInter==T & IncluRacPerInter==T & IncluAgeRacInter==F){
        InterTerm=str_c("(",attributes(df)$AllAgeTerm,")*(",attributes(df)$AllPerTerm,")+(",
                        attributes(df)$AllPerTerm,")*rBlack")
      }else if(IncluAgePerInter==T & IncluRacPerInter==T & IncluAgeRacInter==T){
        InterTerm=str_c("(",attributes(df)$AllAgeTerm,")*(",attributes(df)$AllPerTerm,")+(",
                        attributes(df)$AllPerTerm,")*rBlack+(",attributes(df)$AllAgeTerm,")*rBlack")
      }else if(IncluAgePerInter==F & IncluRacPerInter==F & IncluAgeRacInter==T){
        InterTerm=str_c('(',attributes(df)$AllAgeTerm,')*rBlack')
      }else if(IncluAgePerInter==F & IncluRacPerInter==F & IncluAgeRacInter==F){
        InterTerm=NULL
      }else{
        stop("Wrong value for IncluAgePerInter, IncluRacPerInter, or IncluAgeRacInter!")
      }
      InterTerm
    }
    ###############################
    #create a function to calculate GOF, AIC&BIC when one study factor is dropped; this allows me to get delta_aic/bic later
    ###############################
    Drop1StudyFactorThenCalcGOF=function(data,study_factor_model_nameS,IncluPeriodEffect,FullRegCoefNum){
      #1) if ignore quintile or tertile models, AIC_vec/BIC_vec will initially be a list with len=len(study_factor_model_nameS), 
      #then will be a vec with len=len(study_factor_model_nameS); in such case, RepNum=1
      #2) if there is a quintile model, which happens to be dropped, RepNum=4, AIC_vec/BIC_vec will initially be a list with len=len(study_factor_model_nameS),
      #then will be a vec with longer length
      
      NumOfStudyFactor=length(study_factor_model_nameS)
      AIC_vec=BIC_vec=list()
      for (i in 1:NumOfStudyFactor) {
        if(IncluPeriodEffect==T){
          suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(data)$AllAgePerRacTerm,SigInterTerm,
                                                  "+",str_c(study_factor_model_nameS[-i],collapse="+"),"+offset(log(pop)),data=data)"))))
        }else{
          suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(data)$AllAgeTerm,"+",attributes(data)$AllRacTerm,
                                                  "+",str_c(study_factor_model_nameS[-i],collapse="+"),"+offset(log(pop)),data=data)"))))
        }
        irr_out=exp(coef(res_tmp))
        sam_siz=nrow(data)
        k_model=(length(irr_out)+1)#+1 free param for error
        RepNum=FullRegCoefNum-(k_model-1)#whether the study factor is amc or quintile or tertile
        BIC_vec[[i]]=rep(k_model*log(sam_siz)-summary(res_tmp)$twologlik,RepNum)
        AIC_vec[[i]]=rep(2*k_model-summary(res_tmp)$twologlik,RepNum)
      }
      list(AIC_vec=unlist(AIC_vec),BIC_vec=unlist(BIC_vec))
    }
    
    ######################
    #merge rf_data_use and seer_data_use
    ######################
    DataForReg1=suppressMessages(inner_join(seer_data_use,rf_data_use))
    
    ########
    #run!!!
    ########
    #dummy coding for age, period, and race
    DataForReg2=dummy_code_AgePerRac(DataForReg1)
    
    #assume IncluPeriodEffect==F
    suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgeTerm,"+",attributes(DataForReg2)$AllRacTerm,
                                            "+",str_c(study_factor_model_nameS,collapse="+"),"+",str_c(InterTerm,collapse='+'),
                                            "+offset(log(pop)),data=DataForReg2)"))))
    
    #output
    irr_out=exp(coef(res_tmp))
    irr_ci_out=suppressMessages(exp(confint(res_tmp)))
    irr_p_out=summary(res_tmp)$coefficients[,4]
    irr_p_out=if(length(irr_p_out)<length(irr_out)) NA_real_ else irr_p_out
    #when one pvalue is NA, the model is too saturated to be useful, so output NA for all pvalue
    
    sam_siz=nrow(DataForReg2)
    k_model=(length(irr_out)+1)#+1 free param for error
    BIC=k_model*log(sam_siz)-summary(res_tmp)$twologlik
    AIC=2*k_model-summary(res_tmp)$twologlik
    
    periods_available=str_c(PeriodForOutput,collapse="+")
    num_periods_available=length(PeriodForOutput)
    
    #need to store results of regression
    res_reg=list()
    res_reg[['res_reg']]=res_tmp
    res_reg[['des_mat']]=model.matrix(res_tmp)
    
    #get predicted crc, and adjusted predicted crc
    count_pred_mu=res_tmp$fitted.values %>% unname
    disp_size=res_tmp$theta
    
    out=NULL
    for (boot_i in 1:boot_n) {
      count_pred_adj=list()
      count_pred_adj$count_pred=rnbinom(n=40,mu=count_pred_mu,size=disp_size)
      
      for (study_factor_model_i in study_factor_model_nameS) {
        design_mat=model.matrix(res_tmp)
        col_nam_des_mat=colnames(design_mat)
        lgl_vec_T_rf_mod=(col_nam_des_mat==study_factor_model_i)
        lgl_vec_T_contain_rf_mod=str_detect(col_nam_des_mat,pattern=study_factor_model_i)
        
        design_mat[,lgl_vec_T_rf_mod]=min(design_mat[,lgl_vec_T_rf_mod])#set rf to minimal
        design_mat[,(!lgl_vec_T_rf_mod) & lgl_vec_T_contain_rf_mod]=0#set interaction term to zero
        res_reg[[str_c('des_mat_',study_factor_model_i)]]=design_mat
        
        count_pred_adj_tmp=(design_mat %*% matrix(coef(res_tmp),ncol=1))+res_tmp$model[,'offset(log(pop))']
        count_pred_adj_tmp=rnbinom(n=40,mu=exp(as.vector(count_pred_adj_tmp)),size=disp_size)
        #count_pred_adj$study_factor_model_i=count_pred_adj_tmp
        eval(parse(text=str_c('count_pred_adj$count_pred_adj_',study_factor_model_i,'=count_pred_adj_tmp')))
      }
      
      out=bind_rows(out,bind_cols(gender=TargetGender,DataForReg1,boot_i=boot_i,as_tibble(count_pred_adj)))
      #if(boot_i %% 100 ==0) print(str_c('Finished bootstrap iteraction ',boot_i,'!'))
    }
    attr(out,'res_reg')=res_reg
    out
  }
  
  MultiReg(rf_data_use=DataAllStudyFactorModel,seer_data_use=seer_data_use,TargetGender=TargetGender,
           study_factor_model_nameS=all_study_factor_model_nameS,InterTerm=all_inter_term,boot_n=boot_n) 
  
}

