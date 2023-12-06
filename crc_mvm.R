
crc_mvm=function(rf_data,seer_data,gender,
                 StudyFactorModelCombo,
                 age_grp_bound=c(5,9),
                 IncluPeriodEffect=F,IncluAgePerInter=F,IncluRacPerInter=F,IncluAgeRacInter=F,P_THRESHOLD=0.05,
                 FixedPeriod=NULL,FixedAge=NULL,ref_age,ref_race,
                 CorCoefThreshold=NULL,
                 MinNumOfStudyFactor=2,
                 n_core,core_id,
                 dir_utils
){
  ##NOTE
  #need to install tidyverse, MASS
  ##ARGUMENTS:
  #n_core: the number of cores used to run; n_core=if_else(loc_run==T,1,20)
  #core_id: the nth core used to run; core_id=if_else(loc_run,1,NTH_CORE); NTH_CORE assigned by bash script
  #dir_utils=".../utils.R"
  
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
  #get seer data for one gender
  ###########
  seer_data_use=seer_data %>% dplyr::filter(gender==TargetGender) %>% dplyr::select(-gender) %>% as_tibble
  
  ###########
  #get all risk factor data for one gender
  ###########
  AllRF=GetAllRF(StudyFactorModelCombo,TargetGender)#exclude thyroid for Male; only study factor name returned here
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
  all_study_factor_model_nameS=unlist(StudyFactorModelCombo)
  all_study_factor_model_nameS=all_study_factor_model_nameS[str_detect(all_study_factor_model_nameS,"lag")]#filter out elements with "-" only
  if(TargetGender=="Female"){
    all_study_factor_model_nameS=all_study_factor_model_nameS
  }else{#don't do thyroid for Male
    all_study_factor_model_nameS=all_study_factor_model_nameS[!str_detect(all_study_factor_model_nameS,"thy_m")]
  }
  
  all_study_factor_nameS=as.data.frame(str_split(all_study_factor_model_nameS,"-"))[1,] %>% unlist %>% {x=.;x=as.character(x);names(x)=NULL;x}#may have duplicates
  all_model_nameS=as.data.frame(str_split(all_study_factor_model_nameS,"-"))[2,] %>% unlist %>% {x=.;x=as.character(x);names(x)=NULL;x}#may have duplicates
  #note: unlist rtrn char on Mac and fctc on cluster; so coerce to char
  #the above two var are meant to contain duplicates, eg,
  #adj_avg_alc_m.brfss, adj_avg_alc_m.brfss, fibe_m.fill
  #amc_lag2, amc_cum10y5ylag2, amc_cum10y5ylag2
  
  all_study_factor_model_nameS=str_replace(all_study_factor_model_nameS,"-","_")#replace the symbol between study
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
    print(str_c('Finished exposure-model number ',i,'!'))
    print(str_c(round(as.numeric(Sys.time()-tm_before,units="mins"),2), " minutes have past!"))
  }
  rm(TmpData)
  DataAllStudyFactorModel=DataAllStudyFactorModel %>% drop_na
  
  #############
  #create pairwise correlation coef matrix among study factors
  #############
  CorMat=DataAllStudyFactorModel %>% select(str_c(all_study_factor_model_nameS)) %>% mutate_all(.funs=as.numeric) %>% cor(.,method="spearman")
  diag(CorMat)=0
  
  #############
  #create a matrix holding correlation coef of study factors with period (1st col), age (2nd col) and with raceBlack (3rd col)
  #############
  CorMatDem=DataAllStudyFactorModel %>% mutate_all(.funs=as.numeric) %>% {
    x=.
    out=matrix(0,nrow=length(all_study_factor_model_nameS),ncol=3)
    colnames(out)=colnames(x)[1:3]#period, age_grp, race
    rownames(out)=all_study_factor_model_nameS
    #out stores correlation b/w a study factor and a demo variables
    for (DemoVar in colnames(x)[1:3]) {
      for (RfMod in all_study_factor_model_nameS) {
        out[RfMod,DemoVar]=cor(x %>% dplyr::select(str_c(DemoVar)) %>% unlist, x %>% dplyr::select(str_c(RfMod)) %>% unlist, method='spearman')
      }
    }
    out
  }
  
  #############
  #build a multi var reg model for the for-loop for all exposure-model combinations
  #############
  MultiReg=function(rf_data_use,seer_data_use,TargetGender,study_factor_model_nameS,ref_age,ref_race,
                    IncluPeriodEffect=F,IncluAgePerInter=F,IncluRacPerInter=F,IncluAgeRacInter=F,P_THRESHOLD=0.05,RtrnDemoEst=F){
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
    DataForReg2=dummy_code_AgePerRac(DataForReg1,ref_age=ref_age,ref_race=ref_race)
    
    if(IncluPeriodEffect==T){
      #if do interaction: do regression to test significant interaction terms
      InterTerm=get_inter_term(DataForReg2,IncluAgePerInter,IncluRacPerInter,IncluAgeRacInter)#get all interaction terms for testing
      if(any(c(IncluAgePerInter,IncluAgeRacInter,IncluRacPerInter))==T){
        suppressWarnings(eval(parse(text=str_c("res_test_inter=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgePerRacTerm,"+",InterTerm,
                                               "+offset(log(pop)),data=DataForReg2)"))))
        #select significant interaction terms
        #grab the matrix holding param names (as rownames) and p-values, etc
        SigInterTerm=summary(res_test_inter)$coefficients %>% {x=.;as_tibble(x) %>% mutate(param_name=rownames(x))} %>%
          filter(`Pr(>|z|)`<P_THRESHOLD) %>% filter(grepl(":",param_name)) %>% .$param_name %>%
          str_replace(.,"1:","*") %>% str_sub(.,1,-2) %>% str_c(.,collapse="+") %>% str_c("+",.)
      }else SigInterTerm=NULL
      
      #regression
      suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgePerRacTerm,
                                              SigInterTerm,
                                              "+",str_c(study_factor_model_nameS,collapse="+"),"+offset(log(pop)),data=DataForReg2)"))))
    }else{#IncluPeriodEffect==F
      #regression
      suppressWarnings(eval(parse(text=paste0("res_tmp=MASS::glm.nb(count~",attributes(DataForReg2)$AllAgeTerm,"+",attributes(DataForReg2)$AllRacTerm,
                                              "+",str_c(study_factor_model_nameS,collapse="+"),"+offset(log(pop)),data=DataForReg2)"))))
    }
    
    #output
    irr_out=exp(coef(res_tmp))
    irr_ci_out=suppressMessages(exp(confint(res_tmp)))
    irr_p_out=summary(res_tmp)$coefficients[,4]
    irr_p_out=if(length(irr_p_out)<length(irr_out)) NA_real_ else irr_p_out
    #when one pvalue is NA, the model is too saturated to be useful, so output NA for all pvalue
    
    GOF_vec=Drop1StudyFactorThenCalcGOF(data=DataForReg2,study_factor_model_nameS=study_factor_model_nameS,IncluPeriodEffect=IncluPeriodEffect,
                                        FullRegCoefNum=length(irr_out))
    
    sam_siz=nrow(DataForReg2)
    k_model=(length(irr_out)+1)#+1 free param for error
    BIC=k_model*log(sam_siz)-summary(res_tmp)$twologlik
    AIC=2*k_model-summary(res_tmp)$twologlik
    
    periods_available=str_c(PeriodForOutput,collapse="+")
    num_periods_available=length(PeriodForOutput)
    
    if(IncluPeriodEffect==T){
      if(any(c(IncluAgePerInter,IncluAgeRacInter,IncluRacPerInter))==F) SigInterTerm=" None"#overwrite NULL for output
    }else{
      SigInterTerm=" None"
    }
    out_demo_studyfactor=tibble(gender=TargetGender,covariate_name=names(irr_out),
                                IRR=irr_out,IRR_low=irr_ci_out[,1],IRR_up=irr_ci_out[,2],`p-value`=irr_p_out,BIC=BIC,AIC=AIC,
                                `Periods available`=periods_available,`Num of periods available`=num_periods_available,
                                `Interaction term`=str_sub(SigInterTerm,2,-1))
    if(RtrnDemoEst==F){
      out_demo_studyfactor %>% dplyr::filter(grepl("lag",covariate_name)) %>%#b/c all study factor name contain "lag"!!
        mutate(delta_AIC=GOF_vec$AIC_vec-AIC,delta_BIC=GOF_vec$BIC_vec-BIC) %>%
        dplyr::select(gender,covariate_name,IRR,IRR_low,IRR_up,`p-value`,contains("delta"),everything())
    }else{
      out_demo_studyfactor
    }
  }
  
  
  #############
  #get all study_factor_model_nameS in a list
  #############
  if(gender=="Male") StudyFactorModelCombo[6]=list("-")#don't do thyroid
  max_n_of_rf_mod=length(StudyFactorModelCombo)#both =11 for men and women
  StudyFactorAndModel=do.call(crossing,StudyFactorModelCombo)
  StudyFactorAndModel=t(as.matrix(StudyFactorAndModel))
  if(all(StudyFactorAndModel[,1]==rep('-',max_n_of_rf_mod))){
    #if the first col contain no exposure-model, ie, all being '-'
    StudyFactorAndModel=StudyFactorAndModel[,-1]
  }else{
    stop('Error! The 1st col of StudyFactorAndModel are not all = nothing.')
  }
  print('Got the big matrix StudyFactorAndModel!')
  print(str_c(round(as.numeric(Sys.time()-tm_before,units="mins"),2), " minutes have past!"))
  NumOfModel=ncol(StudyFactorAndModel)
  print(str_c('The total number of models or combinations is ',NumOfModel,'!'))
  models=str_c("model",1:NumOfModel)
  colnames(StudyFactorAndModel)=models
  
  #distribute computing
  if(n_core==1){
    ;#don't overwrite models & NumOfModel
  }else{#n_core>=1
    NumOfModelEachCore=NumOfModel %/% n_core#the last core will run a bit more models
    if(core_id==n_core){#core_id is the last one, eg, 20 out of 1:20
      models=models[((core_id-1)*NumOfModelEachCore+1):length(models)]
      NumOfModel=length(models)
    }else{
      models=models[((core_id-1)*NumOfModelEachCore+1):(core_id*NumOfModelEachCore)]
      NumOfModel=length(models)
    }
    
  }
  
  #subset and REDUCE WEIGHT!
  StudyFactorAndModel=StudyFactorAndModel[,models]
  #replace '-' with NA_character
  StudyFactorAndModel[StudyFactorAndModel=='-']=NA_character_
  #replace exposure-model with exposure_model
  StudyFactorAndModel=matrix(str_replace(StudyFactorAndModel,'-','_'),nrow=max_n_of_rf_mod)
  colnames(StudyFactorAndModel)=models
  study_factor_model_nameS_ls=as.list(as_tibble(StudyFactorAndModel))
  names(study_factor_model_nameS_ls)=models
  rm(StudyFactorAndModel)#save space
  print('Got the big list study_factor_model_nameS_ls and start to loop through models!')
  print(str_c(round(as.numeric(Sys.time()-tm_before,units="mins"),2), " minutes have past!"))
  
  #############
  #loop thr all the combos for testing
  #############
  OutCombo=NULL
  for (model_i in models) {
    if(as.numeric(str_sub(model_i,6,-1)) %% 100==0) {
      print(str_c("Finished ",which(models==model_i)," for ",gender,"; ",NumOfModel," in total on this core; ",n_core," core(s)."))
    }
    #for each combo/model, test whether I need to run regression
    study_factor_model_nameS=na.omit(study_factor_model_nameS_ls[[model_i]]) %>% {x=.;attributes(x)=NULL;x}
    max_abs_cor_coef=max(abs(CorMat[study_factor_model_nameS,study_factor_model_nameS]))
    if(max_abs_cor_coef>CorCoefThreshold |
       max(abs(CorMatDem[study_factor_model_nameS,col_nam_PerAgeRac[2:3]]))>CorCoefThreshold) next
    
    out=MultiReg(rf_data_use=DataAllStudyFactorModel,seer_data_use=seer_data_use,TargetGender=TargetGender,
                 study_factor_model_nameS=study_factor_model_nameS,ref_age=ref_age,ref_race=ref_race)
    
    OutCombo=bind_rows(OutCombo,out %>% mutate(max_abs_cor_coef=max_abs_cor_coef,model=model_i))
  }
  OutCombo
}



