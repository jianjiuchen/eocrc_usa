###################
#used in crc_mvm_v* and crc_mvm_tst_interact_v*
###################
GetAllRF=function(StudyFactorModelCombo,TargetGender){
  all_study_factor_model_nameS=unlist(StudyFactorModelCombo)
  all_study_factor_model_nameS=all_study_factor_model_nameS[str_detect(all_study_factor_model_nameS,"lag")]#filter out elements with "-" only
  all_study_factor_nameS=as.data.frame(str_split(all_study_factor_model_nameS,"-"))[1,] %>% unlist %>% {x=.;names(x)=NULL;x}
  uniqu_study_factor_nameS=unique(all_study_factor_nameS)
  #consider gender
  if(TargetGender=="Female") uniqu_study_factor_nameS else uniqu_study_factor_nameS[!str_detect(uniqu_study_factor_nameS,"thy_m")]
}

drop_complex_mod=function(all_cov,max_dof){
  #ARGUMENTS
  #all_cov (all covariates) could be c('obesity','obesity*rBlack','obesity*a71+obesity*a81+obesity*a91')
  #max_dof: maximum # of degree-of-freedom for a model; dof=# of covariates + 1 (intercept)
  #
  #NOTE
  #dof here is a synanym of # of covariate; so I ignore the intercept when calculating dof
  
  dof=(str_c(all_cov,collapse='+') %>% str_count(.,'\\+'))+5
  #+1 b/c there is 1 less + than cov; 
  #+4 b/c there is a race effect and 3 age effect
  
  if(dof<=max_dof){
    all_cov
  }else NA_character_
}

#create a function for dummy coding age, period, and race
dummy_code_AgePerRac=function(df,ref_age,ref_race){
  #ARUGMENT
  #df: obj by joining seer data and rf data; before dummy coding for age, period, and race
  #ref_age: if =1, ref_grp is age 30-34; if =4, ref_grp is age 45-49
  #ref_race: set to 'Black' or 'White'
  
  periods=sort(unique(df$period))[-1]#no dummy coding for reference level
  age_grps=sort(unique(df$age_grp))[-ref_age]
  for (i_periods in periods) {
    #df$pi_periods=if_else(df$period==i_periods,'1','0') %>% factor
    eval(parse(text=str_c("df$p",i_periods,"=if_else(df$period==i_periods,'1','0') %>% factor")))
  }
  for (i_age_grps in age_grps) {
    #df$ai_age_grps=if_else(df$age_grp==i_age_grps,'1','0')
    eval(parse(text=str_c("df$a",i_age_grps,"=if_else(df$age_grp==i_age_grps,'1','0') %>% factor")))
  }
  ind_race=if_else(ref_race=='Black','White','Black')
  if(ind_race=='Black'){
    df$rBlack=if_else(df$race=="Black",'1','0') %>% factor
  }else{
    df$rWhite=if_else(df$race=="White",'1','0') %>% factor
  }
  
  AllAgeTerm=str_c(str_c("a",age_grps),collapse="+")#works even if length(age_grps) is 1
  AllPerTerm=str_c(str_c("p",periods ),collapse="+")
  AllRacTerm=str_c("r",ind_race)
  AllAgePerRacTerm=str_c(AllAgeTerm,"+",AllPerTerm,"+",AllRacTerm)
  attr(df,"AllAgePerRacTerm")=AllAgePerRacTerm
  attr(df,"AllAgeTerm")=AllAgeTerm
  attr(df,"AllPerTerm")=AllPerTerm
  attr(df,"AllRacTerm")=AllRacTerm
  df
}

###########
#misc
###########
turn_irr_to_str=function(irr){
  out=as.character(round(irr,2))
  #if out[i]='1.1'
  out=if_else(str_length(out)==3,str_c(out,'0'),out)
  #if out[i]='1'
  if_else(str_length(out)==1,str_c(out,'.00'),out)
}

round_n_2_str=function(x,digits=2){
  #AIM: round to given digits and turn to string
  #NOTE: the function is vectorized
  #
  #ARGUMENTS:
  #x is input vector may contain '-' sign
  #
  
  x1=round(x,digits=digits)
  sign_x=if_else(str_detect(x1,'-'),'-','')
  x2=str_replace(x1,'-','')
  if(digits==2){
    x2=if_else(str_length(x2)==3,str_c(x2,'0'),x2)
    x2=if_else(str_length(x2)==1,str_c(x2,'.00'),x2)
  }
  str_c(sign_x,x2)
}

recode_cov_nam_2_fct=function(covariate_name,inclu_line_break=T){
  if(inclu_line_break==T){
    recode_factor(covariate_name,
                  'cs_m.brfss_amc_lag2'='Smoking (10yL)','cs_m.brfss_amc_cum10y5ylag2'='Smoking (C)',
                  'cs_m.nhis_amc_lag2'='Smoking (10yL)','cs_m.nhis_amc_cum10y5ylag2'='Smoking (C)','cs_m.nhanes.fill_amc_cum10y5ylag2'='Smoking (C)',
                  'adj_avg_alc_m.brfss_amc_lag2'='Alcohol (10yL)','adj_avg_alc_m.brfss_amc_cum10y5ylag2'='Alcohol (C)',
                  'adj_avg_alc_m.fill_amc_cum10y5ylag2'='Alcohol (C)',
                  'drinkhea.m2w2.brfss_amc_lag2'='Heavy alcohol\nm2w2 (10yL)','drinkhea.m2w2.brfss_amc_cum10y5ylag2'='Heavy alcohol\nm2w2 (C)',
                  'drinkhea.m2w1.brfss_amc_lag2'='Heavy alcohol\nm2w1 (10yL)','drinkhea.m2w1.brfss_amc_cum10y5ylag2'='Heavy alcohol\nm2w1 (C)',
                  'drinkhea.m2w2.fill_amc_cum10y5ylag2'='Heavy alcohol\nm2w2 (C)','drinkhea.m2w1.fill_amc_cum10y5ylag2'='Heavy alcohol\nm2w1 (C)',
                  'drinkbin.brfss_amc_lag2'='Binge alcohol (10yL)','drinkbin.brfss_amc_cum10y5ylag2'='Binge alcohol (C)',
                  'calcium_m.fill_amc_cum10y5ylag2'='Calcium (C)','fibe_m.fill_amc_cum10y5ylag2'='Fiber (C)','wholegrain_m.nhanes_amc_cum10y5ylag2'='Wholegrain (C)',
                  'dairy_m.nhanes_amc_cum10y5ylag2'='Dairy (C)','rednproc_m.nhanes_amc_cum10y5ylag2'='Red or\nprocessed\nmeat (C)',
                  'overweight_up_m.brfss_amc_no_lag'='Overweight or\nobesity (N)','overweight_up_m_amc_no_lag'='Overweight or\nobesity (N)',
                  'overweight_up_m.nhis_amc_no_lag'='Overweight or\nobesity (N)',
                  'overweight_up_m.brfss.fill_amc_no_lag'='Overweight or\nobesity (N)',
                  'overweight_up_m.nhis_amc_lag52'='Overweight or\nobesity (5yL)','overweight_up_m_amc_lag52'='Overweight or\nobesity (5yL)',
                  'overweight_up_m.brfss_amc_lag52'='Overweight or\nobesity (5yL)',
                  'overweight_up_m.brfss.fill_amc_lag52'='Overweight or\nobesity (5yL)',
                  'overweight_up_m.nhis_amc_lag2'='Overweight or\nobesity (10yL)',
                  'overweight_up_m.brfss.fill_amc_cum10y5ylag2'='Overweight or\nobesity (C)','overweight_up_m.fill_amc_cum10y5ylag2'='Overweight or\nobesity (C)',
                  'overweight_up_m.nhis_amc_cum10y5ylag2'='Overweight or\nobesity (C)',
                  'obese_m.brfss_amc_no_lag'='Obesity (N)','obese_m_amc_no_lag'='Obesity (N)','obese_m.nhis_amc_no_lag'='Obesity (N)',
                  'obese_m.brfss.fill_amc_no_lag'='Obesity (N)',
                  'obese_m.nhis_amc_lag52'='Obesity (5yL)','obese_m_amc_lag52'='Obesity (5yL)',
                  'obese_m.brfss_amc_lag52'='Obesity (5yL)',
                  'obese_m.brfss.fill_amc_lag52'='Obesity (5yL)',
                  'obese_m.nhis_amc_lag2'='Obesity (10yL)',
                  'obese_m.brfss.fill_amc_cum10y5ylag2'='Obesity (C)','obese_m.fill_amc_cum10y5ylag2'='Obesity (C)','obese_m.nhis_amc_cum10y5ylag2'='Obesity (C)',
                  'hbp_m.brfss_amc_no_lag'='Hypertension (N)','hbp_m.nhis_amc_no_lag'='Hypertension (N)',
                  'hbp_m.nhanes_amc_no_lag'='Hypertension (N)','hbp_m.nhanes.fill_amc_no_lag'='Hypertension (N)',
                  'hbp_m_no_preg.brfss_amc_no_lag'='Hypertension (N)',
                  'hbp_m.brfss_amc_lag52'='Hypertension (5yL)','hbp_m.nhis_amc_lag52'='Hypertension (5yL)','hbp_m.nhanes.fill_amc_lag52'='Hypertension (5yL)',
                  'hbp_m_no_preg.brfss_amc_lag52'='Hypertension (5yL)',
                  'hbp_m.brfss_amc_lag2'='Hypertension (10yL)','hbp_m.nhis_amc_lag2'='Hypertension (10yL)','hbp_m.nhanes.fill_amc_lag2'='Hypertension (10yL)',
                  'hbp_m_no_preg.brfss_amc_lag2'='Hypertension (10yL)',
                  'hbp_m.brfss_amc_cum10y5ylag2'='Hypertension (C)','hbp_m.nhis_amc_cum10y5ylag2'='Hypertension (C)',
                  'hbp_m.nhanes.fill_amc_cum10y5ylag2'='Hypertension (C)',
                  'hbp_m_no_preg.brfss_amc_cum10y5ylag2'='Hypertension (C)',
                  'dia_m.brfss_amc_no_lag'='Diabetes (N)','dia_m.nhis_amc_no_lag'='Diabetes (N)',
                  'dia_m.nhanes_amc_no_lag'='Diabetes (N)',
                  'dia_m.brfss.fill_amc_no_lag'='Diabetes (N)','dia_m.nhanes.fill_amc_no_lag'='Diabetes (N)',
                  'dia_m.nhanes_amc_lag52'='Diabetes (5yL)','dia_m.nhis_amc_lag52'='Diabetes (5yL)',
                  'dia_m.brfss.fill_amc_lag52'='Diabetes (5yL)','dia_m.nhanes.fill_amc_lag52'='Diabetes (5yL)',
                  'dia_m.nhis_amc_lag2'='Diabetes (10yL)','dia_m.nhanes.fill_amc_lag2'='Diabetes (10yL)',
                  'dia_m.nhis_amc_cum10y5ylag2'='Diabetes (C)',
                  'dia_m.brfss.fill_amc_cum10y5ylag2'='Diabetes (C)','dia_m.nhanes.fill_amc_cum10y5ylag2'='Diabetes (C)',
                  'thy_m_amc_no_lag'='Thyroid disease (N)','thy_m_amc_lag52'='Thyroid disease (5yL)',
                  'thy_m.fill_amc_no_lag'='Thyroid disease (N)','thy_m.fill_amc_lag52'='Thyroid disease (5yL)',
                  'thy_m.fill_amc_lag2'='Thyroid disease (10yL)','thy_m.fill_amc_cum10y5ylag2'='Thyroid disease (C)')
  }else{
    recode_factor(covariate_name,
                  'cs_m.brfss_amc_lag2'='Smoking (10yL)','cs_m.brfss_amc_cum10y5ylag2'='Smoking (C)',
                  'cs_m.nhis_amc_lag2'='Smoking (10yL)','cs_m.nhis_amc_cum10y5ylag2'='Smoking (C)','cs_m.nhanes.fill_amc_cum10y5ylag2'='Smoking (C)',
                  'adj_avg_alc_m.brfss_amc_lag2'='Alcohol (10yL)','adj_avg_alc_m.brfss_amc_cum10y5ylag2'='Alcohol (C)',
                  'adj_avg_alc_m.fill_amc_cum10y5ylag2'='Alcohol (C)',
                  'drinkhea.m2w2.brfss_amc_lag2'='Heavy alcohol m2w2 (10yL)','drinkhea.m2w2.brfss_amc_cum10y5ylag2'='Heavy alcohol m2w2 (C)',
                  'drinkhea.m2w1.brfss_amc_lag2'='Heavy alcohol m2w1 (10yL)','drinkhea.m2w1.brfss_amc_cum10y5ylag2'='Heavy alcohol m2w1 (C)',
                  'drinkhea.m2w2.fill_amc_cum10y5ylag2'='Heavy alcohol m2w2 (C)','drinkhea.m2w1.fill_amc_cum10y5ylag2'='Heavy alcohol m2w1 (C)',
                  'drinkbin.brfss_amc_lag2'='Binge alcohol (10yL)','drinkbin.brfss_amc_cum10y5ylag2'='Binge alcohol (C)',
                  'calcium_m.fill_amc_cum10y5ylag2'='Calcium (C)','fibe_m.fill_amc_cum10y5ylag2'='Fiber (C)','wholegrain_m.nhanes_amc_cum10y5ylag2'='Wholegrain (C)',
                  'dairy_m.nhanes_amc_cum10y5ylag2'='Dairy (C)','rednproc_m.nhanes_amc_cum10y5ylag2'='Red or processed meat (C)',
                  'overweight_up_m.brfss_amc_no_lag'='Overweight or obesity (N)','overweight_up_m_amc_no_lag'='Overweight or obesity (N)',
                  'overweight_up_m.nhis_amc_no_lag'='Overweight or obesity (N)',
                  'overweight_up_m.brfss.fill_amc_no_lag'='Overweight or obesity (N)',
                  'overweight_up_m.nhis_amc_lag52'='Overweight or obesity (5yL)','overweight_up_m_amc_lag52'='Overweight or obesity (5yL)',
                  'overweight_up_m.brfss_amc_lag52'='Overweight or obesity (5yL)',
                  'overweight_up_m.brfss.fill_amc_lag52'='Overweight or obesity (5yL)',
                  'overweight_up_m.nhis_amc_lag2'='Overweight or obesity (10yL)',
                  'overweight_up_m.brfss.fill_amc_cum10y5ylag2'='Overweight or obesity (C)','overweight_up_m.fill_amc_cum10y5ylag2'='Overweight or obesity (C)',
                  'overweight_up_m.nhis_amc_cum10y5ylag2'='Overweight or obesity (C)',
                  'obese_m.brfss_amc_no_lag'='Obesity (N)','obese_m_amc_no_lag'='Obesity (N)','obese_m.nhis_amc_no_lag'='Obesity (N)',
                  'obese_m.brfss.fill_amc_no_lag'='Obesity (N)',
                  'obese_m.nhis_amc_lag52'='Obesity (5yL)','obese_m_amc_lag52'='Obesity (5yL)',
                  'obese_m.brfss_amc_lag52'='Obesity (5yL)',
                  'obese_m.brfss.fill_amc_lag52'='Obesity (5yL)',
                  'obese_m.nhis_amc_lag2'='Obesity (10yL)',
                  'obese_m.brfss.fill_amc_cum10y5ylag2'='Obesity (C)','obese_m.fill_amc_cum10y5ylag2'='Obesity (C)','obese_m.nhis_amc_cum10y5ylag2'='Obesity (C)',
                  'hbp_m.brfss_amc_no_lag'='Hypertension (N)','hbp_m.nhis_amc_no_lag'='Hypertension (N)',
                  'hbp_m.nhanes_amc_no_lag'='Hypertension (N)','hbp_m.nhanes.fill_amc_no_lag'='Hypertension (N)',
                  'hbp_m_no_preg.brfss_amc_no_lag'='Hypertension (N)',
                  'hbp_m.brfss_amc_lag52'='Hypertension (5yL)','hbp_m.nhis_amc_lag52'='Hypertension (5yL)','hbp_m.nhanes.fill_amc_lag52'='Hypertension (5yL)',
                  'hbp_m_no_preg.brfss_amc_lag52'='Hypertension (5yL)',
                  'hbp_m.brfss_amc_lag2'='Hypertension (10yL)','hbp_m.nhis_amc_lag2'='Hypertension (10yL)','hbp_m.nhanes.fill_amc_lag2'='Hypertension (10yL)',
                  'hbp_m_no_preg.brfss_amc_lag2'='Hypertension (10yL)',
                  'hbp_m.brfss_amc_cum10y5ylag2'='Hypertension (C)','hbp_m.nhis_amc_cum10y5ylag2'='Hypertension (C)',
                  'hbp_m.nhanes.fill_amc_cum10y5ylag2'='Hypertension (C)',
                  'hbp_m_no_preg.brfss_amc_cum10y5ylag2'='Hypertension (C)',
                  'dia_m.brfss_amc_no_lag'='Diabetes (N)','dia_m.nhis_amc_no_lag'='Diabetes (N)',
                  'dia_m.nhanes_amc_no_lag'='Diabetes (N)',
                  'dia_m.brfss.fill_amc_no_lag'='Diabetes (N)','dia_m.nhanes.fill_amc_no_lag'='Diabetes (N)',
                  'dia_m.nhanes_amc_lag52'='Diabetes (5yL)','dia_m.nhis_amc_lag52'='Diabetes (5yL)',
                  'dia_m.brfss.fill_amc_lag52'='Diabetes (5yL)','dia_m.nhanes.fill_amc_lag52'='Diabetes (5yL)',
                  'dia_m.nhis_amc_lag2'='Diabetes (10yL)','dia_m.nhanes.fill_amc_lag2'='Diabetes (10yL)',
                  'dia_m.nhis_amc_cum10y5ylag2'='Diabetes (C)',
                  'dia_m.brfss.fill_amc_cum10y5ylag2'='Diabetes (C)','dia_m.nhanes.fill_amc_cum10y5ylag2'='Diabetes (C)',
                  'thy_m_amc_no_lag'='Thyroid disease (N)','thy_m_amc_lag52'='Thyroid disease (5yL)',
                  'thy_m.fill_amc_no_lag'='Thyroid disease (N)','thy_m.fill_amc_lag52'='Thyroid disease (5yL)',
                  'thy_m.fill_amc_lag2'='Thyroid disease (10yL)','thy_m.fill_amc_cum10y5ylag2'='Thyroid disease (C)')
  }
}

recode_interact_cov_nam=function(covariate_name,inclu_line_break=T){
  #########
  #prepare
  cov_nam_print=c('Smoking (10yL)','Smoking (C)',
                  'Alcohol daily # (10yL)','Alcohol daily # (C)',
                  'Alcohol heavy (10yL)','Alcohol heavy (C)','Alcohol binge (10yL)','Alcohol binge (C)',
                  'Calcium (C)','Fiber (C)','Wholegrain (C)','Dairy (C)','Red & processed meat (C)',
                  'Overweight and obesity (N)','Overweight and obesity (5yL)','Overweight and obesity (10yL)','Overweight and obesity (C)',
                  'Obesity (N)','Obesity (5yL)','Obesity (10yL)','Obesity (C)',
                  'Hypertension (N)','Hypertension (5yL)','Hypertension (10yL)','Hypertension (C)',
                  'Diabetes (N)','Diabetes (5yL)','Diabetes (10yL)','Diabetes (C)',
                  'Thyroid diseases (N)','Thyroid diseases (5yL)','Thyroid diseases (10yL)','Thyroid diseases (C)')
  sep_use=if(inclu_line_break==T) '-' else '*'
  #########
  #main
  covariate_name %>% str_replace_all(.,'cs_m.brfss_amc_lag2','Smoking (10yL)') %>% str_replace_all(.,'cs_m.brfss_amc_cum10y5ylag2','Smoking (C)') %>%
    str_replace_all(.,'cs_m.nhis_amc_lag2','Smoking (10yL)') %>% str_replace_all(.,'cs_m.nhis_amc_cum10y5ylag2','Smoking (C)') %>%
    str_replace_all(.,'cs_m.nhanes.fill_amc_cum10y5ylag2','Smoking (C)') %>%
    str_replace_all(.,'adj_avg_alc_m.brfss_amc_lag2','Alcohol daily # (10yL)') %>% str_replace_all(.,'adj_avg_alc_m.brfss_amc_cum10y5ylag2','Alcohol daily # (C)') %>%
    str_replace_all(.,'adj_avg_alc_m.fill_amc_cum10y5ylag2','Alcohol daily # (C)') %>%
    str_replace_all(.,'drinkhea.m2w2.brfss_amc_lag2','Alcohol heavy (10yL)') %>% 
    str_replace_all(.,'drinkhea.m2w2.brfss_amc_cum10y5ylag2','Alcohol heavy (C)') %>%
    str_replace_all(.,'drinkhea.m2w2.fill_amc_cum10y5ylag2','Alcohol heavy (C)') %>%
    str_replace_all(.,'drinkbin.brfss_amc_lag2','Alcohol binge (10yL)') %>% str_replace_all(.,'drinkbin.brfss_amc_cum10y5ylag2','Alcohol binge (C)') %>%
    str_replace_all(.,'calcium_m.fill_amc_cum10y5ylag2','Calcium (C)') %>% str_replace_all(.,'fibe_m.fill_amc_cum10y5ylag2','Fiber (C)') %>%
    str_replace_all(.,'wholegrain_m.nhanes_amc_cum10y5ylag2','Wholegrain (C)') %>% str_replace_all(.,'dairy_m.nhanes_amc_cum10y5ylag2','Dairy (C)') %>%
    str_replace_all(.,'rednproc_m.nhanes_amc_cum10y5ylag2','Red & processed meat (C)') %>% 
    str_replace_all(.,'overweight_up_m.brfss_amc_no_lag','Overweight and obesity (N)') %>%
    str_replace_all(.,'overweight_up_m_amc_no_lag','Overweight and obesity (N)') %>% 
    str_replace_all(.,'overweight_up_m.nhis_amc_no_lag','Overweight and obesity (N)') %>%
    str_replace_all(.,'overweight_up_m.brfss.fill_amc_no_lag','Overweight and obesity (N)') %>%
    str_replace_all(.,'overweight_up_m.brfss_amc_lag52','Overweight and obesity (5yL)') %>%
    str_replace_all(.,'overweight_up_m_amc_lag52','Overweight and obesity (5yL)') %>% 
    str_replace_all(.,'overweight_up_m.nhis_amc_lag52','Overweight and obesity (5yL)') %>%
    str_replace_all(.,'overweight_up_m.brfss.fill_amc_lag52','Overweight and obesity (5yL)') %>%
    str_replace_all(.,'overweight_up_m.nhis_amc_lag2','Overweight and obesity (10yL)') %>%
    str_replace_all(.,'overweight_up_m.brfss.fill_amc_cum10y5ylag2','Overweight and obesity (C)') %>%
    str_replace_all(.,'overweight_up_m.fill_amc_cum10y5ylag2','Overweight and obesity (C)') %>% 
    str_replace_all(.,'overweight_up_m.nhis_amc_cum10y5ylag2','Overweight and obesity (C)') %>%
    str_replace_all(.,'obese_m.brfss_amc_no_lag','Obesity (N)') %>% str_replace_all(.,'obese_m_amc_no_lag','Obesity (N)') %>% 
    str_replace_all(.,'obese_m.nhis_amc_no_lag','Obesity (N)') %>%
    str_replace_all(.,'obese_m.brfss.fill_amc_no_lag','Obesity (N)') %>%
    str_replace_all(.,'obese_m.brfss_amc_lag52','Obesity (5yL)') %>% str_replace_all(.,'obese_m_amc_lag52','Obesity (5yL)') %>% 
    str_replace_all(.,'obese_m.nhis_amc_lag52','Obesity (5yL)') %>% 
    str_replace_all(.,'obese_m.brfss.fill_amc_lag52','Obesity (5yL)') %>% 
    str_replace_all(.,'obese_m.nhis_amc_lag2','Obesity (10yL)') %>%
    str_replace_all(.,'obese_m.brfss.fill_amc_cum10y5ylag2','Obesity (C)') %>% str_replace_all(.,'obese_m.fill_amc_cum10y5ylag2','Obesity (C)') %>% 
    str_replace_all(.,'obese_m.nhis_amc_cum10y5ylag2','Obesity (C)') %>%
    str_replace_all(.,'hbp_m.brfss_amc_no_lag','Hypertension (N)') %>% str_replace_all(.,'hbp_m.nhis_amc_no_lag','Hypertension (N)') %>%
    str_replace_all(.,'hbp_m.nhanes_amc_no_lag','Hypertension (N)') %>% str_replace_all(.,'hbp_m.nhanes.fill_amc_no_lag','Hypertension (N)') %>% 
    str_replace_all(.,'hbp_m_no_preg.brfss_amc_no_lag','Hypertension (N)') %>% 
    str_replace_all(.,'hbp_m.nhanes.fill_amc_lag52','Hypertension (5yL)') %>% 
    str_replace_all(.,'hbp_m.brfss_amc_lag52','Hypertension (5yL)') %>% str_replace_all(.,'hbp_m.nhis_amc_lag52','Hypertension (5yL)') %>%
    str_replace_all(.,'hbp_m_no_preg.brfss_amc_lag52','Hypertension (5yL)') %>%
    str_replace_all(.,'hbp_m.brfss_amc_lag2','Hypertension (10yL)') %>% str_replace_all(.,'hbp_m.nhis_amc_lag2','Hypertension (10yL)') %>%
    str_replace_all(.,'hbp_m_no_preg.brfss_amc_lag2','Hypertension (10yL)') %>%
    str_replace_all(.,'hbp_m.brfss_amc_cum10y5ylag2','Hypertension (C)') %>% str_replace_all(.,'hbp_m.nhis_amc_cum10y5ylag2','Hypertension (C)') %>%
    str_replace_all(.,'hbp_m_no_preg.brfss_amc_cum10y5ylag2','Hypertension (C)') %>%
    str_replace_all(.,'hbp_m.nhanes.fill_amc_cum10y5ylag2','Hypertension (C)') %>%
    str_replace_all(.,'dia_m.brfss_amc_no_lag','Diabetes (N)') %>% str_replace_all(.,'dia_m.nhanes_amc_no_lag','Diabetes (N)')  %>% 
    str_replace_all(.,'dia_m.nhis_amc_no_lag','Diabetes (N)')  %>% str_replace_all(.,'dia_m.brfss.fill_amc_no_lag','Diabetes (N)') %>%
    str_replace_all(.,'dia_m.nhanes_amc_lag52','Diabetes (5yL)')  %>% str_replace_all(.,'dia_m.nhis_amc_lag52','Diabetes (5yL)')  %>% 
    str_replace_all(.,'dia_m.brfss.fill_amc_lag52','Diabetes (5yL)')  %>% 
    str_replace_all(.,'dia_m.nhis_amc_lag2','Diabetes (10yL)')  %>% 
    str_replace_all(.,'dia_m.nhanes.fill_amc_no_lag','Diabetes (N)')  %>% str_replace_all(.,'dia_m.nhanes.fill_amc_lag52','Diabetes (5yL)')  %>% 
    str_replace_all(.,'dia_m.brfss.fill_amc_cum10y5ylag2','Diabetes (C)')  %>% str_replace_all(.,'dia_m.nhanes.fill_amc_cum10y5ylag2','Diabetes (C)')  %>% 
    str_replace_all(.,'dia_m.nhis_amc_cum10y5ylag2','Diabetes (C)')  %>% 
    str_replace_all(.,'thy_m_amc_no_lag','Thyroid diseases (N)') %>% str_replace_all(.,'thy_m_amc_lag52','Thyroid diseases (5yL)') %>%
    str_replace_all(.,'thy_m.fill_amc_no_lag','Thyroid diseases (N)') %>% str_replace_all(.,'thy_m.fill_amc_lag52','Thyroid diseases (5yL)') %>%
    str_replace_all(.,'thy_m.fill_amc_cum10y5ylag2','Thyroid diseases (C)') %>%
    str_replace_all(.,'rBlack1','Black race') %>% str_replace_all(.,'rWhite1','White race') %>% 
    str_replace_all(.,'a61','Age30-34') %>% str_replace_all(.,'a71','Age35-39') %>% str_replace_all(.,'a81','Age40-44') %>% str_replace_all(.,'a91','Age45-49') %>%
    str_replace_all(.,'rBlack','Black race') %>% str_replace_all(.,'rWhite','White race') %>% 
    str_replace_all(.,'a6','Age30-34') %>% str_replace_all(.,'a7','Age35-39') %>% str_replace_all(.,'a8','Age40-44') %>% str_replace_all(.,'a9','Age45-49') %>% {
      #the below chunk ensure seq of variables are right, by creating levels of factors
      x=.
      
      #get levels for exposure_exposure_InterTerm
      exposure_exposure_InterTerm=expand.grid(cov_nam_print,cov_nam_print) %>% as_tibble %>% 
        transmute(exposure_exposure_InterTerm=str_c(Var1,Var2,sep=sep_use)) %>% .$exposure_exposure_InterTerm
      
      #get levels for race_exposure_InterTerm
      if(any(str_detect(x,pattern='Black race'))){#reference grp was set to White
        race_exposue_InterTerm=str_c('Black race',cov_nam_print,sep=sep_use)
      }else{#ref grp was set to Black
        race_exposue_InterTerm=str_c('White race',cov_nam_print,sep=sep_use)
      }
      
      #get levels for age_exposure_InterTerm
      if(any(str_detect(x,pattern='Age45-49'))){#ref grp for age set to 1, ie, Age30-34
        age_exposure_InterTerm=character(3*length(cov_nam_print))
        for (cov_nam_print_i in cov_nam_print) {
          i_strt=(which(cov_nam_print==cov_nam_print_i)-1)*3+1
          i_end=i_strt+2
          age_exposure_InterTerm[i_strt:i_end]=str_c(c('Age35-39','Age40-44','Age45-49'),cov_nam_print_i,sep=sep_use)
        }
      }else{#ref grp for age set to 4, ie, Age45-49
        age_exposure_InterTerm=character(3*length(cov_nam_print))
        for (cov_nam_print_i in cov_nam_print) {
          i_strt=(which(cov_nam_print==cov_nam_print_i)-1)*3+1
          i_end=i_strt+2
          age_exposure_InterTerm[i_strt:i_end]=str_c(c('Age40-44','Age35-39','Age30-34'),cov_nam_print_i,sep=sep_use)
        }
      }
      
      #use these levels to create factor
      x %>% str_replace_all(.,':',sep_use) %>% factor(.,levels=c(cov_nam_print,exposure_exposure_InterTerm,race_exposue_InterTerm,age_exposure_InterTerm))
      
    }
  
}

get_cutoff=function(vTag){
  if(vTag %in% c('v1a','v1b','v2a','v3a')){
    '0.6'
  }else if(vTag %in% c('v1c')){
    '0.7'
  }else if(vTag %in% c('v1','v2','v3')){
    0.7
  }
}

get_RfSource=function(vTag){
  if(vTag %in% c('v1a','v1b','v1c')){
    c('Most'='BRFSS')
  }else if(vTag %in% c('v2a')){
    c('Most'='NHIS')
  }else if(vTag %in% c('v3a')){
    c('Most'='NHANES')
  }
}

get_Var2Drop=function(vTag,gender){
  if(vTag %in% c('v1a')){
    c('drinkhea.m2w1.brfss_amc_lag2','drinkhea.m2w1.brfss_amc_cum10y5ylag2')
  }else if(vTag %in% c('v1b','v1c','v2a','v3a')){
    c('drinkhea.m2w1.brfss_amc_lag2','drinkhea.m2w1.brfss_amc_cum10y5ylag2','drinkbin.brfss_amc_lag2','drinkbin.brfss_amc_cum10y5ylag2')
  }
}

get_ref_grp=function(vTag){
  if(vTag %in% c('v1','v1a','v1b','v1c','v2','v2a','v3','v3a')){
    list(age=1,race='White')
  }
}

rank_mod=function(GOF){
  #AIM
  #need to rank models in result of mvm when each model can take more than 1 row (b/c of multiple exposures)
  #
  #create a vectorized function that rank the models within each main effect combo; the same BIC should be together
  #eg, take in BIC=c(100,100,99,99,99,300,300,300,300) and return c(2,2,1,1,1,3,3,3,3)
  #need unique vector
  #count num of duplicates
  #get rank of unique vector
  #expand the rank given duplicates
  
  uniq_vec=unique(GOF)
  num_of_dup=GOF %>% as.character %>% factor %>% fct_inorder %>% table %>% {x=.;attributes(x)=NULL;x}
  rank_uniq_vec=rank(uniq_vec)
  out=list()
  for (i in 1:length(uniq_vec))  out[[i]]=rep(rank_uniq_vec[i],num_of_dup[i])
  unlist(out)
}


##########
#find nested models
##########
find_nest_mod=function(df){
  #VALUES
  ##preserve the original df & all it's attr
  ##add a col: parsim_mod: T/F, indicating whether a model (on >=1 rows) is the best among a nest of models (a parsimonious model)
  ##add as attr a list, with each element/vector being the cov_nam_raw of those parsimonious models (ie, no model ranked higher than it, and also nested in it)
  #
  #ARGUMENT
  #df is returned by proc_mvm_res/proc_mvm_interact_res
  
  #save attribute; otherwise it will be lost on cluster
  gender=attributes(df)$gender
  RankBy=attributes(df)$RankBy
  cutoff=attributes(df)$cutoff
  RfSource=attributes(df)$RfSource
  vTag=attributes(df)$vTag
  ref_grp=attributes(df)$ref_grp
  
  #preserve df, get a copy of it for playing, ie, df_sub
  #go through each model of df_sub (each model contain >=1 rows)
  #use a count variable, as long as count<the seq num of the last model, keep running
  #get cov_nam for current model
  #group_by models, check whether a model is a superset of the current model
  #filter df_sub to get 1) the models whose rank are equal/higher than count OR 2) the models whose rank are lower and is a superset
  df_sub=df
  cnt=1
  while (cnt<length(unique(df_sub$model_rank))) {
    model_rank_i=unique(df_sub$model_rank)[cnt]
    cov_nam_i=df_sub[df_sub$model_rank==model_rank_i,]$cov_nam_raw
    df_sub=df_sub %>% group_by(model_rank) %>% mutate(is_superset=all(cov_nam_i %in% cov_nam_raw)) %>% ungroup %>% 
      filter(model_rank<=model_rank_i | (model_rank>model_rank_i & !is_superset))
    cnt=cnt+1
  }
  
  #turn df_sub to the list I want
  model_rank_keep=unique(df_sub$model_rank)#not consecutive
  df=suppressMessages(left_join(df,tibble(model_rank=model_rank_keep,parsim_mod=1),by='model_rank'))
  df$parsim_mod=if_else(is.na(df$parsim_mod),F,T)
  
  mod_ls=list()
  for (model_rank_keep_i in model_rank_keep) {
    mod_ls[[which(model_rank_keep==model_rank_keep_i)]]=df_sub %>% filter(model_rank==model_rank_keep_i) %>% .$cov_nam_raw
  }
  
  #######
  #output
  #######
  attr(df,'parsim_mod_ls')=mod_ls
  attr(df,"gender")=gender
  attr(df,'RankBy')=RankBy
  attr(df,'cutoff')=cutoff
  attr(df,'RfSource')=RfSource
  attr(df,'vTag')=vTag
  attr(df,'ref_grp')=ref_grp
  
  df
}


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