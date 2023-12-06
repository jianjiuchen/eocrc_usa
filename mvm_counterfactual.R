dir=""
source(str_c(dir,'crc_mvm_other.R'))
source(str_c(dir,'utils.R'))

##############
#get data for men and women
main_fun=function(vTag='v1b',gender,boot_n=5000){
  load(str_c(dir,'OutComb_',vTag,'_',gender,'_interact_5FoldCV.RData'))
  cov_nam_raw=AllModels %>% filter(model_rank==first(model_rank)) %>% .$cov_nam_raw
  print(str_c('For ',gender,', the raw covariate names are:'))
  print(cov_nam_raw)
  crc_mvm4adj_crc(gender=gender,cov_nam_raw=cov_nam_raw,boot_n=boot_n,dir_utils=str_c(dir,'utils.R'),dir=dir)
}
data_m=main_fun(gender='Male')
data_w=main_fun(gender='Female')

##############
#calculate the % of EOCRC averted (Std error)
#men 'overweight_up_m.brfss.fill_amc_lag52'
pct_m=data_m %>% 
  group_by(period,age_grp,race) %>% 
  summarise(ir_pred=mean(count_pred/pop),
            ir_pred_adj_overweight_up_m.brfss.fill_amc_lag52=mean(count_pred_adj_overweight_up_m.brfss.fill_amc_lag52/pop)) %>% ungroup %>%
  summarise(sum(ir_pred-ir_pred_adj_overweight_up_m.brfss.fill_amc_lag52)/sum(ir_pred)*100) %>% unname %>% unlist

des_mat_ls_m=list()
des_mat_ls_m[['true']]=attributes(data_m)$res_reg$des_mat
for (i in 1:nrow(des_mat_ls_m[['true']])) {
  des_mat_ls_m[['true']][i,]=str_c(as.character(des_mat_ls_m[['true']][i,]),
                                   str_c('x',1:7),
                                   sep='*')
}
fun_agg_row=function(x){
  tmp=str_c('exp(',str_c(x,collapse='+'),')')
  str_replace_all(tmp,pattern='\\+-',replacement='-')
}
crc_pred_m=apply(des_mat_ls_m[['true']],1,FUN=fun_agg_row) %>% str_c(.,collapse='+')
des_mat_ls_m[['overweight_up_m.brfss.fill_amc_lag52']]=attributes(data_m)$res_reg$des_mat_overweight_up_m.brfss.fill_amc_lag52
for (i in 1:nrow(des_mat_ls_m[['overweight_up_m.brfss.fill_amc_lag52']])) {
  des_mat_ls_m[['overweight_up_m.brfss.fill_amc_lag52']][i,]=str_c(as.character(des_mat_ls_m[['overweight_up_m.brfss.fill_amc_lag52']][i,]),
                                                                   str_c('x',1:7),sep='*')
}
crc_pred_counterf_m=apply(des_mat_ls_m[['overweight_up_m.brfss.fill_amc_lag52']],1,FUN=fun_agg_row) %>% str_c(.,collapse='+')


SE_m=msm::deltamethod(g=as.formula(str_c('~100*((',crc_pred_m,')-(',crc_pred_counterf_m,'))/(',crc_pred_m,')')),
                      mean=unname(coef(attributes(data_m)$res_reg$res_reg)),
                      cov=vcov(attributes(data_m)$res_reg$res_reg))

#women "adj_avg_alc_m.brfss_amc_cum10y5ylag2", "overweight_up_m.brfss.fill_amc_no_lag"
pct_w=data_w %>% 
  group_by(period,age_grp,race) %>% 
  summarise(ir_pred=mean(count_pred/pop),
            ir_pred_adj_adj_avg_alc_m.brfss_amc_cum10y5ylag2=mean(count_pred_adj_adj_avg_alc_m.brfss_amc_cum10y5ylag2/pop),
            ir_pred_adj_overweight_up_m.brfss.fill_amc_no_lag=mean(count_pred_adj_overweight_up_m.brfss.fill_amc_no_lag/pop)) %>% ungroup %>%
  summarise(alc_per=sum(ir_pred-ir_pred_adj_adj_avg_alc_m.brfss_amc_cum10y5ylag2)/sum(ir_pred)*100,
            ove_per=sum(ir_pred-ir_pred_adj_overweight_up_m.brfss.fill_amc_no_lag)/sum(ir_pred)*100) %>% unlist

des_mat_ls_w=list()
des_mat_ls_w[['true']]=attributes(data_w)$res_reg$des_mat
for (i in 1:nrow(des_mat_ls_w[['true']])) {
  des_mat_ls_w[['true']][i,]=str_c(as.character(des_mat_ls_w[['true']][i,]),
                                   str_c('x',1:length(coef(attributes(data_w)$res_reg$res_reg))),
                                   sep='*')
}
crc_pred_w=apply(des_mat_ls_w[['true']],1,FUN=fun_agg_row) %>% str_c(.,collapse='+')
des_mat_ls_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']]=attributes(data_w)$res_reg$des_mat_adj_avg_alc_m.brfss_amc_cum10y5ylag2
for (i in 1:nrow(des_mat_ls_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']])) {
  des_mat_ls_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']][i,]=str_c(as.character(des_mat_ls_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']][i,]),
                                                                   str_c('x',1:length(coef(attributes(data_w)$res_reg$res_reg))),sep='*')
}
crc_pred_counterf_w=list()
crc_pred_counterf_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']]=apply(des_mat_ls_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']],1,FUN=fun_agg_row) %>%str_c(.,collapse='+')
des_mat_ls_w[['overweight_up_m.brfss.fill_amc_no_lag']]=attributes(data_w)$res_reg$des_mat_overweight_up_m.brfss.fill_amc_no_lag
for (i in 1:nrow(des_mat_ls_w[['overweight_up_m.brfss.fill_amc_no_lag']])) {
  des_mat_ls_w[['overweight_up_m.brfss.fill_amc_no_lag']][i,]=str_c(as.character(des_mat_ls_w[['overweight_up_m.brfss.fill_amc_no_lag']][i,]),
                                                                    str_c('x',1:length(coef(attributes(data_w)$res_reg$res_reg))),sep='*')
}
crc_pred_counterf_w[['overweight_up_m.brfss.fill_amc_no_lag']]=apply(des_mat_ls_w[['overweight_up_m.brfss.fill_amc_no_lag']],1,FUN=fun_agg_row) %>%
  str_c(.,collapse='+')

SE_w=list()
SE_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']]=
  msm::deltamethod(g=as.formula(str_c('~100*((',crc_pred_w,')-(',crc_pred_counterf_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']],'))/(',crc_pred_w,')')),
                   mean=unname(coef(attributes(data_w)$res_reg$res_reg)),
                   cov=vcov(attributes(data_w)$res_reg$res_reg))
SE_w[['overweight_up_m.brfss.fill_amc_no_lag']]=
  msm::deltamethod(g=as.formula(str_c('~100*((',crc_pred_w,')-(',crc_pred_counterf_w[['overweight_up_m.brfss.fill_amc_no_lag']],'))/(',crc_pred_w,')')),
                   mean=unname(coef(attributes(data_w)$res_reg$res_reg)),
                   cov=vcov(attributes(data_w)$res_reg$res_reg))

#########
# result
#########
str_c('For men, increases in overweight and obesity prevalence contributed to ',round(pct_m),'% (SD: ',round(SE_m),'%) of EOCRC incidence during 1992-2016')
str_c('For women, increases in overweight and obesity prevalence contributed to ',round(pct_w['ove_per']),
      '% (SD: ',round(SE_w[['overweight_up_m.brfss.fill_amc_no_lag']]),'%) of EOCRC incidence during 1992-2016')
str_c('For women, increases in alcohol drinking contributed to ',round(pct_w['alc_per']),
      '% (SD: ',round(SE_w[['adj_avg_alc_m.brfss_amc_cum10y5ylag2']]),'%) of EOCRC incidence during 1992-2016')

