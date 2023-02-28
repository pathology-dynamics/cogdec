%macro quantreg(); 

%let dep_list =
CFI_total
FCSRT_freerecall 
FCSRT_cuedrecall 
FCSRT_total 
RCFT_copy
RCFT_immedrecall
RCFT_delayrecall 
RCFT_recogtotal
RAVLT_a5_score
RAVLT_immed
RAVLT_delay 
MOCA_total
nsfd_longest
nsbk_longest
trailA
trailB
VFPT_FAStotal
animal
MINT
JoLO_total
; 

%let ind_list =
ctr_ttau 
ctr_ptau
ctr_AB42 
ctr_tTau_AB42
ctr_pTau_AB42
;

%local i next_dep;
%local j next_ind;
%do i=1 %to %sysfunc(countw(&dep_list));
%do j=1 %to %sysfunc(countw(&ind_list));
   %let next_dep = %scan(&dep_list, &i);
   %let next_ind = %scan(&ind_list, &j);
       
        proc quantreg data=Elecsys.EHBS_CSF_Cog_CV_BL_FU1 algorithm=simplex ci=resampling namelen=32;
		where race in('Caucasian or White','Black or African American');
        model 
        &next_dep = &next_ind ASYM_AD2 ctr_age5 female race ctr_edu &next_ind*&next_ind &next_ind*ASYM_AD2 ctr_age5*ASYM_AD2 female*ASYM_AD2 ctr_edu*ASYM_AD2/ quantile=0.5 seed=1268;
		ods output ParameterEstimates=Est_quant_curv
                   ModelInfo=ModelInfo_quant_curv;
run;
        data ModelInfo_quant_curv;
		   set ModelInfo_quant_curv(where=(label1='Dependent Variable'));
		   rename cvalue1=dependent;
		   keep cvalue1;
        run; 

		data Est_quant_curv;
		   merge ModelInfo_quant_curv Est_quant_curv(where=(parameter^='Intercept'));
		run;

%end;
%end;
%mend;
