clear all
clc

disp('Computing the results...')

datasets=1:100; %%
p_thres=10.^linspace(-10,log10(0.05),100);
q_thres=0.001:0.001:0.1;
types = {'LD','HD','CYTO'};
a = 2;

% mean FDR
results_FDR_nostablenoambig=zeros(length(p_thres),length(types),length(dataset));
results_FDR_nostable=zeros(length(p_thres),length(types),length(dataset));
results_FDR_noambig=zeros(length(p_thres),length(types),length(dataset));
results_FDR_pcp=zeros(length(p_thres),length(types),length(dataset));
results_FDR_pcpnorobust=zeros(length(p_thres),length(types),length(dataset));
results_FDR_original=zeros(length(p_thres),length(types),length(dataset));
results_FDR_mmpc=zeros(length(p_thres),length(types),length(dataset));

% control bias
results_CB_nostablenoambig=zeros(length(p_thres),length(types),length(dataset));
results_CB_nostable=zeros(length(p_thres),length(types),length(dataset));
results_CB_noambig=zeros(length(p_thres),length(types),length(dataset));
results_CB_pcp=zeros(length(p_thres),length(types),length(dataset));
results_CB_pcpnorobust=zeros(length(p_thres),length(types),length(dataset));
results_CB_original=zeros(length(p_thres),length(types),length(dataset));
results_CB_mmpc=zeros(length(p_thres),length(types),length(dataset));

% estimation bias
results_EB_nostablenoambig=zeros(length(p_thres),length(types),length(dataset));
results_EB_nostable=zeros(length(p_thres),length(types),length(dataset));
results_EB_noambig=zeros(length(p_thres),length(types),length(dataset));
results_EB_pcp=zeros(length(p_thres),length(types),length(dataset));
results_EB_pcpnorobust=zeros(length(p_thres),length(types),length(dataset));
results_EB_original=zeros(length(p_thres),length(types),length(dataset));
results_EB_mmpc=zeros(length(p_thres),length(types),length(dataset));

load('CYTO_Dlingam.mat')
for d=1:length(datasets)
    disp(d)
    for t=1:length(types)
    
        %% PCp
        if strcmp(types{t}, 'LD')
            load(['PC_pcp_graphLD1_synth' num2str(d) '_' num2str(a) '.mat'])
            graph=B>0;
        elseif strcmp(types{t}, 'HD')
            load(['PC_pcp_graphHD1_synth' num2str(d) '_' num2str(a) '.mat'])
            graph=B>0;
        elseif strcmp(types{t}, 'CYTO')
            load(['PC_pcp_graphCYTO1_real' num2str(d) '_' num2str(a) '.mat'])
            graph=Bsig>0;
        end
        p=correct_p_values(cell_p,num_struc);
        for pt1=1:length(p_thres)
            results_FDR_pcp(pt1,t,d)=get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
            results_CB_pcp(pt1,t,d)=get_proc_bias2_learned(q_thres(pt1),pdag,graph,cell_p,num_struc);
            results_EB_pcp(pt1,t,d)=get_BY_FDR(p,p_thres(pt1))-get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
        end
        
        %% PC
        if strcmp(types{t}, 'LD')
            load(['PC_original_graphLD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'HD')
            load(['PC_original_graphHD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'CYTO')
            load(['PC_original_graphCYTO1_real' num2str(d) '_' num2str(a) '.mat'])
            graph=Bsig>0;
        end   
        p=correct_p_values(cell_p,num_struc);
        for pt1=1:length(p_thres)
            results_FDR_original(pt1,t,d)=get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
            results_CB_original(pt1,t,d)=get_proc_bias2_learned(q_thres(pt1),pdag,graph,cell_p,num_struc);
            results_EB_original(pt1,t,d)=get_BY_FDR(p,p_thres(pt1))-get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
        end
        
        %% PC-p without robust p-values
        if strcmp(types{t}, 'LD')
            load(['PC_pcpnorobust_graphLD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'HD')
            load(['PC_pcpnorobust_graphHD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'CYTO')
            load(['PC_pcpnorobust_graphCYTO1_real' num2str(d) '_' num2str(a) '.mat'])
            graph=Bsig>0;
        end       
        p=correct_p_values(cell_p,num_struc);
        for pt1=1:length(p_thres)
            results_FDR_pcpnorobust(pt1,t,d)=get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
            results_CB_pcpnorobust(pt1,t,d)=get_proc_bias2_learned(q_thres(pt1),pdag,graph,cell_p,num_struc);
            results_EB_pcpnorobust(pt1,t,d)=get_BY_FDR(p,p_thres(pt1))-get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
        end
        
        %% PC-p without stabilization
        if strcmp(types{t}, 'LD')
            load(['PC_nostable_graphLD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'HD')
            load(['PC_nostable_graphHD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'CYTO')
            load(['PC_nostable_graphCYTO1_real' num2str(d) '_' num2str(a) '.mat'])
            graph=Bsig>0;
        end       
        p=correct_p_values(cell_p,num_struc);
        for pt1=1:length(p_thres)
            results_FDR_nostable(pt1,t,d)=get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
            results_CB_nostable(pt1,t,d)=get_proc_bias2_learned(q_thres(pt1),pdag,graph,cell_p,num_struc);
            results_EB_nostable(pt1,t,d)=get_BY_FDR(p,p_thres(pt1))-get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
        end
        
        %% PC-p without ambiguation
        if strcmp(types{t}, 'LD')
            load(['PC_noambig_graphLD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'HD')
            load(['PC_noambig_graphHD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'CYTO')
            load(['PC_noambig_graphCYTO1_real' num2str(d) '_' num2str(a) '.mat'])
            graph=Bsig>0;
        end       
        p=correct_p_values(cell_p,num_struc);
        for pt1=1:length(p_thres)
            results_FDR_noambig(pt1,t,d)=get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
            results_CB_noambig(pt1,t,d)=get_proc_bias2_learned(q_thres(pt1),pdag,graph,cell_p,num_struc);
            results_EB_noambig(pt1,t,d)=get_BY_FDR(p,p_thres(pt1))-get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
        end
        
        %% PC-p without stabilization and ambiguation
        if strcmp(types{t}, 'LD')
            load(['PC_nostablenoambig_graphLD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'HD')
            load(['PC_nostablenoambig_graphHD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'CYTO')
            load(['PC_nostablenoambig_graphCYTO1_real' num2str(d) '_' num2str(a) '.mat'])
            graph=Bsig>0;
        end       
        p=correct_p_values(cell_p,num_struc);
        for pt1=1:length(p_thres)
            results_FDR_nostablenoambig(pt1,t,d)=get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
            results_CB_nostablenoambig(pt1,t,d)=get_proc_bias2_learned(q_thres(pt1),pdag,graph,cell_p,num_struc);
            results_EB_nostablenoambig(pt1,t,d)=get_BY_FDR(p,p_thres(pt1))-get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
        end

        %% MMPC
        if strcmp(types{t}, 'LD')
            load(['MMPC_graphLD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'HD')
            load(['MMPC_graphHD1_synth' num2str(d) '_' num2str(a) '.mat'])
        elseif strcmp(types{t}, 'CYTO')
            load(['MMPC_graphCYTO1_real' num2str(d) '_' num2str(a) '.mat'])
            graph=Bsig>0;
        end       
        p=correct_p_values(cell_p,num_struc);
        for pt1=1:length(p_thres)
            results_FDR_mmpc(pt1,t,d)=get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
            results_CB_mmpc(pt1,t,d)=get_proc_bias2_learned(q_thres(pt1),pdag,graph,cell_p,num_struc);
            results_EB_mmpc(pt1,t,d)=get_BY_FDR(p,p_thres(pt1))-get_emp_fdr_learned(pdag,graph,cell_p,p_thres(pt1),num_struc);
        end
       
    end
end
