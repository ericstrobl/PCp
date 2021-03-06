%% this file runs all of the algorithms and saves their results

%% low dimensional experiments
disp('running low dimensional experiments...');
n= 1000;
k=10;
alphas=[0.01 0.2/sqrt(n/250)];

for g=1:100
    disp(g)
    B= create_dag(2,20);
    data=create_dataset_dag(B>0,n,B);
    data_t=data;
    Ct = corr(data_t);
    for a=1:length(alphas)
        tic;
        [Gt, sep, cell_pt, num_struct] = get_skeleton_lin_stable(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        time_p=toc;

        tic;
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_pcp_graphLD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p_norobust(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_pcpnorobust_graphLD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [pdag,cell_p2,num_struc,G]=PC_p_noambig(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_noambig_graphLD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        %% no stable
        tic;
        [Gt, sep, cell_pt, num_struct] = get_skeleton_lin(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        time_p=toc;

        tic;
        [pdag,cell_p,num_struc,G]=PC_original(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_original_graphLD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_nostable_graphLD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p_noambig(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_nostablenoambig_graphLD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [Gt,sep,cell_pt,num_struct]=MMPC_p(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=toc;
        save(['MMPC_graphLD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

    end
    
end


%% high dimensional experiments
disp('running high dimensional experiments...');

k=3;
for g=1:100
    disp(g)
    B= create_dag(2,100);
    data=create_dataset_dag(B>0,n,B);
    data_t=data;
    Ct = corr(data_t);
    for a=1:length(alphas)
        tic;
        [Gt, sep, cell_pt, num_struct] = get_skeleton_lin_stable(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        time_p=toc;

        tic;
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_pcp_graphHD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p_norobust(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_pcpnorobust_graphHD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [pdag,cell_p2,num_struc,G]=PC_p_noambig(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_noambig_graphHD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        %% no stable
        tic;
        [Gt, sep, cell_pt, num_struct] = get_skeleton_lin(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        time_p=toc;

        tic;
        [pdag,cell_p,num_struc,G]=PC_original(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_original_graphHD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_nostable_graphHD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p_noambig(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_nostablenoambig_graphHD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')

        tic;
        [Gt,sep,cell_pt,num_struct]=MMPC_p(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=toc;
        save(['MMPC_graphHD1_synth' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time','B')
    end
    
end


%% CYTO experiments
disp('running CYTO experiments...');

% make sure that these files are in the working directory
load('CYTO_data_cd2cd28icam2.mat')
load('CYTO_Dlingam.mat')
%

graph=Bsig>0;

n=size(data,1);
k=10;
alphas=[0.01 0.2/sqrt(n/250)];


[r,~]=size(data);

for g=1:100
    disp(g)
    data_t=create_dataset_dag(Bsig>0,size(data,1),Bsig);
    Ct = corr(data_t);
    for a=1:length(alphas)
        tic;
        [Gt, sep, cell_pt, num_struct] = get_skeleton_lin_stable(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        time_p=toc;

        tic;
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_pcp_graphCYTO1_real' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p_norobust(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_pcpnorobust_graphCYTO1_real' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time')

        tic;
        [pdag,cell_p2,num_struc,G]=PC_p_noambig(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_noambig_graphCYTO1_real' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time')

        %% no stable
        tic;
        [Gt, sep, cell_pt, num_struct] = get_skeleton_lin(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        time_p=toc;

        tic;
        [pdag,cell_p,num_struc,G]=PC_original(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_original_graphCYTO1_real' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_nostable_graphCYTO1_real' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time')

        tic;
        [pdag,cell_p,num_struc,G]=PC_p_noambig(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=time_p+toc;
        save(['PC_nostablenoambig_graphCYTO1_real' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time')

        tic;
        [Gt,sep,cell_pt,num_struct]=MMPC_p(@gaussCItest, size(data,2), k, alphas(a), Ct, n);
        [pdag,cell_p,num_struc,G]=PC_p(Gt, sep, cell_pt, num_struct, k, @gaussCItest, Ct, n);
        time=toc;
        save(['MMPC_graphCYTO1_real' num2str(g) '_' num2str(a) '.mat'],'pdag','cell_p','num_struc','Gt','cell_pt','num_struct','time')
                
    end
    
end
