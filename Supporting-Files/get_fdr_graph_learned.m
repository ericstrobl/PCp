function fdr=get_fdr_graph_learned(est_graph,true_graph,num_struc)

R=max([length(unique(num_struc))-1,1]);

idx_d=find(est_graph==-1);
idx_d_r=find(true_graph==1);
[~,idx_dt]=setdiff(idx_d,idx_d_r);
V1=length(unique(num_struc(idx_d(idx_dt))));

idx_u=find(est_graph==1);
idx_u_r=find((true_graph+true_graph')>=1);
[~,idx_ut]=setdiff(idx_u,idx_u_r);
V2=length(unique(num_struc(idx_u(idx_ut))));


fdr=(V1+V2)/R;
end