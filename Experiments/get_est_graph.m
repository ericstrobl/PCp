function est_graph=get_est_graph(est_graph, cell_p, p_thres)

idx=find(~cellfun(@isempty,cell_p));
ps=zeros(1,length(idx));
for t=1:length(idx),
    ps(t)=cell_p{idx(t)};
end

idx_p=find(ps>p_thres);
est_graph(idx(idx_p))=0;


end