function fdr = get_emp_fdr_learned(est_graph,true_graph,cell_p,p_thres,num_struc)

est_graph=get_est_graph(est_graph, cell_p, p_thres);

fdr=get_fdr_graph_learned(est_graph,true_graph,num_struc);

end