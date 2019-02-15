function proc_bias=get_proc_bias2_learned(q,pdag,graph,cell_p,num_struc)

p=correct_p_values(cell_p,num_struc);

t=binary_search(1E-100,1,p,q);
% fdr_t=get_BY_FDR(p,t);
% fdr_e=get_emp_fdr(pdag,graph,cell_p,fdr_t,num_struc);
fdr_e=get_emp_fdr_learned(pdag,graph,cell_p,t,num_struc);
proc_bias=q-fdr_e;
end

function a2=binary_search(start,last,p,q)
cent=mean([start,last]);
for t=1:50,
    a1=mean([start,cent]);
    ans1=get_BY_FDR(p,a1);
    a2=mean([cent,last]);
    ans2=get_BY_FDR(p,a2);
    
    if ans1>q,
        last=a1;
    elseif ans2>q,
        last=a2;
    else
        start=a2;
    end
    cent=mean([start,last]);
        
end

end