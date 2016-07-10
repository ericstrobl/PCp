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
