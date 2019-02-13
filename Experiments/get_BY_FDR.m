function FDR=get_BY_FDR(p,alpha)

nt=length(p);
den=sum(p<=alpha);
if den==0, den=1; end
FDR=(nt*alpha*sum(1./(1:nt)))/den;