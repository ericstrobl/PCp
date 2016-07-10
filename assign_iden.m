function matr=assign_iden(G)

matr=tril(G);
idx=find(matr==1);
matr(matr==1)=1:length(idx);

matr=matr+tril(matr)';