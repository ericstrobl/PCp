function SS = nchoosek_vec(vec,n)

if n==0 || isempty(vec)
    SS = [];
elseif length(vec)==1 && n==1
    SS = vec;
elseif length(vec)>1
    SS = nchoosek(vec,n);
end