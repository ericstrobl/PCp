function [p,r]=gaussCItest(x,y,z,C,n)

if isempty(z)
    r = C(x,y);
elseif length(z)==1
    r = (C(x,y) - C(x,z) * C(y,z))/sqrt((1 - C(y,z)^2)*(1 - C(x,z)^2));
else
    PM = pinv(C([x, y, z], [x, y, z]));
    r = PM(1, 2)/sqrt(PM(1, 1) * PM(2, 2));
end

r = sqrt(n - length(z) - 3) * 0.5 * log(1+(2*r/(1-r)));
if isnan(r)
    r=0;
end
p = 2*(1-normcdf(abs(r)));