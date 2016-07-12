function [data,coeff]=create_dataset_dag(G,samps,coeff)

r=size(G,1);
% coeff=G.*(rand(r,r)+0.5);
if isempty(coeff),
    coeff=G.*randn(r,r);
end


err_std=rand(1,r);  
err_std=repmat(err_std,samps,1);
err=err_std.*randn(samps,r);

data=randn(samps,r);

done=find(sum(G,1)==0);
stop=0;
while stop==0,
    for s=done,
        ch=find(G(s,:)==1);
        for c=ch,
            pa=find(G(:,c)==1);
            pa=sort(pa,'ascend');
            
            h=intersect(pa,done);
            
            data(:,c)=(data(:,h)*coeff(h,c))+err(:,c);
            done=unique([done, c]);
             
        end
    end
    
    if length(done)==r,
        stop=1;
    end
    
end
    