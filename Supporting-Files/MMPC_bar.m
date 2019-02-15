function [CPC,sep,cell_p] = MMPC_bar(T,cond_indep, n, k, alpha, sep, cell_p, varargin)

CPC=[];
while ( length(CPC)>-1 )
    CPCo = CPC;
    [F,~,p,sep,cell_p]=MaxMinHeuristic(T,CPC,k,n,cond_indep,alpha,sep,cell_p,varargin{:});
    if p<alpha,
        CPC = [CPC,F];
    end
    if ( isequal(CPC,CPCo) )
        break
    end
end
for x = CPC
    ord = 0;
    k1=min(length(CPC),k);
    while ord<=k1,
         SS = nchoosek_vec(CPC,ord);
         for si=1:size(SS,1),
            S = SS(si,:);
            if (ismember(x,S))
                continue
            end
            [p,~]=feval(cond_indep, x, T, S, varargin{:});
            if p>alpha,
                if isempty(sep{x,T}),
                   sep{x,T} = [p,S];
                   sep{T,x} = [p,S];  
               elseif (length(sep{x,T})-1)>length(S) || ...
                       (((length(sep{x,T})-1)==length(S)) && sep{x,T}(1)<p),
                   sep{x,T} = [p,S];
                   sep{T,x} = [p,S];
                end
                CPC = setdiff(CPC,x);
                break;
            else
                cell_p{x,T}=[cell_p{x,T};p]; cell_p{T,x}=cell_p{x,T}; %*******
            end
        end
        ord = ord+1;
    end
end

end


function [F,assocF,p,sep,cell_p]=MaxMinHeuristic(T,CPC,k,n,cond_indep,alpha,sep,cell_p,varargin)

assocF=[];
assocF_p = [];
idx_x = setdiff(1:n,[T,CPC]);
for x=idx_x
    [r,p,sep,cell_p]=MinAssocEach(x,T,CPC,k,cond_indep,alpha,sep,cell_p,varargin{:}); 
    assocF=[assocF, r];
    assocF_p=[assocF_p, p];
end
F=find(assocF==max(assocF));
p=assocF_p(F);
assocF = assocF(F);
F=idx_x(F);

end

function [r,p,sep,cell_p]=MinAssocEach(x,T,CPC,k,cond_indep,alpha,sep,cell_p,varargin)
k1 = min(k,length(CPC));
ord=0;
rs=[];
ps=[];
while ord<=k1,
     if (isempty(CPC) || ord==0),
         [pt,rt]=feval(cond_indep, x, T, [], varargin{:});
         rs=[rs,rt];
         ps=[ps,pt];
         if pt>alpha,
             if isempty(sep{x,T}),
                 sep{x,T} = pt;
                 sep{T,x} = pt;  
             elseif (sep{x,T}(1)<pt)
                 sep{x,T} = pt;
                 sep{T,x} = pt;
             end
         else
             cell_p{x,T}=[cell_p{x,T};pt]; cell_p{T,x}=cell_p{x,T}; %*******
         end
     else,
         SS = nchoosek_vec(CPC,ord);
         for si=1:size(SS,1),
            S = SS(si,:);
            [pt,rt]=feval(cond_indep, x, T, S, varargin{:});
            rs=[rs,rt];
            ps=[ps,pt];
            if pt>alpha,
               if isempty(sep{x,T}),
                   sep{x,T} = [pt,S];
                   sep{T,x} = [pt,S];  
               elseif (length(sep{x,T})-1)>length(S) || ...
                       (((length(sep{x,T})-1)==length(S)) && sep{x,T}(1)<pt),
                   sep{x,T} = [pt,S];
                   sep{T,x} = [pt,S];
               end
            else
                cell_p{x,T}=[cell_p{x,T};pt]; cell_p{T,x}=cell_p{x,T}; %*******
            end
         end
     end
     ord=ord+1;
end
rs=abs(rs);
ri=find(rs==min(rs));
r=rs(ri);
p=ps(ri);
end

