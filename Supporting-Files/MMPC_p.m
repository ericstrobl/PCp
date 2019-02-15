function [G,sep,cell_p,num_struc] = MMPC_p(cond_indep, n, k, alpha, varargin)

CPCs={};
sep = cell(n,n);
cell_p = sep;
for T=1:n,
    [CPCs{T},sep,cell_p] = MMPC_bar(T,cond_indep, n, k, alpha, sep, cell_p, varargin{:});
end

CPCf={};
for T=1:n,
    CPC = CPCs{T};
    for v = CPC
        if ~ismember(T,CPCs{v})
            CPC = setdiff(CPC,v);
        end
    end
    CPCf{T} = CPC;
end

G=zeros(n,n);
G=sparse(setdiag(G,0));

for T=1:n,
    G(T,CPCf{T})=1;
    G(CPCf{T},T)=1;
end

% cell_p = cell(n,n); %*******
% for T=1:n,
%     CPC = CPCf{T};
%     for ne=CPC,
%         CPC2 = setdiff(CPC,[T,ne]);
%         for x=CPC2
%             ord=0;
%             CPC2t = setdiff(CPC,x);
%             k1=min(k,length(CPC2t));
%             while ord<=k1,
%                 SS = nchoosek_vec(CPC2t,ord);
%                 for si=1:size(SS,1),
%                    S = SS(si,:);
%                    [p,~]=feval(cond_indep, x, T, S, varargin{:});
%                    cell_p{x,T}=[cell_p{x,T};p]; cell_p{T,x}=cell_p{x,T}; %*******
%                 end
%                 ord = ord+1;
%             end
%         end
%                 
%         CPC3 = setdiff(CPCs{ne},[T,ne]);       
%         for x=CPC3
%             ord=0;
%             CPC3t = setdiff(CPC,x);
%             k1=min(k,length(CPC3t));
%             while ord<=k1,
%                 SS = nchoosek_vec(CPC3t,ord);
%                 if ismember(x,SS)
%                     continue
%                 end
%                 for si=1:size(SS,1),
%                    S = SS(si,:);
%                    [p,~]=feval(cond_indep, x, T, S, varargin{:});
%                    cell_p{x,T}=[cell_p{x,T};p]; cell_p{T,x}=cell_p{x,T}; %*******
%                 end
%                 ord = ord+1;
%             end
%         end
%         
%         
%     end
% end


idx=find(~cellfun(@isempty,cell_p)); %max of p-values
for t=1:length(idx),
    [i,j]= ind2sub(size(G),idx(t));
    if G(i,j)==1,
        cell_p{i,j}=max(cell_p{idx(t)});
    else
        cell_p{i,j}=[];
    end
end

idx=find(~cellfun(@isempty,sep)); %remove p-values from sep
for t=1:length(idx),
    [i,j]= ind2sub(size(G),idx(t));
    sep{i,j}(1)=[];
end


num_struc=assign_iden(G);