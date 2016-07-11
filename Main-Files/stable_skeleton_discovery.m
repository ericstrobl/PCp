function [G, sep, cell_p, num_struc] = stable_skeleton_discovery(cond_indep, d, k, alpha, varargin)
%
%
% Skeleton discovery procedure of PC-p
%
% This code adapted from the Bayes Net Toolbox
% (https://github.com/bayesnet/bnt) by Eric V. Strobl
%
% Inputs:
% 1) cond_indep = conditional independence test
% 2) n = number of variables
% 3) d = maximum conditioning set size.  If empty, uses default value of 3.
% 4) alpha = alpha threshold. If empty, uses default value of 0.20.
% 5) varargin = stuff necessary for the conditional independence test
% specified with cond_indep (e.g., the data, hyperparameters)
%
% Outputs:
% 1) G = skeleton. G(i,j)=1 and G(j,i)=1 iff adjacency between i and j
% 2) sep = separating sets
% 3) cell_p = p-values of each edge
% 4) num_struc = cell that keeps track of number of hypothesis tests
%
% Example call: 
% [Gt, sep, cell_pt, num_struct] = get_skeleton_lin_stable(@rho_test_PC, size(data,2), 3, [], data);
%
%

if isempty(k), k=3; end
if isempty(alpha), alpha=0.20; end
  
sep = cell(d,d);
cell_p = cell(d,d); % cell of p-values

ord = 0;
done = 0;

G=ones(d,d);
G=sparse(setdiag(G,0));

while ~done
  done = 1;
  [X,Y] = find(G); 
  
  for k1=1:size(G,1), % changing conditioning sets according to PC-stable's procedure
      nbrs2{k1} = neighbors2(G, k1);  %*******
  end %*******
  
  for i=1:length(X)
    x = X(i); y = Y(i);
    nbrs = mysetdiff(nbrs2{y}, x);      
    if length(nbrs) >= ord && G(x,y) ~= 0
      done = 0;
      SS = nchoosek(nbrs,ord);
      
        for si=1:size(SS,1),
            S = SS(si,:);
            [p,~]=feval(cond_indep, x, y, S, varargin{:});
            
            if p<=alpha, % record p-value if below threshold
                cell_p{x,y}=[cell_p{x,y};p]; cell_p{y,x}=cell_p{x,y}; %*******
            end %*******
            
            if p>alpha, % if p-value above threshold, 
              G(x,y) = 0;
              G(y,x) = 0;
              sep{x,y} = myunion(sep{x,y}, S);
              sep{y,x} = myunion(sep{y,x}, S);
              break;
            end
        end
        
    end 
  end
  ord = ord + 1;
  if ord > k,
      break;
  end
end

idx=find(~cellfun(@isempty,cell_p)); %takes max of surviving p-values (i.e. equation 16 in the paper)
for t=1:length(idx),
    [i,j]= ind2sub(size(G),idx(t));
    if G(i,j)==1,
        cell_p{i,j}=max(cell_p{idx(t)});
    else
        cell_p{i,j}=[];
    end
end

num_struc=assign_iden(G);
