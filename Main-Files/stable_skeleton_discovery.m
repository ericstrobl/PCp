function [G, sep, p_val, IDs] = stable_skeleton_discovery(cond_indep, alpha, k, d, n, varargin)
%
% PC-stable's skeleton discovery procedure
%
% Inputs:
% 1) cond_indep = conditional independence test function handle. The test
% should be of the form CItest(x,y,Z,varargin).
% 2) alpha = alpha threshold. If empty, uses heuristic default value that 
% starts at 0.20 and then decreases at rate root n, where n denotes sample 
% size; this is again just a heuristic but decreasing the alpha value
% ensures asymptotic consistency. If you want to be extra conservative, set 
% alpha=0.20 no matter what sample size. The theory says that the alpha
% threshold should be high *for the problem at hand*, so the Type II error 
% rate is low.
% 3) k = maximum conditioning set size. If empty, algorithm dynamically 
% adjusts the conditioning set size for each node according to whatever 
% size is necessary. Set this if you know the maximum in-degree.
% 4) d = number of variables
% 5) n = sample size
% 6) varargin = stuff necessary for the conditional independence test
% specified with cond_indep (e.g., the data, hyperparameters)
%
% Outputs:
% 1) G = inferred skeleton. G(i,j)=G(j,i)=1 iff adjacency between i and j
% 2) sep = separating sets
% 3) p_val = p-values of each edge
% 4) IDs = cell that keeps track of number of hypothesis tests
%
% Example call: 
% [Gt, sep, cell_pt, num_struct] = get_skeleton_lin_stable(@rho_test_PC, size(data,2), 3, [], data);
%
% Code written by Eric V. Strobl using the Bayes Net Toolbox
% (https://github.com/bayesnet/bnt) as a base.
%

if isempty(alpha),
%     alpha=0.20;
    if n<=500, alpha=0.20;
    else alpha=0.20/sqrt(n/500);
    end
end
  
sep = cell(d,d);
cell_p = cell(d,d); % cell of p-values

ord = 0;
done = 0;

G=ones(d,d);
G=sparse(setdiag(G,0));

while ~done
  done = 1;
  [X,Y] = find(G); 
  
%   sort(r_idx(X))
  
  for k1=1:size(G,1), % changing conditioning sets according to PC-stable's procedure
      nbrs2{k1} = neighbors2(G, k1);  %*******
  end %*******

%  sort(r_idx(nbrs2{r_idx==10}))
  
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
  
  
  if isempty(k),
     if ord > max(sum(G,2)),
        break;
    end
  else
    if ord > k,
        break;
    end
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


p_val=cell_p;
IDs=num_struc;