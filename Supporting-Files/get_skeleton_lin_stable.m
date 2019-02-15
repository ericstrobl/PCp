function [G, sep, cell_p, num_struc] = get_skeleton_lin_stable(cond_indep, n, k, alpha, varargin)


  
sep = cell(n,n);
cell_p = cell(n,n); %*******

ord = 0;
done = 0;

G=ones(n,n);
G=sparse(setdiag(G,0));

while ~done
  done = 1;
  [X,Y] = find(G); 
  
  for k1=1:size(G,1), %*******
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
            
            if p<=alpha, %*******
                cell_p{x,y}=[cell_p{x,y};p]; cell_p{y,x}=cell_p{x,y}; %*******
            end %*******
            
            if p>alpha,
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

idx=find(~cellfun(@isempty,cell_p)); %max of p-values
for t=1:length(idx),
    [i,j]= ind2sub(size(G),idx(t));
    if G(i,j)==1,
        cell_p{i,j}=max(cell_p{idx(t)});
    else
        cell_p{i,j}=[];
    end
end

num_struc=assign_iden(G);