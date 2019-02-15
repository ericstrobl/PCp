function [pdag, cell_p, num_struc] = get_v_structures2_original(G, sep, cell_p, num_struc, k, cond_indep, varargin)

pdag = G;
[X, Y] = find(G);
cell_p2=cell(size(G,1),size(G,1));

max_num=max(max(num_struc));
for i=1:length(X)
  x = X(i);
  y = Y(i);
  Z = find(G(y,:));
  Z = mysetdiff(Z, x);
  for z=Z(:)'
    if G(x,z)==0 && ~ismember(y, sep{x,z}) && ~ismember(y, sep{z,x}),
        
      if (pdag(x,y) ~= -1 || pdag(z,y) ~= -1),

          pdag(x,y) = -1; 
          pdag(z,y) = -1; 
          
          pdag(y,x) = 0;
          num_struc(y,x)=0;
          pdag(y,z) = 0;
          num_struc(y,z)=0;
          
          % p-value
          p2=[];
          ii1 = find(G(x,:)>0); % neighbors of x
          ii1 = mysetdiff(ii1,x);
          SS1a=[];
          for kk=1:k,
              if length(ii1)>1,
                 SS1=nchoosek(ii1,kk);
                 SS1=[SS1 zeros(size(SS1,1),k-kk)];
                 SS1a=[SS1a; SS1];
              else
                 SS1=ii1;
                 SS1=[SS1 zeros(size(SS1,1),k-1)];
                 SS1a=[SS1a; SS1];
                 break;
              end
          end
  
          ii2 = find(G(z,:)>0); % neighbors of z
          ii2 = mysetdiff(ii2,z);
          SS2a=[];
          for kk=1:k,
              if length(ii2)>1,
                 SS2=nchoosek(ii2,kk);
                 SS2=[SS2 zeros(size(SS2,1),k-kk)];
                 SS2a=[SS2a; SS2];
              else
                 SS2=ii2;
                 SS2=[SS2 zeros(size(SS2,1),k-1)];
                 SS2a=[SS2a; SS2];
                 break;
              end
          end
          
          SS=[SS1a;SS2a]; % combine all subsets
          for t=1:size(SS,1),
              cond=SS(t);
              cond=cond(cond>0);
              [p2t,~]=feval(cond_indep, x, z, [z,cond], varargin{:});
%               [~,p2t]=partialcorr(data(:,x),data(:,z),data(:,[z,cond]));
              p2=[p2; p2t];
          end
          p2=max(p2);
         
          
          cell_p2{z,y}=max([cell_p{x,y},p2]); %*******
          cell_p2{x,y}=max([cell_p{z,y},p2]); %*******
          
          cell_p2{y,z}={}; %*******
          cell_p2{y,x}={}; %*******
%           cell_p2{y,x}=[]; cell_p2{y,z}=[]; %*******

          if length(cell_p2{z,y})==1 && length(cell_p2{x,y})==1,
            num_struc(x,y)=max_num+1;
            num_struc(z,y)=max_num+1;
          else
            num_struc(x,y)=max_num+1;
            num_struc(z,y)=max_num+2;
          end
          max_num=max(max(num_struc));
          
      
      end
      
    end
  end
end

%erase p-values of removed edges
idx=find(pdag==0);
for t=1:length(idx),
    cell_p{idx(t)}=[];
end

% modify cell_p
idx=find(~cellfun(@isempty,cell_p2));
for t=1:length(idx),
    cell_p{idx(t)}=max([cell_p{idx(t)},cell_p2{idx(t)}]);
end
