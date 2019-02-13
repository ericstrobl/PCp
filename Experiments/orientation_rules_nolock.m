function [pdag, cell_p, num_struc] = orientation_rules_nolock(G, pdag, cell_p, num_struc)
% Pearl (2000) p51

old_pdag = zeros(size(G,1));
iter = 0;

max_num=max(max(num_struc));
while ~isequal(pdag, old_pdag)
  iter = iter + 1;
  old_pdag = pdag;
  
  pdag2=pdag;%******
  cell_p_d=cell(size(pdag,1),size(pdag,1));
  
  % rule 1
  [A,B] = find(pdag==-1); % a -> b
  for i=1:length(A)
    a = A(i); b = B(i);
    C = find(pdag(b,:)==1 & G(a,:)==0); % all nodes adj to b but not a
    if ~isempty(C)
      pdag2(b,C) = -1; 
      %fprintf('rule 1: a=%d->b=%d and b=%d-c=%d implies %d->%d\n', a, b, b, C, b, C);
       
      for t=1:length(C),
         cell_p_d{b,C(t)}=sum([cell_p_d{b,C(t)},cell_p{a,b}]);
%             cell_p2{b,C(t)}=cell_p{a,b};
      end

    end
  end
  
  % rule 2
  [A,B] = find(pdag==1); % unoriented a-b edge
  for i=1:length(A)
    a = A(i); b = B(i);
    if any( (pdag(a,:)==-1) & (pdag(:,b)==-1)' );
      pdag2(a,b) = -1; 
      %fprintf('rule 2: %d -> %d\n', a, b);

      idx= find( (pdag(a,:)==-1) & (pdag(:,b)==-1)' );

      p_t=[];
      for t=1:length(idx),
          p_t=[p_t; cell_p{a,idx(t)}, cell_p{idx(t),b}];
      end
      p_t=sum(max(p_t,[],2));

      cell_p_d{a,b}=max([0,cell_p_d{a,b}])+p_t; %******
    end
  end
  
  % rule 3
  [A,B] = find(pdag==1); % a-b
  for i=1:length(A)
    a = A(i); b = B(i);
    C = find( (pdag(a,:)==1) & (pdag(:,b)==-1)' );
    % C contains nodes c s.t. a-c->b
    G2 = setdiag(G(C, C), 1);
    if any(G2(:)==0) % there are 2 different non adjacent elements of C
      pdag2(a,b) = -1; 

      p_t=[];
      for t=1:length(C),
         p_t=[p_t; cell_p{a,C(t)}, cell_p{C(t),b}];
      end
      p_t=max(p_t,[],2);
      r1=nchoosek(1:length(C),2);
      
      p_t_s=0;
      for t=1:size(r1,1),
          p_t_t=max([p_t(r1(t,1)),p_t(r1(t,2))]);
          p_t_s=sum([p_t_s,p_t_t]);
      end

      cell_p_d{a,b}=max([0,cell_p_d{a,b}])+p_t_s;  %******
      
      %fprintf('rule 3: %d -> %d\n', a, b);
    end
  end
  


%     %% clamp bidirected edges
%     idx=find(pdag2==-1);
%     for t=1:length(idx),
%         [i,j]= ind2sub(size(G),idx(t));
%         if pdag2(i,j)==-1 && pdag2(j,i)==-1,  
%             pdag2(i,j)=2;
%             pdag2(j,i)=2;
%             G(i,j)=2;
%             G(j,i)=2;
%         end
%     end
  
  idxLa=find(~cellfun(@isempty,cell_p_d));  %******
  for t=1:length(idxLa),  %******
      [r,c]=ind2sub(size(cell_p),idxLa(t));
          
      cell_p{r,c}=max([cell_p_d{r,c}, cell_p{r,c}]);  %****** transfer all temporary p-values to permanent p-values
      max_num=max_num+1;
      num_struc(r,c)=max_num;
      pdag(r,c)=pdag2(r,c);
%       if ~ismember(sub2ind(size(cell_p),c,r),idxLa), % erase other direction
%           cell_p{c,r}=[];
%           num_struc(c,r)=0;
%           pdag(c,r)=0;
%       end
  end  %******

end