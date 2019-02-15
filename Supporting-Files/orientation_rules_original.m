function [pdag, cell_p, num_struc] = orientation_rules_original(G, pdag, cell_p, num_struc)
% Pearl (2000) p51

old_pdag = zeros(size(G,1));
iter = 0;

max_num=max(max(num_struc));
while ~isequal(pdag, old_pdag)
  iter = iter + 1;
  old_pdag = pdag;
  
  % rule 1
  [A,B] = find(pdag==-1); % a -> b
  for i=1:length(A)
    a = A(i); b = B(i);
    C = find(pdag(b,:)==1 & G(a,:)==0); % all nodes adj to b but not a
    if ~isempty(C)
      pdag(b,C) = -1; 
      pdag(C,b) = 0;
      %fprintf('rule 1: a=%d->b=%d and b=%d-c=%d implies %d->%d\n', a, b, b, C, b, C);
       
      for t=1:length(C),
         cell_p{b,C(t)}=max([cell_p{b,C(t)},cell_p{a,b}]);
         max_num=max_num+1;
         num_struc(b,C(t))=max_num;
         num_struc(C(t),b)=0;
         cell_p{C(t),b}=[];
      end

    end
  end
  
%   eidx=find(~cellfun(@isempty,cell_p_d));
%   for t=1:length(eidx),
%       cell_p{eidx(t)}=max([cell_p{eidx(t)},cell_p_d{eidx(t)}]); %******
%   end
  
  % rule 2
  [A,B] = find(pdag==1); % unoriented a-b edge
  for i=1:length(A)
    a = A(i); b = B(i);
    if any( (pdag(a,:)==-1) & (pdag(:,b)==-1)' );
      pdag(a,b) = -1; 
      pdag(b,a) = 0;
      %fprintf('rule 2: %d -> %d\n', a, b);

      idx= find( (pdag(a,:)==-1) & (pdag(:,b)==-1)' );

      p_t=[];
      for t=1:length(idx),
          p_t=[p_t; cell_p{a,idx(t)}, cell_p{idx(t),b}];
      end
      p_t=sum(max(p_t,[],2));

      cell_p{a,b}=max([cell_p{a,b},p_t]); %******
      max_num=max_num+1;
      num_struc(a,b)=max_num;
      num_struc(b,a)=0;
      cell_p{b,a}=[];
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
      pdag(a,b) = -1; 
      pdag(b,a) = 0;

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

      cell_p{a,b}=max([cell_p{a,b},p_t_s]);  %******
      max_num=max_num+1;
      num_struc(a,b)=max_num;
      num_struc(b,a)=0;
      cell_p{b,a}=[];
      
      %fprintf('rule 3: %d -> %d\n', a, b);
    end
  end


end