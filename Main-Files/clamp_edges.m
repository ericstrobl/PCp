function [pdag,G,cell_p2,num_struc]=clamp_edges(pdag,G,cell_p2,cell_p,num_struc)

%clamp all bidirected edges
%also, replace p-values with undirected edges
idx=find(pdag==-1);
for t=1:length(idx),
    [i,j]= ind2sub(size(G),idx(t));
    if pdag(i,j)==-1 && pdag(j,i)==-1,  
        
        pdag(i,j)=2;
        pdag(j,i)=2;
        
        G(i,j)=2;
        G(j,i)=2;
        
        cell_p2{i,j}=cell_p{i,j};
        cell_p2{j,i}=cell_p{j,i};
        
    end
end

%clamp all directed edges connected to the bidirected edges
%also, replace p-values with undirected edges
idx=find(pdag==2);
for t=1:length(idx),
    [i,j]= ind2sub(size(G),idx(t));
    
    idx1=find(pdag(:,j)==-1)'; %connected to j
%     idx2=find(pdag(i,:)==-1); %connected from i
    idx2=find(pdag(:,i)==-1)'; %connected to i
    
    G(idx1,j)=2; G(j,idx1)=2;
    G(i,idx2)=2; G(idx2,i)=2;
    
    if ~isempty(idx1),
        for t1=idx1,
            cell_p2{t1,j}=cell_p{t1,j};
            cell_p2{j,t1}=cell_p2{t1,j};
        end
    end
    
    if ~isempty(idx2),
        for t2=idx2,
            cell_p2{i,t2}=cell_p{i,t2};
            cell_p2{t2,i}=cell_p2{i,t2};
        end
    end
    
    pdag(idx1,j)=2; pdag(j,idx1)=2;
    pdag(i,idx2)=2; pdag(idx2,i)=2;
end

%adjust number of structures
max_num=max(max(num_struc));
idx=find(pdag==2);
for t=1:length(idx),
    [r,c]=ind2sub(size(G),idx(t));
    max_num=max_num+1;
    num_struc(r,c)=max_num;
    num_struc(c,r)=max_num;
end