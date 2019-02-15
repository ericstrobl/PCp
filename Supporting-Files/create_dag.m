function B= create_dag(en,p)

  N = (p*p - p)/2;
  
  samplesB = binornd(1, en/(p-1),[N,1] );
  
  A=ones(p,p);
  B=zeros(p,p);
  B(triu(A, 1)==1) = samplesB;
  B=B.*((0.9.*rand(p)+0.1).*vec2mat(randsample([-1,1],p^2,true),p));
