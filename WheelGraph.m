function G=WheelGraph(n)

K=1; beta=0; G = WattsStrogatz(n-1,K,beta); A=adjacency(G);  

A(:,n)=ones(n-1,1); A(n,:)=ones(1,n); A(n,n)=0;
G=graph(A);
