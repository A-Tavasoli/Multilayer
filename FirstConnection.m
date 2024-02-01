clc;clear
%Optimization input
n=50; m=50; N=n+m; nm=n*m;  c=10;
B0=ones(n,m);
% load('A1'); load('A2'); n=max(size(A1)); m=max(size(A2)); N=n+m; nm=n*m; B0=ones(n,m);
% G1 = graph(A1); L0_1=laplacian(G1); G2 = graph(A2); L0_2=laplacian(G2);
%% BA Layers
% n1=15; seed1 = rand(n1,n1) < 1; seed1 = triu(seed1,1); seed1 = seed1 + seed1';
% n2=10; seed2 = rand(n2,n2) < 1; seed2 = triu(seed2,1); seed2 = seed2 + seed2';
% mlink1=5; mlink2=2; A1 = SFNG(n, mlink1, seed1); A2 = SFNG(m, mlink2, seed2);
% G1 = graph(A1); L0_1=laplacian(G1);
% G2 = graph(A2); L0_2=laplacian(G2);
% e1 = eig(L0_1); e2 = eig(L0_2);
% [e1(2) e2(2)],  
%% ER layers
% p1=.3; p2=.7; G1=ErdosRenyi(n,p1); G2=ErdosRenyi(m,p2); 
% L0_1=laplacian(G1); L0_2=laplacian(G2);
% A1=adjacency(G1); A2=adjacency(G2); 
%% Geo layers
% r1=.8*sqrt(2*log(n)/n);    
% r2=.7*sqrt(2*log(m)/m); 
% A1=NetGen_Geo(n,r1);
% A2 =NetGen_Geo(m,r2);
% G1 = graph(A1); L0_1=laplacian(G1);
% G2 = graph(A2); L0_2=laplacian(G2);
% e1 = eig(L0_1); e2 = eig(L0_2); [e1(2) e2(2)],
%% WS layers
K1=2; K2=2; beta1=.1; beta2=1; 
G1 = WattsStrogatz(n,K1,beta1); L0_1=laplacian(G1);
G2 = WattsStrogatz(m,K2,beta2); L0_2=laplacian(G2);
A1=adjacency(G1); A2=adjacency(G2); 
e1 = eig(L0_1); e2 = eig(L0_2);
[e1(2) e2(2)]
%% Combined layers
% r1=1*sqrt(2*log(n)/n); A1=NetGen_Geo(n,r1);   
% p2=.4; G2=ErdosRenyi(m,p2);   L0_2=laplacian(G2);   A2=adjacency(G2);
% e1 = eig(L0_1); e2 = eig(L0_2);
% [e1(2) e2(2)]
% [e1(2)/n e2(2)/m]
%% Computing the laplacian L0
L0=[L0_1 zeros(n,m);zeros(m,n) L0_2];
% computing ee'
ee=ones(N,N); un=ones(n,1); um=ones(m,1);
%% cvx
% [u,v]=eig(full(L0_1));  u2_1=u(:,2);   lam2_1=v(2,2),
% [u,v]=eig(full(L0_2));  u2_2=u(:,2);   lam2_2=v(2,2),
% [k,l]=min(abs(u2_1)); i=l;
% [k,l]=min(abs(u2_2)); j=l+n;
k=0; lam2=zeros(n,m);
for i=1:n
    for j=1:m
        k=k+1;
        a=zeros(N,1); a(i)=1; a(j+n)=-1;
        L=L0+c*a*a';  
        u=eig(L); 
        lam2(i,j)=u(2); 
    end
end
[k,i]=max(max(lam2));  
[k,j]=max(max(lam2')); [j i],
[cent1,CentNode1]=max(centrality(G1,'eigenvector')); 
[cent2,CentNode2]=max(centrality(G2,'eigenvector')); [CentNode1 CentNode2],
figure; plot(G1)
figure; plot(G2,'r')