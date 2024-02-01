clc;clear
%Optimization input
n=10; m=15; N=n+m; nm=n*m; B0=ones(n,m);
% load('A1Geo'); load('A2Geo'); n=max(size(A1)); m=max(size(A2)); N=n+m; nm=n*m; B0=ones(n,m);
% G1 = graph(A1); L0_1=laplacian(G1); G2 = graph(A2); L0_2=laplacian(G2);
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)],  [e1(2)/n e2(2)/m],
%% BA Layers
% n1=15; seed1 = rand(n1,n1) < 1; seed1 = triu(seed1,1); seed1 = seed1 + seed1';
% n2=6; seed2 = rand(n2,n2) < 1; seed2 = triu(seed2,1); seed2 = seed2 + seed2';
% mlink1=2; mlink2=1; A1 = SFNG(n, mlink1, seed1); A2 = SFNG(m, mlink2, seed2);
% G1 = graph(A1); L0_1=laplacian(G1);
% G2 = graph(A2); L0_2=laplacian(G2);
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)],  [e1(2)/n e2(2)/m],
%% ER layers
% p1=.2; p2=.15; G1=ErdosRenyi(n,p1); G2=ErdosRenyi(m,p2); 
% L0_1=laplacian(G1); L0_2=laplacian(G2);
% A1=adjacency(G1); A2=adjacency(G2); 
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)], [e1(2)/n e2(2)/m],
%% Geo layers
r1=.8*sqrt(2*log(n)/n);    
r2=.7*sqrt(2*log(m)/m); 
A1=NetGen_Geo(n,r1);
A2 =NetGen_Geo(m,r2);
G1 = graph(A1); L0_1=laplacian(G1);
G2 = graph(A2); L0_2=laplacian(G2);
e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
[e1(2) e2(2)],
[e1(2)/n e2(2)/m],
%% WS layers
% K1=2; beta1=.5; K2=1;  beta2=1; 
% G1 = WattsStrogatz(n,K1,beta1); L0_1=laplacian(G1);
% G2 = WattsStrogatz(m,K2,beta2); L0_2=laplacian(G2);
% A1=adjacency(G1); A2=adjacency(G2); 
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)]
% [e1(2)/n e2(2)/m]
%% Combined layers
% r1=.8*sqrt(2*log(n)/n); A1=NetGen_Geo(n,r1);   G1 = graph(A1); L0_1=laplacian(G1);
% p2=.2; G2=ErdosRenyi(m,p2);   L0_2=laplacian(G2);   A2=adjacency(G2);
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)]
% [e1(2)/n e2(2)/m]
%% Computing the laplacian L0
L0=[L0_1 zeros(n,m);zeros(m,n) L0_2];
% computing ee'
ee=ones(N,N); un=ones(n,1); um=ones(m,1);
%% cvx
l2=10; %c=100;

    cvx_begin quiet sdp
        cvx_solver mosek
        % cvx_solver sedumi
        % cvx_solver SDPT3

        variables w(n,m) nu;
        expression L(N,N); 
        L=[diag(w*um) -w;-w' diag(w'*un)];
        minimize( norm(vec(w),1) );
        subject to
            L + L0 +nu*ones(N,N)-l2*eye(N)>=0;
            vec(w) >=0;
%             sum(vec(w))==c;
     cvx_end
  w   