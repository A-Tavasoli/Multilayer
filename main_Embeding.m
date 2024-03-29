clc;clear
%Optimization input
c=50; 
%n=30; 
% m=15;
load('A1GeoCase1'); 
load('A2GeoCase1'); 
n=max(size(A1)); 
m=max(size(A2)); 
N=n+m; nm=n*m; 
B0=ones(n,m); 
I=eye(N,N);
G1 = graph(A1); 
L0_1=laplacian(G1); 
G2 = graph(A2); 
L0_2=laplacian(G2);
e1 = eig(L0_1); 
e2 = eig(L0_2); 
[eigen0_1,e] = eig(full(L0_1)); 
[eigen0_2,e] = eig(full(L0_2));
%% BA Layers
% n1=9; 
% seed1 = rand(n1,n1) < 1; 
% seed1 = triu(seed1,1); 
% seed1 = seed1 + seed1';
% n2=2; 
% seed2 = rand(n2,n2) < 1; 
% seed2 = triu(seed2,1); 
% seed2 = seed2 + seed2';
% mlink1=3; 
% mlink2=1; 
% A1 = SFNG(n, mlink1, seed1); 
% A2 = SFNG(m, mlink2, seed2);
% G1 = graph(A1); 
% L0_1=laplacian(G1);
% G2 = graph(A2); 
% L0_2=laplacian(G2);
% e1 = eig(L0_1); 
% e2 = eig(L0_2);
%% ER layers
% p1=.2; p2=.3; 
% G1=ErdosRenyi(n,p1); 
% G2=ErdosRenyi(m,p2); 
% L0_1=laplacian(G1); 
% L0_2=laplacian(G2);
% A1=adjacency(G1); 
% A2=adjacency(G2); 
% e1 = eig(L0_1); 
% e2 = eig(L0_2); 
%% Geo layers
% r1=.8*sqrt(2*log(n)/n);    
% r2=.6*sqrt(2*log(m)/m); 
% A1=NetGen_Geo(n,r1);
% A2 =NetGen_Geo(m,r2);
% G1 = graph(A1); 
% L0_1=laplacian(G1);
% G2 = graph(A2); 
% L0_2=laplacian(G2);
% e1 = eig(L0_1); 
% e2 = eig(L0_2); 
% [eigen0_1,e] = eig(full(L0_1)); 
% [eigen0_2,e] = eig(full(L0_2));
%% WS layers
% K1=2; K2=1; 
% beta1=.1; 
% beta2=1; 
% G1 = WattsStrogatz(n,K1,beta1); 
% L0_1=laplacian(G1);
% G2 = WattsStrogatz(m,K2,beta2); 
% L0_2=laplacian(G2);
% A1=adjacency(G1); 
% A2=adjacency(G2); 
% e1 = eig(L0_1); 
% e2 = eig(L0_2);
%% Combined layers
% r1=1*sqrt(2*log(n)/n); 
% A1=NetGen_Geo(n,r1);   
% G1 = graph(A1); 
% L0_1=laplacian(G1);
% p2=.4; G2=ErdosRenyi(m,p2);   
% L0_2=laplacian(G2);   
% A2=adjacency(G2);
% p1=.15; G1=ErdosRenyi(n,p1);   
% L0_1=laplacian(G1);   
% A1=adjacency(G1);
% K2=1; beta2=0; 
% G2 = WattsStrogatz(m,K2,beta2); 
% L0_2=laplacian(G2); 
% A2=adjacency(G2);  
% e1 = eig(L0_1); 
% e2 = eig(L0_2);
%% Computing the laplacian L0
L0=[L0_1 zeros(n,m);zeros(m,n) L0_2];
% computing ee'
ee=ones(N,N); 
un=ones(n,1); 
um=ones(m,1);
%% 
A3=[zeros(n) B0;B0' zeros(m)];
G3=graph(A3);
c1=n*e2(2); 
c2=1/(1/m-1/n)*(e1(2)-e2(2)); 
%% cvx
cvx_begin quiet sdp
    %cvx_solver mosek
%     cvx_solver sedumi
    % cvx_solver SDPT3

    variables zeta X(N,N);
   
    maximize( c*zeta - trace(L0'*X) );
    subject to
        for i=1:n
            for j=1:m
            E = zeros(N);    
            E(i,i)=1;    
            E(j+n,j+n)=1;    
            E(i,j+n)=-1;    
            E(j+n,i)=-1;
            trace(E'*X)<=-zeta;
            end
        end
        X>=0;
        trace(I'*X) ==1;
        trace(ee'*X) ==0;
        
 cvx_end
 %% Embedding problem
 [P,D]=eig(X);
 U=D^.5*P';  %X=U'*U
 %% Plots
% figure;
% g1=plot(G1,'NodeColor', 'b', 'EdgeColor', [0.80 0.80 0.80],'MarkerSize', 5); 
% title('G1'); 
% figure;
% g2=plot(G2,'NodeColor', 'r', 'EdgeColor', [0.80 0.80 0.80]); 
% title('G2'); 
% labelnode(g1,1:numnodes(G1),''); 
% labelnode(g2,1:numnodes(G2),'');

x1=U(N,1:n)'; 
x2=U(N,n+1:N)'; 
y1=U(N-1,1:n)'; 
y2=U(N-1,n+1:N)'; 
z1=U(N-2,1:n)'; 
z2=U(N-2,n+1:N)';

figure; 
h=plot(G1,'XData',x1,'YData',y1,'ZData',z1,'NodeColor', [0.47,0.67,0.19], 'EdgeColor', [0.47,0.67,0.19], 'MarkerSize', 4, 'LineWidth', .5);
labelnode(h,1:numnodes(G1),'')
hold on; 
g=plot(G2,'XData',x2,'YData',y2,'ZData',z2,'NodeColor', [0.93,0.69,0.13], 'EdgeColor', [0.93,0.69,0.13], 'MarkerSize', 4, 'LineWidth', .5);
labelnode(g,1:numnodes(G2),'');

x = [x1;x2]; 
y=[y1;y2]; 
z=[z1;z2];
hold on; 
f=plot(G3,'-','XData',x,'YData',y,'ZData',z,'NodeColor', [0.80 0.80 0.80], 'EdgeColor', [0.80 0.80 0.80], 'MarkerSize', 1, 'LineWidth', .1)
labelnode(f,1:numnodes(G3),'');
legend('Graph G_1','Graph G_2','Interlayer graph')
%% 
figure; 
h=plot(G1,'XData',x1,'YData',y1,'NodeColor', [0.47,0.67,0.19], 'EdgeColor', [0.47,0.67,0.19], 'MarkerSize', 5, 'LineWidth', .5);
labelnode(h,1:numnodes(G1),'')
hold on; 
g=plot(G2,'XData',x2,'YData',y2,'NodeColor', [0.93,0.69,0.13], 'EdgeColor', [0.93,0.69,0.13], 'MarkerSize', 6, 'LineWidth', .5);
labelnode(g,1:numnodes(G2),'');
x = [x1;x2]; 
y=[y1;y2]; 
hold on; 
f=plot(G3,'-','XData',x,'YData',y,'NodeColor', [0.80 0.80 0.80], 'EdgeColor', [0.80 0.80 0.80], 'MarkerSize', 1, 'LineWidth', .1);
labelnode(f,1:numnodes(G3),'');
legend('Graph G_1','Graph G_2','Interlayer graph')
if c==1; 
    ylim([-0.1 0.1]); 
end