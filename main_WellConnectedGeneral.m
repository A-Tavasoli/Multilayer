clc;clear
%Optimization input
n=6; %number of G_1 nodes
m=6; %number of G_2 nodes
cf=15; %final total budget 
mc = 6;
q=cf/mc; %mc number of intrlinks
B0=ones(n,m);
load('A_Small1'); 
load('A_Small2'); 
n=max(size(A1)); 
m=max(size(A2)); 
N=n+m; nm=n*m; 
B0=ones(n,m);
G1 = graph(A1); 
L0_1=laplacian(G1); 
G2 = graph(A2); 
L0_2=laplacian(G2);
%% BA Layers
% n1=3; 
% seed1 = rand(n1,n1) < 1; 
% seed1 = triu(seed1,1); 
% seed1 = seed1 + seed1';
% n2=2; 
% seed2 = rand(n2,n2) < 1; 
% seed2 = triu(seed2,1); 
% seed2 = seed2 + seed2';
% mlink1=1; 
% mlink2=1; 
% A1 = SFNG(n, mlink1, seed1); 
% A2 = SFNG(m, mlink2, seed2);
% G1 = graph(A1); L0_1=laplacian(G1);
% G2 = graph(A2); L0_2=laplacian(G2);
%% ER layers
% p1=.3; 
% p2=.5; 
% G1=ErdosRenyi(n,p1); 
% G2=ErdosRenyi(m,p2); 
% L0_1=laplacian(G1); 
% L0_2=laplacian(G2);
% A1=adjacency(G1);
% A2=adjacency(G2);
%% Geo layers
% r1=1*sqrt(2*log(n)/n);    
% r2=.6*sqrt(2*log(m)/m); 
% A1=NetGen_Geo(n,r1);
% A2 =NetGen_Geo(m,r2);
% G1 = graph(A1); 
% L0_1=laplacian(G1);
% G2 = graph(A2); 
% L0_2=laplacian(G2);
%% WS layers
% K1=2; 
% K2=3; 
% beta1=.2; 
% beta2=.2; 
% G1 = WattsStrogatz(n,K1,beta1); 
% L0_1=laplacian(G1);
% G2 = WattsStrogatz(m,K2,beta2); 
% L0_2=laplacian(G2);
% A1=adjacency(G1); 
% A2=adjacency(G2); 
%% Combined layers
% r2=.8*sqrt(2*log(m)/m); 
% A2=NetGen_Geo(n,r2); 
% G2 = graph(A2); 
% L0_2=laplacian(G2); 
% K1=1; beta1=0; 
% G1 = WattsStrogatz(n,K1,beta1); 
% L0_1=laplacian(G1); 
% A1=adjacency(G1);  
% G1=WheelGraph(n); 
% L0_1=laplacian(G1); 
% A1=adjacency(G1);
% G1=MeshGraph3(n); 
% L0_1=laplacian(G1); 
% A1=adjacency(G1);
% load('A_Small'); 
% A2=A_Small; 
% G2 = graph(A2); 
% L0_2=laplacian(G2);
% p2=.4; 
% G2=ErdosRenyi(m,p2);   
% L0_2=laplacian(G2);   
% A2=adjacency(G2);
% n2=4; 
% seed2 = rand(n2,n2) < 1; 
% seed2 = triu(seed2,1); 
% seed2 = seed2 + seed2';
% mlink2=2; 
% A2 = SFNG(m, mlink2, seed2); 
% G2 = graph(A2); 
% L0_2=laplacian(G2);
%% Computing the laplacian L0
L0=[L0_1 zeros(n,m);zeros(m,n) L0_2];
% computing ee'
if n==m
    L_ave=(L0_1+L0_2)/2; 
    e_ave=eig(L_ave); 
    lam_ave=e_ave(2);
end
ee=ones(N,N); 
un=ones(n,1); 
um=ones(m,1);
%% 
[EigVec1,eig1]=eig(full(L0_1)); [EigVec2,eig2]=eig(full(L0_2));
A=zeros(N);
[cent1,CentNode1]=max(centrality(G1,'degree')); 
[cent2,CentNode2]=max(centrality(G2,'degree')); 
% [cent1,CentNode1]=max(abs(EigVec1(:,2))); 
% [cent2,CentNode2]=min(abs(EigVec2(:,2))); 
a=zeros(N,1); 
a(CentNode1)=1; 
a(CentNode2+n)=-1;
L=L0+q*a*a'; 
[u,v]=eig(full(L)); 
u=u(:,2); 
lam2=zeros(mc,1); 
lam2(1)=v(2,2);
U=zeros(n,n); 
A=A+a*a';
%% 
e1 = eig(L0_1); 
e2 = eig(L0_2); 
c1=n*e2(2); 
c2=1/(1/m-1/n)*(e1(2)-e2(2));
em = max([e1(2) e2(2)]);
%% 
for k=2:mc
    [u,v]=eig(full(L)); 
    u=u(:,2); 
    lam2(k)=v(2,2);
    for i=1:n
        for j=1:m
            U(i,j)=(u(i)-u(j+n))^2;
        end
    end
    [p,i]=max(max(U'));  
    [p,j]=max(max(U));
    a=zeros(N,1); 
    a(i)=1; 
    a(j+n)=-1;
    L=L+q*a*a';
    A=A+a*a';
end
% figure; plot(lam2)
%% 
Au=[eye(n) -eye(n);-eye(n) eye(n)];
i=0; dc=1; lam_unif=zeros(cf/dc+1,N); lam_well=zeros(cf/dc+1,N); C=zeros(cf/dc+1,1);
for c=0:dc:cf
    i=i+1;
    q=c/n;
    Lw=L0+q*A;  
    lam_well(i,:)=eig(Lw); 
    C(i)=c;
end
%% 
Gw=graph(A-diag(diag(A))); 
Gu=graph(Au-diag(diag(Au)));
figure; 
p1=plot(G1); 
x1=p1.XData; 
y1=p1.YData; 
p2=plot(G2); 
x2=p2.XData; 
y2=p2.YData; 
z1=ones(n,1); 
z2=-ones(m,1);
%=====================================================================================================================
Color1=[0 0.45 0.74]; 
Color2=[0.47 0.67 0.19]; 
Color3='r'; 
h=plot(G1,'XData',x1,'YData',y1,'ZData',z1,'NodeColor', Color1, 'EdgeColor', Color1, 'MarkerSize', 8, 'LineWidth', 2);
% labelnode(h,1:numnodes(G1),'')
hold on; 
g=plot(G2,'XData',x2,'YData',y2,'ZData',z2,'NodeColor', Color2, 'EdgeColor', Color2, 'MarkerSize', 8, 'LineWidth', 2);
% labelnode(g,1:numnodes(G2),'');
x = [x1';x2']; 
y=[y1';y2']; 
z=[z1;z2];
hold on; 
f=plot(Gw,'-','XData',x,'YData',y,'ZData',z,'NodeColor', Color3, 'EdgeColor', Color3, 'MarkerSize', 4, 'LineWidth', 2)
labelnode(f,1:numnodes(Gw),'');
view(-107.1,39.6)
%=====================================================================================================================
AA1=diag(A(1:n,1:n));
AA2=diag(A(1+n:N,n+1:N));
figure; 
plot(EigVec1(:,2),AA1,'o','Color',[.3 .3 .3])
figure; 
plot(EigVec2(:,2),AA2,'o','Color',[.3 .3 .3])
%=====================================================================================================================
%% 
Au=[eye(n) -eye(n);-eye(n) eye(n)];
i=0; dc=1; 
lam_unif=zeros(cf/dc+1,N); 
lam_well=zeros(cf/dc+1,N); 
C=zeros(cf/dc+1,1);
for c=0:dc:cf
    i=i+1;
    q=c/n;
    Lu=L0+q*Au; 
    lam_unif(i,:)=eig(Lu);
    Lw=L0+q*A;  
    lam_well(i,:)=eig(Lw); 
    C(i)=c;
end
figure; 
plot(C,lam_well(:,2),'k',C,lam_unif(:,2),'r','linewidth',2)
hold on; 
plot([0 cf],[em em],'-.',[0 cf],[lam_ave lam_ave],'-.','linewidth',2)
legend('Well-interconnected','One-to-one','max(\lambda_2[L_1],\lambda_2[L_2])','\lambda_2[L_{ave}]')
