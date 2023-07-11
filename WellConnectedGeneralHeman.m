clc;clear
%Optimization input
n=5; nGrid=3; mGrid=3; m=nGrid*mGrid;
N=n+m; nm=n*m;   mc=7; cf=mc; q=cf/mc; %mc number of intrlinks
B0=ones(n,m);
%% Geo layer
r1=.6*sqrt(2*log(n)/n);    
A1=NetGen_Geo(n,r1);
G1 = graph(A1); L0_1=laplacian(G1);
%% Grid layer
GridEdges = load('grid_3_3.txt');
sGrid = GridEdges(:,1);
tGrid = GridEdges(:,2);
G2 = graph(sGrid,tGrid);
L0_2=laplacian(G2);
%% Computing the laplacian L0
L0=[L0_1 zeros(n,m);zeros(m,n) L0_2];
% computing ee'
ee=ones(N,N); un=ones(n,1); um=ones(m,1);
%% 
[EigVec1,eig1]=eig(full(L0_1)); [EigVec2,eig2]=eig(full(L0_2));
A=zeros(N);
[cent1,CentNode1]=max(centrality(G1,'degree')); 
[cent2,CentNode2]=max(centrality(G2,'degree')); 
% [cent1,CentNode1]=max(abs(EigVec1(:,2))); 
% [cent2,CentNode2]=min(abs(EigVec2(:,2))); 
a=zeros(N,1); a(CentNode1)=1; a(CentNode2+n)=-1;
L=L0+q*a*a'; [u,v]=eig(full(L)); u=u(:,2); lam2=zeros(mc,1); lam2(1)=v(2,2);
U=zeros(n,n); A=A+a*a';
%% 
for k=2:mc
    [u,v]=eig(full(L)); u=u(:,2); lam2(k)=v(2,2);
    for i=1:n
        for j=1:m
            if L(j+n,1:n)==0
                U(i,j)=(u(i)-u(j+n))^2;
            else
                U(i,j)=0;
            end
        end
    end
    [p,i]=max(max(U'));  
    [p,j]=max(max(U));
    a=zeros(N,1); a(i)=1; a(j+n)=-1;
    L=L+q*a*a';
    A=A+a*a';
end
% figure; plot(lam2)
%% 
Gw=graph(diag(diag(A))-A); 
figure; p1=plot(G1); x1=p1.XData; y1=p1.YData; p2=plot(G2); x2=p2.XData; y2=p2.YData; z1=ones(n,1); z2=-ones(m,1);
%=====================================================================================================================
Color1=[0 0.45 0.74]; Color2=[0.47 0.67 0.19]; Color3='r'; 
h=plot(G1,'XData',x1,'YData',y1,'ZData',z1,'NodeColor', Color1, 'EdgeColor', Color1, 'MarkerSize', 8, 'LineWidth', 2);
% labelnode(h,1:numnodes(G1),'')
hold on; g=plot(G2,'XData',x2,'YData',y2,'ZData',z2,'NodeColor', Color2, 'EdgeColor', Color2, 'MarkerSize', 8, 'LineWidth', 2);
% labelnode(g,1:numnodes(G2),'');
x = [x1';x2']; y=[y1';y2']; z=[z1;z2];
hold on; f=plot(Gw,'-','XData',x,'YData',y,'ZData',z,'NodeColor', Color3, 'EdgeColor', Color3, 'MarkerSize', 4, 'LineWidth', 2);
labelnode(f,1:numnodes(Gw),'');
%=====================================================================================================================
AA1=diag(A(1:n,1:n));
AA2=diag(A(1+n:N,n+1:N));
figure; plot(EigVec1(:,2),AA1,'o','Color',[.3 .3 .3])
figure; plot(EigVec2(:,2),AA2,'o','Color',[.3 .3 .3])
%=====================================================================================================================
%% Generating random connections
L_rand = L0;
A_rand=zeros(N);
i_rand = [ randperm(n,n) randi(n,[1,mc-n]) ];
j_rand = randperm(m,mc);
for k = 1:mc
    a_rand=zeros(N,1); 
    a_rand(i_rand(k))=1; 
    a_rand(j_rand(k)+n)=-1;
    L_rand=L_rand+q*a_rand*a_rand';
    A_rand=A_rand+a_rand*a_rand';
end
Gw_rand=graph(A_rand-diag(diag(A_rand)));
figure; 
p1=plot(G1); x1=p1.XData; y1=p1.YData; 
p2=plot(G2); x2=p2.XData; y2=p2.YData; z1=ones(n,1); z2=-ones(m,1);
Color1=[0 0.45 0.74]; Color2=[0.47 0.67 0.19]; Color3='r'; 
h=plot(G1,'XData',x1,'YData',y1,'ZData',z1,'NodeColor', Color1, 'EdgeColor', Color1, 'MarkerSize', 8, 'LineWidth', 2);
% labelnode(h,1:numnodes(G1),'')
hold on; g=plot(G2,'XData',x2,'YData',y2,'ZData',z2,'NodeColor', Color2, 'EdgeColor', Color2, 'MarkerSize', 8, 'LineWidth', 2);
% labelnode(g,1:numnodes(G2),'');
x = [x1';x2']; y=[y1';y2']; z=[z1;z2];
hold on; f=plot(Gw_rand,'-','XData',x,'YData',y,'ZData',z,'NodeColor', Color3, 'EdgeColor', Color3, 'MarkerSize', 4, 'LineWidth', 2);
labelnode(f,1:numnodes(Gw_rand),'');

e1 = eig(L0_1); e2 = eig(L0_2); 
[e1(2) e2(2)],
eig_opt  = eig(L);    
eig_rand = eig(L_rand);
eig_opt(2),
eig_rand(2),