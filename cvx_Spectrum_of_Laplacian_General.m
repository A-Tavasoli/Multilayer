clc;clear
%Optimization input
% n=30; m=30; N=n+m; B0=zeros(n,m);  nm=n*m; 
load('A1GeoCase1'); load('A2GeoCase1'); n=max(size(A1)); m=max(size(A2)); N=n+m; nm=n*m; B0=zeros(n,m);
G1 = graph(A1); L0_1=laplacian(G1); G2 = graph(A2); L0_2=laplacian(G2);
e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
[e1(2) e2(2)],  [e1(2)/n e2(2)/m],
%% 
% load('B0RandK3_2'); 
K=3; nc=10;
Nodes1=randi([1,n],nc,1);  Nodes2=randi([1,m],nc,1);
Links=[Nodes1 Nodes2];
for i=1:nc-1
    for j=i+1:nc
        if Links(i,1)==Links(j,1) && Links(i,2)==Links(j,2)
            nc=nc-1;
        end
    end
end
nc
for i=1:n
    for j=1:m
        l=0;
        for k=1:nc
            if Links(k,1)==i && Links(k,2)==j
               B0(i,j)=1;
            end
        end
    end
end
% 
%% 
% K=5; nc=K*n;
% for i=1:K, B0(1,i)=1; end
% for i=2:n
%     B0(i,:)=circshift(B0(i-1,:),1);
%% 
% end
%% 
% load('A1'); load('A2'); n=max(size(A1)); m=max(size(A2)); N=n+m; nm=n*m; 
% G1 = graph(A1); L0_1=laplacian(G1); G2 = graph(A2); L0_2=laplacian(G2);
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)]
%% BA Layers
% n1=15; seed1 = rand(n1,n1) < 1; seed1 = triu(seed1,1); seed1 = seed1 + seed1';
% n2=6; seed2 = rand(n2,n2) < 1; seed2 = triu(seed2,1); seed2 = seed2 + seed2';
% mlink1=2; mlink2=1; A1 = SFNG(n, mlink1, seed1); A2 = SFNG(m, mlink2, seed2);
% G1 = graph(A1); L0_1=laplacian(G1);
% G2 = graph(A2); L0_2=laplacian(G2);
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)],  [e1(2)/n e2(2)/m],
%% ER layers
% p1=.12; p2=.1; G1=ErdosRenyi(n,p1); G2=ErdosRenyi(m,p2); 
% L0_1=laplacian(G1); L0_2=laplacian(G2);
% A1=adjacency(G1); A2=adjacency(G2); 
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)], [e1(2)/n e2(2)/m],
%% Geo layers
% r1=.8*sqrt(2*log(n)/n);    
% r2=.6*sqrt(2*log(m)/m); 
% A1=NetGen_Geo(n,r1);
% A2 =NetGen_Geo(m,r2);
% G1 = graph(A1); L0_1=laplacian(G1);
% G2 = graph(A2); L0_2=laplacian(G2);
% e1 = eig(L0_1); e2 = eig(L0_2); [eigen0_1,e] = eig(full(L0_1)); [eigen0_2,e] = eig(full(L0_2));
% [e1(2) e2(2)],
% [e1(2)/n e2(2)/m],
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
c1=n*e2(2); c2=1/(1/m-1/n)*(e1(2)-e2(2));
cf=50; dc=1; i=0; C=zeros(cf/dc+1,1);  lam_u=zeros(cf/dc+1,N);  lam=zeros(cf/dc+1,N); 
W=zeros(cf/dc+1,nm); WW=cell(cf/dc+1,1); wSum1=zeros(cf/dc+1,n); wSum2=zeros(cf/dc+1,m); 
lam_u1=zeros(cf/dc+1,N);
for c=0:dc:cf
    i=i+1,
    cvx_begin quiet sdp
        cvx_solver mosek
        % cvx_solver sedumi
        % cvx_solver SDPT3

        variables w(n,m) nu l2;
        expression L(N,N); 
        L=[diag(B0.*w*um) -B0.*w;-B0'.*w' diag(B0'.*w'*un)];
        maximize(l2);
        subject to
            L + L0 +nu*ones(N,N)-l2*eye(N)>=0;
            sum(sum(w)) ==c; 
            
            vec(w) >=0;
     cvx_end
     lam(i,:)=eig(L+L0)';
     W(i,:)=vec(w);
     wSum1(i,:)=sum(w,2);
     wSum2(i,:)=sum(w);
     WW{i}=w; C(i)=c;
     Bu=(c/nc)*B0;
     Lu=[diag(Bu*um) -Bu;-Bu' diag(Bu'*un)];
     lam_u(i,:)=eig(Lu+L0)';
     [x,y]=max(eigen0_1(:,2).^2);
     Wu1=zeros(n,m); Wu1(y,:)=c/m;
     Lu1=[diag(Wu1*um) -Wu1;-Wu1' diag(Wu1'*un)];
     lam_u1(i,:)=eig(Lu1+L0)';
end

%% 
figure; hold on;  plot(C,lam(:,2),C,lam_u(:,2),'--','linewidth',2,'Color',[1 0 0]);
b0=(1/n+1/m)*C; plot(C,b0,'linewidth',1);
% for j=3:N
%     plot(C,lam(:,j),'linewidth',2);
% end
figure; hold on; B0_vec=vec(B0);
for j=1:nm
    if B0_vec(j)==1
        plot(C,W(:,j),'linewidth',.1,'Color',[.8 .8 .8]);
    end
end
figure; hold on; 
for j=1:n
    plot(C,wSum1(:,j),'linewidth',.1,'Color',[.8 .8 .8]);
end
figure; hold on; 
for j=1:m
    plot(C,wSum2(:,j),'linewidth',.1,'Color',[.8 .8 .8]);
end
figure; hold on; 
for j=1:4
    plot(C,lam_u(:,j),'linewidth',.1,'Color',[.8 .8 .8]);
end
%%
% c=30; Bu=(c/nm)*B0; Lu=[diag(Bu*um) -Bu;-Bu' diag(Bu'*un)]; Lu=full(Lu+L0); [uu,vu]=eig(Lu); uu(:,2);
% w=WW{20}; w=sum(w,2); [u,v]=eig(full(L0_1)); figure; plot(u(:,2).^2,w,'o')
w=WW{50}; L=[diag(w*um) -w;-w' diag(w'*un)]; L=full(L+L0); [u,v]=eig(full(L));  u=u(:,2); %u=u(1:n);
figure; plot(abs(eigen0_1(:,2)),sum(w,2),'o','Color',[.3 .3 .3])
figure; plot(abs(eigen0_2(:,2)),sum(w,1),'o','Color',[.3 .3 .3])
%% 
w=WW{10}; L=[diag(w*um) -w;-w' diag(w'*un)]; L=full(L+L0);
x0=[ones(n,1);-ones(m,1)]; tspan=[0 2]; x0=rand(N,1);
[t,x] = ode15s(@(t,x) -L*x, tspan, x0); x1=x(:,1:n); x2=x(:,n+1:N);
figure; hold on; plot(t,x1,'r');  plot(t,x2,'b');