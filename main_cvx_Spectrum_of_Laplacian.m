clc;clear
%Optimization input
% n=30; 
% m=10; 
% N=n+m; 
% nm=n*m; 
% B0=ones(n,m);
load('A1GeoCase1'); 
load('A2GeoCase1'); 
n=max(size(A1)); 
m=max(size(A2)); 
N=n+m; nm=n*m; 
B0=ones(n,m);
G1 = graph(A1); 
L0_1=laplacian(G1); 
G2 = graph(A2); 
L0_2=laplacian(G2);
e1 = eig(L0_1); 
e2 = eig(L0_2); 
[eigen0_1,e] = eig(full(L0_1)); 
[eigen0_2,e] = eig(full(L0_2));
%% BA Layers
% n1=15; 
% seed1 = rand(n1,n1) < 1; 
% seed1 = triu(seed1,1); 
% seed1 = seed1 + seed1';
% n2=6; 
% seed2 = rand(n2,n2) < 1; 
% seed2 = triu(seed2,1); 
% seed2 = seed2 + seed2';
% mlink1=2; mlink2=1; 
% A1 = SFNG(n, mlink1, seed1); 
% A2 = SFNG(m, mlink2, seed2);
% G1 = graph(A1); 
% L0_1=laplacian(G1);
% G2 = graph(A2); 
% L0_2=laplacian(G2);
% e1 = eig(L0_1); 
% e2 = eig(L0_2); 
% [eigen0_1,e] = eig(full(L0_1)); 
% [eigen0_2,e] = eig(full(L0_2));
%% ER layers
% p1=.2; 
% p2=.15; 
% G1=ErdosRenyi(n,p1); 
% G2=ErdosRenyi(m,p2); 
% L0_1=laplacian(G1); 
% L0_2=laplacian(G2);
% A1=adjacency(G1); 
% A2=adjacency(G2); 
% e1 = eig(L0_1); 
% e2 = eig(L0_2); 
% [eigen0_1,e] = eig(full(L0_1)); 
% [eigen0_2,e] = eig(full(L0_2));
%% Geo layers
% r1=.8*sqrt(2*log(n)/n);    
% r2=.7*sqrt(2*log(m)/m); 
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
% K1=2; 
% beta1=.5; 
% K2=1;  
% beta2=1; 
% G1 = WattsStrogatz(n,K1,beta1); 
% L0_1=laplacian(G1);
% G2 = WattsStrogatz(m,K2,beta2); 
% L0_2=laplacian(G2);
% A1=adjacency(G1); 
% A2=adjacency(G2); 
% e1 = eig(L0_1); 
% e2 = eig(L0_2); 
% [eigen0_1,e] = eig(full(L0_1)); 
% [eigen0_2,e] = eig(full(L0_2));
%% Combined layers
% r1=.8*sqrt(2*log(n)/n); 
% A1=NetGen_Geo(n,r1);   
% G1 = graph(A1); 
% L0_1=laplacian(G1);
% p2=.3; 
% G2=ErdosRenyi(m,p2);   
% L0_2=laplacian(G2);   
% A2=adjacency(G2);
% e1 = eig(L0_1); 
% e2 = eig(L0_2); 
% [eigen0_1,e] = eig(full(L0_1)); 
% [eigen0_2,e] = eig(full(L0_2));
%% Computing the laplacian L0
L0=[L0_1 zeros(n,m);zeros(m,n) L0_2];
% computing ee'
ee=ones(N,N); 
un=ones(n,1); 
um=ones(m,1);
%% cvx
c1=n*e2(2); 
c2=1/(1/m-1/n)*(e1(2)-e2(2)); 
cf=50; 
dc=1; 
i=0; 
C=zeros(cf/dc+1,1);  
lam_u=zeros(cf/dc+1,N);  
lam=zeros(cf/dc+1,N); b1=zeros(cf/dc+1,1);
W=zeros(cf/dc+1,nm); 
WW=cell(cf/dc+1,1); 
wSum1=zeros(cf/dc+1,n); 
wSum2=zeros(cf/dc+1,m); 
b2=zeros(cf/dc+1,1);
lam_u1=zeros(cf/dc+1,N);
for c=0:dc:cf
    i=i+1,
    cvx_begin quiet sdp
        %cvx_solver mosek
        % cvx_solver sedumi
        % cvx_solver SDPT3

        variables w(n,m) nu l2;
        expression L(N,N); 
        L=[diag(w*um) -w;-w' diag(w'*un)];
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
     b1(i)=e1(2)+eigen0_1(:,2)'.^2*sum(w,2);
     b2(i)=e2(2)+eigen0_2(:,2)'.^2*sum(w)';
     WW{i}=w; C(i)=c;
     Bu=(c/nm)*B0; Lu=[diag(Bu*um) -Bu;-Bu' diag(Bu'*un)];
     lam_u(i,:)=eig(Lu+L0)';
     [x,y]=max(eigen0_1(:,2).^2);
     Wu1=zeros(n,m); Wu1(y,:)=c/m;
     Lu1=[diag(Wu1*um) -Wu1;-Wu1' diag(Wu1'*un)];
     lam_u1(i,:)=eig(Lu1+L0)';
end

%% 
figure; 
hold on;  
plot(C,lam(:,2),'k',C,lam_u(:,2),'k--','linewidth',2);
xlabel('c'); 
ylabel('\lambda_2[L]');
legend('Optimal weights', 'Uniform weights')
% figure; hold on;  plot(C,(1/n+1/m)*C,C,e1(2)+C/n,C,e2(2)+C/m,C,lam_u(:,2));
% b0=(1/n+1/m)*C; plot(C,b0,'k',C,b1,'b',C,b2,'g','linewidth',1);
% for j=3:N
%     plot(C,lam(:,j),'linewidth',2);
% end
i=0; 
[u,v]=eig(full(L0_2)); 
v=eig(full(L0_2));
for c=0:dc:cf
    i=i+1;
    aa(i)=v(2)+u(:,2)'*diag(WW{i}'*ones(n,1))*u(:,2);
end
i=0; [u,v]=eig(full(L0_1)); v=eig(full(L0_1));
for c=0:dc:cf
    i=i+1;
    bb(i)=v(2)+u(:,2)'*diag(WW{i}*ones(m,1))*u(:,2);
end  
% figure, plot(C,aa,C,bb,C,lam(:,2),'--','linewidth',2)
% lam0_1=eig(L0_1); 
% lam0_2=eig(L0_2); 
% Bound0=(1/n+1/m)*C; 
% Bound1=sum(lam0_1)/n+C/n; 
% Bound2=sum(lam0_2)/m+C/m;
% plot(C,Bound0,C,Bound1,C,Bound2);
% figure; 
% hold on; 
% for j=1:nm
%     plot(C,W(:,j),'linewidth',1,'Color',[.6 .6 .6]);
% end
figure; hold on; 
xlabel('c'); 
ylabel('W^T1_m')
for j=1:n
    plot(C,wSum1(:,j),'linewidth',1,'Color',[.6 .6 .6]);
end
figure; hold on; 
xlabel('c'); ylabel('W^T1_n')
for j=1:m
    plot(C,wSum2(:,j),'linewidth',1,'Color',[.6 .6 .6]);
end
% figure; hold on; 
% for j=1:4
%     plot(C,lam_u(:,j),'linewidth',1,'Color',[.6 .6 .6]);
% end
%%
% c=5; 
% Bu=(c/nm)*B0; 
% Lu=[diag(Bu*um) -Bu;-Bu' diag(Bu'*un)]; 
% Lu=full(Lu+L0); 
% [uu,vu]=eig(Lu); 
% uu(:,2);
% w=WW{c}; 
% w=sum(w,2); 
% [u,v]=eig(full(L0_1)); 
% figure; 
% plot(u(:,2).^2,w,'o')
% w=WW{c}; 
% L=[diag(w*um) -w;-w' diag(w'*un)]; L=full(L+L0); [u,v]=eig(full(L));  
% u=u(:,2); 
%u=u(1:n);
c=20; w=WW{c};
figure; plot(abs(eigen0_1(:,2)),sum(w,2),'o','Color',[.3 .3 .3]); xlabel('|u_2^{(1)}(i)|'); ylabel('Optimal weights');
figure; plot(abs(eigen0_2(:,2)),sum(w,1),'o','Color',[.3 .3 .3]); xlabel('|v_2^{(2)}(i)|'); ylabel('Optimal weights'); 
xlim([0 0.5]); ylim([0 1.6])
%% 
% w=WW{2}; 
% L=[diag(w*um) -w;-w' diag(w'*un)]; L=full(L+L0);  
% [u,v]=eig(full(L));  
% u=u(:,2),
% x0=[ones(n,1);-ones(m,1)]; 
% tspan=[0 30]; %x0=rand(N,1);
% [t,x] = ode15s(@(t,x) -L*x, tspan, x0); 
% x1=x(:,1:n); 
% x2=x(:,n+1:N);
% figure; 
% hold on; 
% plot(t,x1(:,1),'color',[0.47,0.67,0.19],'linewidth',2); 
% plot(t,x2(:,1),'--','color',[0.93,0.63,0.13],'linewidth',2); 
% plot(t,x1,'color',[0.47,0.67,0.19],'linewidth',2);  
% plot(t,x2,'--','color',[0.93,0.63,0.13],'linewidth',2);
% xlabel('Time'); 
% ylabel('States of G nodes X(t)');
% legend('States of G_1 nodes','States of G_2 nodes')
%% 
% c=10; 
% Bu=(c/nm)*B0; 
% Lu=[diag(Bu*um) -Bu;-Bu' diag(Bu'*un)]; 
% L=full(Lu+L0);
% x0=[zeros(n,1)+1;rand(m,1)]; 
% tspan=[0 20]; 
% x0=[rand(n,1);-rand(m,1)];
% [t,x] = ode15s(@(t,x) -L*x, tspan, x0); 
% x1=x(:,1:n); 
% x2=x(:,n+1:N);
% figure; 
% hold on; 
% plot(t,x1,'r','linewidth',.1);  
% plot(t,x2,'b','linewidth',.1);
%% 
% figure; 
% plot(C,e1(2)+max(eigen0_1(:,2).^2)*C,C,lam(:,2))
% for i=1:n
%     for j=1:m
%         def(i,j)=(eigen0_1(i,2)-eigen0_2(j,2))^2;
%     end
% end
% w=WW{6}; 
% figure, 
% plot(vec(def),vec(w),'o');
% k=0;
% for i=1:n
% for j=1:m
% k=k+1;
% U(i,j)=u1(i)+u2(j);
% end
% end
% figure, 
% plot(vec(WW{15}),vec(U),'o')
% w=WW{100}; 
% [u,v]=eig(full(L0_1)); 
% u21=u(:,2); 
% [u,v]=eig(full(L0_2)); 
% u22=u(:,2);
% for i=1:n
%     for j=1:m
%         def(i,j)=(u21(i)-u22(j))^2;
%     end
% end
% figure, plot(vec(def),vec(w),'o');
% [u,v]=eig(full(A1)); 
% u1=u(:,n)*sign(u(1,n));
% [u,v]=eig(full(A2)); 
% u2=u(:,m)*sign(u(1,m));