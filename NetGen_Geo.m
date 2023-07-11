function adj=NetGen_Geo(N,r)

x=rand(1,N); y=rand(1,N);

r2=r^2;

l=0;
for i=1:N
    for j=i+1:N
        d2=(x(i)-x(j))^2+(y(i)-y(j))^2;
        if d2<=r2
            l=l+1;
            L1(l)=i; L2(l)=j;
        end
    end
end

adj=sparse(L1,L2,1,N,N); adj=adj+adj'; 


end