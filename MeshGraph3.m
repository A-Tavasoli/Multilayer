function G=MeshGraph3(n)

A=zeros(n);
for i=1:3:n
    for j=0:1
        A(i+j,i+j+1)=1;
        A(i+j+1,i+j)=1;
    end
end
for i=1:3:n-5
    for j=0:2
        A(i+j,i+j+3)=1;
        A(i+j+3,i+j)=1;
    end
end
G=graph(A);