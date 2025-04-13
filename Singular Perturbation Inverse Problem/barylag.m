%重心Lagrange插值微分矩阵
function D=barylag(x,m)
n=length(x);
D=zeros(n,n,m);
b=wvalue_Lag(x);
for i=1:n
    for j=1:n
        if i~=j
            D(i,j,1)=(b(j)/b(i))/(x(i)-x(j));
        end
    end
end
for i=1:n
    D(i,i,1)=-sum(D(i,:,1));
end
for k=2:m
    for i=1:n
        for j=1:n
            if i~=j
                D(i,j,k)=k*(D(i,i,k-1)*D(i,j,1)-D(i,j,k-1)/(x(i)-x(j)));
            end
        end
    end
    for i=1:n
        D(i,i,k)=-sum(D(i,:,k));
    end
end
end