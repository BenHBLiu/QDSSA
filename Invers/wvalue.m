function w = wvalue(d,x)
n=length(x);
for k=1:n
    if k<=d
        Jk=1:k;
    else
        Jk=k-d:k;
    end
    w(k)=0;
    m=length(Jk);
    for i=1:m
        pw(i)=1;
        for j=Jk(i):Jk(i)+d
            if (j~=k & j<=n)
                pw(i)=pw(i)*1/(x(k)-x(j));
            end
            if j>n
                pw(i)=0;
            end
        end
        w(k)=w(k)+(-1)^(Jk(i)+1)*pw(i);
    end
end