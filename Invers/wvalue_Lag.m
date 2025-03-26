function w = wvalue_Lag(x)
n=length(x);
for j=1:n
    w(j)=1;
    for k=1:n
        if j~=k
            w(j)=w(j)/(x(j)-x(k));
        end
    end
end
w=w/abs(w(1));
end