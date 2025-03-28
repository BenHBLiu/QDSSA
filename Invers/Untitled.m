clc
clear
a=0;b=1;%x的范围
epsilon=1e-8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%修改%%%%%%%%%%%%%%%%%%%%


q1 = [2.98872459409823,5.94241883208608,4.99962065388960,1.28902592574191,1.21093700545359,4.06732481471903,5.21820876733293,1.34956939479910,5.51368432658255
    ]; 
q1=sort(q1);

n=length(q1) - 1;
m=n*2;
dx=1;

lambda1=-1;
lambda2=1;
mu=3*epsilon;

%Chebyshev nodes
x=-cos((0:n)*pi/n)';
%equidistant nodes
tau=1/m;
t=0:tau:1;

for i=1:length(x)
    if x(i)<0
        x(i)=0.5*(lambda1+mu*sinh(asinh((1+lambda1)/mu)*x(i)+asinh((1-lambda1)/mu)*(x(i)+1))-1);
    else
        x(i)=0.5*(lambda2+mu*sinh(asinh((1+lambda2)/mu)*(x(i)-1)+asinh((1-lambda2)/mu)*x(i))+1);
    end
end

x=0.5*(b+a)+0.5*(b-a)*x;
x(1)=0;

n=length(x);m=length(t);

%Barycentric rational intepolation
D=baryrat(x,dx,2);D1=D(:,:,1);D2=D(:,:,2);

I=eye(n,n);




q=zeros(1,n);
for i=1:length(x)
    q(i)=qun(x(i));
end

figure(3)
subplot(121)
plot(x,q1,'b*-');
legend('Numerical q','location','northwest')

subplot(122)
plot(x,q,'rp-');
legend('Exact q','location','northwest')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L=-epsilon*D2+diag(q);

F=zeros(n,m);
for i=1:n
    for j=1:m
        % F(i,j)=exp(-t(j))*((exp(-x(i)/epsilon^(1/2)) + exp((x(i) - 1)/epsilon^(1/2)))/(exp(-1/epsilon^(1/2)) + 1) - cos(pi*x(i))^2) + epsilon*(exp(-t(j))...
        %     - 1)*((exp(-x(i)/epsilon^(1/2))/epsilon + exp((x(i) - 1)/epsilon^(1/2))/epsilon)/(exp(-1/epsilon^(1/2)) + 1) - 2*pi^2*sin(pi*x(i))^2 + 2*pi^2*cos(pi*x(i))^2)...
        %     - ((exp(-x(i)/epsilon^(1/2)) + exp((x(i) - 1)/epsilon^(1/2)))/(exp(-1/epsilon^(1/2)) + 1) - cos(pi*x(i))^2)*(exp(-t(j)) - 1)*(x(i) + 1);
        F(i,j)=exp(-t(j))*((exp(-x(i)/epsilon^(1/2))+exp((x(i)-1)/epsilon^(1/2)))/(exp(-1/epsilon^(1/2))+1)-cos(pi*x(i))^2)+...
    epsilon*(exp(-t(j))-1)*((exp(-x(i)/epsilon^(1/2))/epsilon+exp((x(i)-1)/epsilon^(1/2))/epsilon)/(exp(-1/epsilon^(1/2))+1)-2*pi^2*sin(pi*x(i))^2+2*pi^2*cos(pi*x(i))^2)...
            -((exp(-x(i)/epsilon^(1/2))+exp((x(i)-1)/epsilon^(1/2)))/(exp(-1/epsilon^(1/2))+1)-cos(pi*x(i))^2)*(exp(-t(j))-1)*(3*x(i)^2+2*x(i)+1);
    end
end


ue=zeros(n,m);
for i=1:n
    for j=1:m
        ue(i,j)=(1-exp(-t(j)))*((exp(-x(i)/sqrt(epsilon))+exp(-(1-x(i))/sqrt(epsilon)))/(1+exp(-1/sqrt(epsilon)))-cos(pi*x(i))^2);
    end
end

ue_xx=zeros(n,m);
for i=1:n
    for j=1:m
        ue_xx(i,j)=-(exp(-t(j))-1)*((exp(-x(i)/epsilon^(1/2))/epsilon + exp((x(i)-1)/epsilon^(1/2))/epsilon)/(exp(-1/epsilon^(1/2))+1)-2*pi^2*sin(pi*x(i))^2+2*pi^2*cos(pi*x(i))^2);
    end
end

A=1/tau*I+0.5.*L;
B=1/tau*I-0.5.*L;
A(1,:)=I(1,:);
A(:,1)=I(:,1);
A(n,:)=I(n,:);
A(:,n)=I(:,n);

uc=zeros(n,m);
for j=1:m-1
    C=B*uc(:,j)+0.5*(F(:,j)+F(:,j+1));
    C(1)=0;
    C(n)=0;
    uc(:,j+1)=A\C;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%修改%%%%%%%%%%%%%%%%%%%%
% o1=sum(((uc(2:end-1,end)-ue(2:end-1,end))).^2.*(x(2:end-1)-x(1:end-2)));
o1=max(abs(uc(2:end-1,end)-ue(2:end-1,end)));
% o2=0;
% for i=2:n-1
% o2=((F(i,end)-(uc(i,end)-uc(i,end-1))/tau+epsilon*ue_xx(i,end)-q(i)'.*ue(i,end)).^2).*(x(i)-x(i-1)) +o2;
% end
% o=o1+0.5*o2;
qexact=zeros(1,n);
for i=1:length(x)
    qexact(i)=qun(x(i));
end
o2=sqrt(sum((q1-qexact).^2));
% o2=max(abs(q1-qexact));
o=o1+o2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%修改%%%%%%%%%%%%%%%%%%%%


% 创建画布和子图
figure;
% 绘制第一个三维图
subplot(1,2,1);
mesh(t,x,uc);
title('Numerical solution');

% 绘制第二个三维图
subplot(1,2,2);
mesh(t,x,ue);
title('Exact solution');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%修改%%%%%%%%%%%%%%%%%%%%
% 绘制误差图
figure(2)
mesh(t,x,abs(uc-ue));
title('Error');


sum(abs(uc(:,end)-ue(:,end)))

% N=64       1.651391758105714e-02       1.651543730103344e-02
% N=128     3.434612706724294e-03       3.434955027127984e-03
% N=256     6.855466888688877e-04       6.855466888688877e-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%修改%%%%%%%%%%%%%%%%%%%%