function [lb,ub,dim,fobj] = Get_Functions_details(F)
switch F
            
    case 'F29'
        fobj = @F29;
        lb=1;
        ub=6;
        dim=64+1; 
end  
end

function o=F29(q)
a=0;b=1;%x ķ Χ
epsilon=1e-6;

n=length(q) - 1;
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

% for i=1:length(x)
%     if x(i)<0
%         x(i)=0.5*(lambda1+mu*sinh(asinh((1+lambda1)/mu)*x(i)+asinh((1-lambda1)/mu)*(x(i)+1))-1);
%     else
%         x(i)=0.5*(lambda2+mu*sinh(asinh((1+lambda2)/mu)*(x(i)-1)+asinh((1-lambda2)/mu)*x(i))+1);
%     end
% end


asinh1_lambda1 = asinh((1 + lambda1) / mu);
asinh2_lambda1 = asinh((1 - lambda1) / mu);
asinh1_lambda2 = asinh((1 + lambda2) / mu);
asinh2_lambda2 = asinh((1 - lambda2) / mu);


idx_neg = x < 0;
x_neg = x(idx_neg);
x(idx_neg) = 0.5 * (lambda1 + mu * sinh(asinh1_lambda1 * x_neg + asinh2_lambda1 * (x_neg + 1)) - 1);


idx_pos = x >= 0;
x_pos = x(idx_pos);
x(idx_pos) = 0.5 * (lambda2 + mu * sinh(asinh1_lambda2 * (x_pos - 1) + asinh2_lambda2 * x_pos) + 1);

x=0.5*(b+a)+0.5*(b-a)*x;
x(1)=0;

n=length(x);m=length(t);

%Barycentric rational intepolation
D=baryrat(x,dx,2);D1=D(:,:,1);D2=D(:,:,2);

I=eye(n,n);

% q=1+2*x+3*x.^2;
% Q=zeros(1,n);
Q = qun(x)';
% for i=1:length(x)
%     Q(i)=qun(x(i));
% end
L=-epsilon*D2+diag(Q);

%free source
F=zeros(n,m);
% for i=1:n
%     for j=1:m
%         F(i,j)=exp(-t(j))*((exp(-x(i)/epsilon^(1/2))+exp((x(i)-1)/epsilon^(1/2)))/(exp(-1/epsilon^(1/2))+1)-...
%             cos(pi*x(i))^2)+epsilon*(exp(-t(j))-1)*((exp(-x(i)/epsilon^(1/2))/epsilon+exp((x(i)-1)/epsilon^(1/2))/...
%             epsilon)/(exp(-1/epsilon^(1/2))+1)-2*pi^2*sin(pi*x(i))^2+2*pi^2*cos(pi*x(i))^2)-((exp(-x(i)/epsilon^(1/2))+...
%             exp((x(i)-1)/epsilon^(1/2)))/(exp(-1/epsilon^(1/2))+1)-cos(pi*x(i))^2)*(exp(-t(j))-1)*(3*x(i)^2+2*x(i)+1);
%     end
% end


[Xmat,Tmat]=ndgrid(x,t);
F=fun_f(Xmat,Tmat,epsilon);

%Exact solution
% load chi_0.001.mat

% ue=zeros(n,m);
% for i=1:n
%     for j=1:m
%         ue(i,j)=(1-exp(-t(j)))*((exp(-x(i)/sqrt(epsilon))+exp(-(1-x(i))/sqrt(epsilon)))/(1+exp(-1/sqrt(epsilon)))-cos(pi*x(i))^2);
%     end
% end
% ue(:,end)=ue(:,end).*(1+chi)';

ue=fun_ue(Xmat,Tmat,epsilon);

% ue_xx=zeros(n,m);
% for i=1:n
%     for j=1:m
%         ue_xx(i,j)=-(exp(-t(j))-1)*((exp(-x(i)/epsilon^(1/2))/epsilon + exp((x(i)-1)/epsilon^(1/2))/epsilon)/(exp(-1/epsilon^(1/2))+1)-2*pi^2*sin(pi*x(i))^2+2*pi^2*cos(pi*x(i))^2);
%     end
% end
ue_xx=fun_ue_xx(Xmat,Tmat,epsilon);

A=1/tau*I+0.5.*L;
B=1/tau*I-0.5.*L;
A(1,:)=I(1,:);
A(:,1)=I(:,1);
A(n,:)=I(n,:);
A(:,n)=I(:,n);

cond(A);

% inva=inv(A);
%Numerical calculation
uc=zeros(n,m);
for j=1:m-1
    C=B*uc(:,j)+0.5*(F(:,j)+F(:,j+1));
    C(1)=0;
    C(n)=0;
    uc(:,j+1)=A\C;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ޸ %%%%%%%%%%%%%%%%%%%%
% o1=sum(((uc(2:end-1,end)-ue(2:end-1,end))).^2.*(x(2:end-1)-x(1:end-2)));

o1=max(abs(uc(2:end-1,end)-ue(2:end-1,end)));

% o2=0;
% for i=2:n-1
% o2=((F(i,end)-(uc(i,end)-uc(i,end-1))/tau+epsilon*ue_xx(i,end)-q(i)'.*ue(i,end)).^2).*(x(i)-x(i-1)) +o2;
% end
% o=o1+0.5*o2;
% qexact=zeros(1,n);
qexact = qun(x)';
% for i=1:length(x)
%     qexact(i)=qun(x(i));
% end
o2=sqrt(sum((q-qexact).^2));
% o2=max(abs(q-qexact));
o=o1+o2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ޸ %%%%%%%%%%%%%%%%%%%%

% % o=abs(sum(abs(uc(:,end)-ue(:,end)))+(rand*2-1)/1000);
% o1=sum(((uc(2:end-1,end)-ue(2:end-1,end))).^2.*(x(2:end-1)-x(1:end-2)));
% % o1=max(abs(uc(2:end-1,end)-ue(2:end-1,end)));
% 
% o2=0;
% for i=2:n-1
% o2=((F(i,end)-(uc(i,end)-uc(i,end-1))/tau+epsilon*ue_xx(i,end)-q(i)'.*ue(i,end)).^2).*(x(i)-x(i-1)) +o2;
% end
% 
% o=o1+0.5*o2;
end

