%写一个纵横交叉的改进方法

function [X, counter]=cross(X_mat,fobj, counter)
h=size(X_mat,1);%行
z=size(X_mat,2);%列
for i=1:h
   M1=X_mat(i,:);
   counter = counter + 1;
fit(i)=fobj(M1);
end
for i=1:2:h  %行
    r1=rand;
    r2=rand;   %0~1之间均匀分布的随机数
    c1=-1+2*rand;%-1,1之间均匀分布的随机数
    c2=-1+2*rand;
    
    %先按（1，2）（3，4）这样搭档
    for j=1:z  %列
        MJ1=X_mat(i,j);
        MJ2=X_mat(i+1,j);
         X1(i,j)=r1*MJ1+(1-r1)*MJ2+c1*(MJ1-MJ2);
         X1(i+1,j)=r2*MJ2+(1-r2)*MJ1+c2*(MJ2-MJ1);  %横的交叉完了 
    end
XI1=  X1(i,:);
XI2=  X1(i+1,:);
%对比适应度
counter = counter + 2;
fitx1=fobj(XI1);
fitx2=fobj(XI2);

if  fitx1<fit(i)
    X_mat(i,:)=XI1;
    fit(i)=fitx1;
end
if  fitx2<fit(i+1)
    X_mat(i+1,:)=XI2;
    fit(i+1)=fitx2;
end

%end   %end heng


%开始纵交叉
%for p=1:h
q1=randi([1,z]);
q2=setdiff(1:z,q1);
q3=q2(randperm(length(q2)));
r=rand;
MJ3=X_mat(i,q1);
MJ4=X_mat(i+1,q1);
X1(i,q1)=r*MJ3+(1-r)*X_mat(i,q3(1));
X1(i+1,q1)=r*MJ4+(1-r)*X_mat(i+1,q3(1));

XI3=  X1(i,:);
XI4=  X1(i+1,:);
counter = counter + 2;
fitx3=fobj(XI3);
fitx4=fobj(XI4);
%对比适应度
if  fitx3<fit(i)
    X_mat(i,:)=XI3;
    fit(i)=fitx3;
end
if  fitx3<fit(i+1)
    X_mat(i+1,:)=XI4;
    fit(i+1)=fitx3;
end
end%zong


X=X_mat;

end  %end function 