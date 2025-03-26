function [Fbest,Best,gbest]=PSO(N,iter,lb,ub,D,fobj)
c1 = 2;
c2 = 2;
w = 0.9;
counter = 0;
%% 初始化
x=lb+rand(N,D).*(ub-lb);%位置
for i=1:N
    counter = counter + 1;
    fit(i,1)=fobj(x(i,:));%适应值 
end
v=zeros(N,D);%速度
v_min=lb*0.2;
v_max=ub*0.2;
pbest=x;%个体最优解
[best_f,best_index]=min(fit);
gbest=x(best_index,:);%全局最优解
Fbest=best_f;%全局最优
%% 主循环
for it=1:iter
    Div(it) = sum(sum(abs(x - median(x))) / N);
    w=w-1/(2*iter);%惯量权重
    for i=1:N
        v(i,:)=w*v(i,:)+c1*rand*(pbest(i,:)-x(i,:))+c2*rand*(gbest-x(i,:));%更新速度
        v(i,:)=v(i,:)+(v_min-v(i,:)).*(v(i,:)<v_min)+(v_max-v(i,:)).*(v(i,:)>v_max);%边界条件
        x(i,:)=x(i,:)+v(i,:);%更新位置
        x(i,:)=x(i,:)+(lb-x(i,:)).*(x(i,:)<lb)+(ub-x(i,:)).*(x(i,:)>ub);%边界条件
        f=fobj(x(i,:));%计算适应值
        counter = counter + 1;
        if fit(i,1)>f
            fit(i,1)=f;
            pbest(i,:)=x(i,:);%更新个体最优
        end
        if Fbest>f
            Fbest=f;
            gbest=x(i,:);%更新全局最优
        end
    end
    Best(it)=Fbest;%记录最优
    
end
end