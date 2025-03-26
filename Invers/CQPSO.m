function [f_gbest,best,gbest]=CQPSO(N,iter,lb,ub,dim,fobj)
%% 初始化

X=initialization(N,dim,ub,lb);%布谷鸟粒子
for i=1:N
    fit_z(i)=fobj(X(i,:));%适应值
end

[F,index]=sort(fit_z);%排序
f_gbest=F(1);       %全局最优值
gbest=X(index(1),:);%全局最优解

Y=X(N/2+1:end,:);%QPSO粒子
z=X;
X=X(1:N/2,:);
pbest=Y;%个体最优解


f_pbest=fit_z(N/2+1:end);%QPSO个体最优值
fit_CSO=fit_z(1:N/2);
fit=fit_z(N/2+1:end);

mpbest=zeros(N/2,dim);
NS=0;
best=zeros(1,iter);
%% 主循环
for it=1:iter
    alpha_QPSO=1-0.5*(log(it)/log(iter));%快速下降收缩-扩张系数
    alpha_CSO=iter*exp(-4*it/iter)*10^(-4);%列维飞行步长
    for j=1:N/2
        mpbest(j,:)=pbest(j,:)*(abs((f_pbest(j)-max(f_pbest))/(min(f_pbest)-max(f_pbest))));%带权重的平均最优位置
    end
    for i=1:N/2
        %% IQPSO，Update Y by IQPSO
        phi=rand;
        p=phi.*pbest(i,:)+(1-phi).*gbest;%产生随机点
        u=rand;
        r_QPSO=rand;
        Y(i,:)=p+alpha_QPSO*abs((sum(mpbest)/(N/2))-Y(i,:)).*log(1./u).*(1-2*(r_QPSO>=0.5));%更新位置
        Y(i,:)=Y(i,:)+(lb-Y(i,:)).*(Y(i,:)<lb)+(ub-Y(i,:)).*(Y(i,:)>ub);%边界条件
        fit(i)=fobj(Y(i,:));
        %% ICSA Update X by ICSA;
        levy_X=X(i,:)+alpha_CSO.*levy(1,1,1.5).*(X(i,:)-gbest);%列维飞行更新位置
        levy_X=levy_X+(lb-levy_X).*(levy_X<lb)+(ub-levy_X).*(levy_X>ub);%边界条件
        fit_levy_X=fobj(levy_X);
        if fit_CSO(i)>fit_levy_X
            X(i,:)=levy_X;
            fit_CSO(i)=fit_levy_X;
        end
        Pa=0.85+1.3*(it/iter)^3-0.5*(2*it/iter)^2;%动态发现概率
        omega=cos((pi*it)/(2*iter)+pi/2)+1;
        rand_X=omega*X(i,:)+rand.*heaviside(rand-Pa).*(X(randi(N/2),:)-X(randi(N/2),:));%巢寄生机制
        rand_X=rand_X+(lb-rand_X).*(rand_X<lb)+(ub-rand_X).*(rand_X>ub);%边界条件
        fit_rand_X=fobj(rand_X);
        if fit_CSO(i)>fit_rand_X
            X(i,:)=rand_X; 
            fit_CSO(i)=fit_rand_X;
        end
        %% combination，Cooperative mechanism; 
        d=rand(1,dim);
        Z=d.*X(i,:)+(1-d).*Y(i,:);%协同机制
        Z=Z+(lb-Z).*(Z<lb)+(ub-Z).*(Z>ub);%边界条件
        fit_Z=fobj(Z);
        %% Elite selection mechanism
        if fit_Z<fit_z(i)
            z(i,:)=Z;
            fit_z(i)=fit_Z;
        end
        if fit_CSO(i)>fit_z(i)
            X(i,:)=z(i,:);
            fit_CSO(i)=fit_z(i);
        end
        if fit(i)>fit_z(i)
            Y(i,:)=z(i,:);
            fit(i)=fit_z(i);
        end
        if f_pbest(i)>fit(i)
            f_pbest(i)=fit(i);
            pbest(i,:)=Y(i,:);
        end
        if f_gbest>f_pbest(i)
            f_gbest=f_pbest(i);
            gbest=pbest(i,:);
        end
        [xx,index]=min(fit_CSO);
        if f_gbest>xx
            f_gbest=xx;
            gbest=X(index,:);
        end
    end
    best(it)=f_gbest;
     %% Mechanism for Preventing Premature Puberty
    if it>iter/2  
        if best(it)==best(it-1)
            NS=NS+1;
            if NS==10
                [~,index_QPSO]=sort(fit);
                [~,index_CSO]=sort(fit_CSO);
                for i=floor(0.9*N/2):N/2
                    Y(index_QPSO(i),:)=initialization(1,dim,ub,lb);
                    fit(index_QPSO(i))=fobj(Y(index_QPSO(i),:));
                    if fit(index_QPSO(i))<f_pbest(index_QPSO(i))
                        f_pbest(index_QPSO(i))=fit(index_QPSO(i));
                        pbest(index_QPSO(i),:)=Y(index_QPSO(i),:);
                    end
                    if f_gbest>f_pbest(index_QPSO(i))
                        f_gbest=f_pbest(index_QPSO(i));
                        gbest=pbest(index_QPSO(i),:);
                    end
                    X(index_CSO(i),:)=initialization(1,dim,ub,lb);
                    fit_CSO(index_CSO(i))=fobj(X(index_CSO(i),:));
                    [xx,index]=min(fit_CSO);
                    if f_gbest>xx
                        f_gbest=xx;
                        gbest=X(index,:);
                    end
                end
                NS=0;
            end
        end
    end
    best(it)=f_gbest;
%     display(['C-QPSO:The optimum at iteration ' , num2str(it), ' is ', num2str(best(it))]);
end
end