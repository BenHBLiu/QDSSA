function [Best_score,Best_pos,Convergence_curve]=SAO(N,Max_iter,lb,ub,dim,fobj)
if(max(size(ub)) == 1)
   ub = ub.*ones(1,dim);
   lb = lb.*ones(1,dim);  
end
%% Initialize the set of random solutions
X=initialization(N,dim,ub,lb);%种群
% X=repmat(lb, N, 1)+chaos(3,N,dim) .* repmat((ub-lb), N, 1);
% a=rand(1,dim);
% X(1,:)=lb+(ub-lb).*a;
% for i=1:N-1
%     a=4*a.^3-3*a;
%     X(i+1,:)=lb+(ub-lb).*a;
% end
counter = 0;
Objective_values = zeros(1,size(X,1));%种群适应值
Convergence_curve=[];%优化曲线
N1=floor(N*0.5);
Elite_pool=[];
% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    counter = counter + 1;
    Objective_values(i,1)=fobj(X(i,:));
end
[~,idx1]=sort(Objective_values);%对适应值从小到大排序，[n,m]=sort()，n装着排序结果，m装着对应位置
Best_pos=X(idx1(1),:);
Best_score=Objective_values(idx1(1),1);
second_best=X(idx1(2),:);
third_best=X(idx1(3),:);
sum1=0;
for i=1:N1
    sum1=sum1+X(idx1(i),:);
end
half_best_mean=sum1/N1;
gbest = Best_pos;
gbest_f = Best_score;
Elite_pool(1,:)=Best_pos;
Elite_pool(2,:)=second_best;
Elite_pool(3,:)=third_best;
Elite_pool(4,:)=half_best_mean;
Convergence_curve(1) = Best_score;
for i=1:N
    index(i)=i;
end
Na=N/2;
Nb=N/2;
Div = zeros(1, Max_iter);
Exploration = zeros(1, Max_iter);
Exploitation = zeros(1, Max_iter);
%% Main loop
for it=1:Max_iter
    Div(it) = sum(sum(abs(X - median(X))) / N);
    RB=randn(N,dim);          %Brownian random number vector
    T=exp(-it/Max_iter);
    DDF=0.35+0.25*(exp(it/Max_iter)-1)/(exp(1)-1);%eq.(9)
    M=DDF*T;  %eq.(10)
    %% Calculate the centroid position of the entire population
    for j=1:dim
        sum1=0;
        for i=1:N
            sum1=sum1+X(i,j);
        end
        X_centroid(1,j)=sum1/N;
    end   
    %% Select individuals randomly to construct pop1 and pop2
    index1=randperm(N,Na);%在[1,N]中生成Na个不重复的随机整数列
    index2=setdiff(index,index1);  %求index和index1的差集，如：index=[1,2,3,4,5],index1=[1,2,3,4,5,6,7,8,9,10]则setdiff(index,index1)=[6,7,8,9,10]
    for i=1:Na
        r1=rand(1,dim);
        k1=randperm(4,1);
        X(index1(i),:)= Elite_pool(k1,:)+RB(index1(i),:).*(r1.*(Best_pos-X(index1(i),:))+(1-r1).*(X_centroid-X(index1(i),:)));%eq.(3)
    end
    if Na<N
        Na=Na+1;
        Nb=Nb-1;
    end
    if Nb>=1
        for i=1:Nb
            r2=2*rand-1;
            X(index2(i),:)= M*Best_pos+RB(index2(i),:).*(r2*(Best_pos-X(index2(i),:))+(1-r2)*(X_centroid-X(index2(i),:)));%eq.(11)
        end
    end
    % Check if solutions go outside the search spaceand bring them back
    for i=1:size(X,1)
        X(i,:)=X(i,:)+(lb-X(i,:)).*(X(i,:)<lb)+(ub-X(i,:)).*(X(i,:)>ub);%边界条件
        % Calculate the objective values
        Objective_values(i,1)=fobj(X(i,:));
        counter = counter + 1;
        % Update the destination if there is a better solution
        if Objective_values(i,1)<Best_score
            Best_pos=X(i,:);
            Best_score=Objective_values(i,1);
        end
    end
    %% Update the elite pool
    [~,idx1]=sort(Objective_values);
    second_best=X(idx1(2),:);
    third_best=X(idx1(3),:);
    sum1=0;
    for i=1:N1
        sum1=sum1+X(idx1(i),:);
    end
    half_best_mean=sum1/N1;
    Elite_pool(1,:)=Best_pos;
    Elite_pool(2,:)=second_best;
    Elite_pool(3,:)=third_best;
    Elite_pool(4,:)=half_best_mean;
    Convergence_curve(it)=Best_score;
%     display(['The optimum at iteration ' , num2str(it), ' is ', num2str(Convergence_curve(it))]);
%     it=it+1;
end
% figure()
% plot(Exploration, 'r')
% hold on
% plot(Exploitation, 'b')
% ylim([0 100]);
% legend(['Exploration (Avg ',num2str(mean(Exploration)), ')'], ['Exploitation (Avg ',num2str(mean(Exploitation)), ')'])
end