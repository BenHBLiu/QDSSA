function [best_f, trace, gbest] = GA(NP, G, Xx, Xs, D, fobj)
Xx = Xx * ones(1, D);
Xs = Xs * ones(1, D);
nf=zeros(NP,D);                       	%子种群赋空间
Pc=0.8;                              	%交叉概率
Pm=0.1;                             	%变异概率
f=rand(NP,D).*(Xs-Xx)+Xx;             	%随机获得初始种群
counter = 0;
%%%%%%%%%%%%%%按适应度升序排列%%%%%%%%%%%%%%%
for np=1:NP
    MSLL(np)=fobj(f(np,:));
    counter = counter + 1;
end
[SortMSLL,Index]=sort(MSLL);                            
Sortf=f(Index,:);
%%%%%%%%%%%%%%遗传算法循环%%%%%%%%%%%%%%%%%
for gen=1:G
    Div(gen) = sum(sum(abs(f - median(f))) / NP);
    %%%%%%%%%采用君主方案进行选择交叉操作%%%%%%%%%%
    Emper=Sortf(1,:);                      	%君主染色体
    NoPoint=round(D*Pc);                   	%每次交叉点的个数
    PoPoint=randi([1 D],NP/2,NoPoint);   	%交叉基因的位置
    nf=Sortf;
    for i=1:NP/2
        nf(2*i-1,:)=Emper;
        nf(2*i,:)=Sortf(2*i,:);
        for k=1:NoPoint
            nf(2*i-1,PoPoint(i,k))=nf(2*i,PoPoint(i,k));
            nf(2*i,PoPoint(i,k))=Emper(PoPoint(i,k));
        end
    end
    %%%%%%%%%%%%变异操作%%%%%%%%%%%%%%%%%%%
    for m=1:NP
        for n=1:D
            r=rand(1,1);
            if r<Pm
                nf(m,n)=rand(1,1)*(Xs(n)-Xx(n))+Xx(n);
            end
        end
    end
    %%%%%%%%%子种群按适应度升序排列%%%%%%%%%%%%%%
    for np=1:NP 
        counter = counter + 1;
          NMSLL(np)=fobj(nf(np,:));   
    end
    [NSortMSLL,Index]=sort(NMSLL);           
    NSortf=nf(Index,:);
    %%%%%%%%%产生新种群%%%%%%%%%%%%%%%%%%%%%
    f1=[Sortf;NSortf];                	%子代和父代合并
    MSLL1=[SortMSLL,NSortMSLL];       	%子代和父代的适应度值合并
    [SortMSLL1,Index]=sort(MSLL1);    	%适应度按升序排列
    Sortf1=f1(Index,:);               	%按适应度排列个体
    SortMSLL=SortMSLL1(1:NP);         	%取前NP个适应度值
    Sortf=Sortf1(1:NP,:);             	%取前NP个个体
    trace(gen)=SortMSLL(1);           	%历代最优适应度值
end
gbest = Sortf(1,:);                     	%最优个体 
best_f = trace(end);                            	%最优值
end