function [gbest_f,Best,gbest]=DSSA(N,iter,lb,ub,dim,fobj)
%% 参数
Gc=1.9; %飞行系数
Cr=0.5; %交叉概率
dg=0.8; %飞翔距离
Pdp=0.1;%猎人出现概率
Best=zeros(1,iter);
f=ones(1,N)*inf;
counter = 0;
%% 初始化
FS=initialization(N,dim,ub,lb);
for i=1:N
    counter = counter + 1;
    f(i)=fobj(FS(i,:));%计算适应值
end
[sort_f,sort_index]=sort(f);
FSh=FS(sort_index(1),:);    %核桃树上的飞鼠
FSa=FS(sort_index(2:4),:);  %橡树上的飞鼠
FSn=FS(sort_index(5:end),:);%普通树上的飞鼠
n1=round(rand*(N-4));
n2=N-n1-4;
L1=randperm(size(FSn,1));%随机的整数排列
FSn1=FSn(L1(1:end-n2),:);
FSn2=FSn(L1(end-n2+1:end),:);
G=sort_f(1);
Best(1)=sort_f(1);
Dic(1) = sum(sum(abs(FS - median(FS))) / N);
%% 主循环
for it=2:iter
    Div(it) = sum(sum(abs(FS - median(FS))) / N);
    FS_avg=mean(FS,1); %种群平均位置
    n1=round(rand*(N-4));
    n2=N-n1-4;
    L1=randperm(size(FSn,1));
    FSn1=FSn(L1(1:end-n2),:);
    FSn2=FSn(L1(end-n2+1:end),:);
    %临时树
    FSa_t=[];
    FSn1_t=[];
    FSn2_t=[]; 
    %% 橡树->核桃树
    for i=1:3
        if rand>Pdp
            FSa_t(i,:)=FSa(i,:)+dg*Gc*(FSh-FSa(i,:)-FS_avg);%Eq(1)              
        else
            Xrand=initialization(1,dim,ub,lb);
            FSa_t(i,:)=Xrand-rand(1,dim).*abs(Xrand-2*rand*FSa(i,:));
        end
        FSa_t(i,:)=FSa_t(i,:)+(lb-FSa_t(i,:)).*(FSa_t(i,:)<lb)+(ub-FSa_t(i,:)).*(FSa_t(i,:)>ub);%边界条件
        %交叉
        Xa_t=[];
        Jrand=randi([1,dim],1);        
        for k=1:dim
            if ((k==Jrand)||(rand<=Cr))%Eq(2)
                Xa_t(1,k)=FSa_t(i,k);
            else
                Xa_t(1,k)=FSa(i,k);
            end
        end
        counter = counter + 1;
        if fobj(Xa_t)<fobj(FSa_t(i,:))
            FSa_t(i,:)=Xa_t;     
        end    
    end
    %% 普通树->橡树
    for i=1:n1
        if rand>Pdp
            k=randi(3,1);
            FSn1_t(i,:)=FSn1(i,:)+dg*Gc*(FSa(k,:)-FSn1(i,:));%Eq(3)
        else
            Xrand=initialization(1,dim,ub,lb);
            FSn1_t(i,:)=Xrand-rand(1,dim).*abs(Xrand-2*rand*FSn1(i,:));
        end
        FSn1_t(i,:)=FSn1_t(i,:)+(lb-FSn1_t(i,:)).*(FSn1_t(i,:)<lb)+(ub-FSn1_t(i,:)).*(FSn1_t(i,:)>ub);%边界条件
        %交叉
        Yn_t=[];
        Jrand=randi([1,dim],1);        
        for k=1:dim
            if ((k==Jrand)||(rand<=Cr))%Eq(5)
                Yn_t(1,k)=FSn1_t(i,k);
            else
                Yn_t(1,k)=FSn1(i,k);
            end
        end          
        counter = counter + 1;
        if fobj(Yn_t)<fobj(FSn1_t(i,:))
              FSn1_t(i,:)=Yn_t;
        end  
    end
    %% 普通树->核桃树
    for i=1:n2
        if rand>Pdp
            FSn2_t(i,:)=FSn2(i,:)+dg*Gc*(FSh-FSn2(i,:));%Eq(4)
        else
            Xrand=initialization(1,dim,ub,lb);
            FSn2_t(i,:)=Xrand-rand*abs(Xrand-2*rand*FSn2(i,:)); 
        end 
        FSn2_t(i,:)=FSn2_t(i,:)+(lb-FSn2_t(i,:)).*(FSn2_t(i,:)<lb)+(ub-FSn2_t(i,:)).*(FSn2_t(i,:)>ub);%边界条件
        %交叉
        Zn_t=[];
        Jrand=randi([1,dim],1);        
        for k=1:dim
            if ((k==Jrand)||(rand<=Cr))%Eq(5)
                Zn_t(1,k)=FSn2_t(i,k);
            else
                Zn_t(1,k)=FSn2(i,k);
            end
        end
        counter = counter + 1;
        if fobj(Zn_t)<fobj(FSn2_t(i,:))
             FSn2_t(i,:)=Zn_t;
        end          
    end
    counter = counter + 1;
    FSn_t=[FSn1_t;FSn2_t];
    FSh_new=FSh+dg*Gc*(FSh-mean(FSa,1));%Eq(6)
    if fobj(FSh_new)<fobj(FSh)
        FSh=FSh_new;
    end
    FS_new=[FSh;FSa_t;FSn_t];  
    %% 更新种群
    FT_new=[];
    for i=1:N    
        counter = counter + 1;
        FT_new(i)=fobj(FS_new(i,:));
        if FT_new(i)<f(i)
            f(i)=FT_new(i);
            FS(i,:)=FS_new(i,:);
        end
    end
    %% 更新最优
    [sort_f,sort_index]=sort(f);
    if sort_f(1)<G(it-1)
        FSh=FS(sort_index(1),:); 
        G(it)=sort_f(1);
    else
        G(it)=G(it-1);
    end
    FSa=FS(sort_index(2:4),:);
    FSn=FS(sort_index(5:end),:);
    %% 保存循环
    Best(it)=G(it); 
%     display(['DSSA:The optimum at iteration ' , num2str(it), ' is ', num2str(Best(it))]);
end
%%
gbest_f=Best(end);
gbest=FSh;
end