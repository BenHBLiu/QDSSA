function [gbesst_f,gbest,Best]=SSA(N,MaxIT,lb,ub,dim,fobj)
Ns=3;     %����
Gc=1.9;   %����ϵ��
Pdp=0.1;  %�����ߴ��ڸ���
n1=Ns;
sc_count=0;
pdp=0;
%%  ��ʼ��
FS=initialization(N,dim,ub,lb);%����(flying squirrels)Ⱥ
for i=1:N
    sc_count=sc_count+1;
    f(i)=fobj(FS(i,:));%������Ӧֵ
end
[F,index]=sort(f);%���� 
FSh=FS(index(1),:);     %�������ϵķ���
FSa=FS(index(2:4),:);   %�����ϵķ���
FSn=FS(index(5:end),:); %��ͨ���ϵķ���
G=F(1);
Best(1)=F(1);
%% ��ѭ��
for it=2:MaxIT
    Div(it) = sum(sum(abs(FS - median(FS))) / N);
    n2=round(rand*(N-Ns-1));%��ͨ���ɵ��������ķ�������
    n3=N-Ns-1-n2;           %��ͨ���ɵ������ķ�������
    L1=randperm(size(FSn,1));
    FSn2=FSn(L1(1:end-n3),:);
    FSn3=FSn(L1(end-n3+1:end),:);
    %% ����->������
    for t=1:n1
        r1=rand;
        if r1>Pdp
            dg=0.5+(1.11-0.5)*rand;
            FSa(t,:)=FSa(t,:)+dg*Gc*(FSh-FSa(t,:));
            FSa(t,:)=FSa(t,:)+(lb-FSa(t,:)).*(FSa(t,:)<lb)+(ub-FSa(t,:)).*(FSa(t,:)>ub);
        else
            pdp=pdp+1;
            FSa(t,:)=lb+rand(1,dim).*(ub-lb);
        end
    end
    %% ��ͨ��->������
    for t=1:n2
        r2=rand;
        if r2>Pdp
            dg=0.5+(1.11-0.5)*rand;
            t1=randi(3);%���ѡ��һ�ú�����
            FSn2(t,:)=FSn2(t,:)+dg*Gc*(FSa(t1,:)-FSn2(t,:));
            FSn2(t,:)=FSn2(t,:)+(lb-FSn2(t,:)).*(FSn2(t,:)<lb)+(ub-FSn2(t,:)).*(FSn2(t,:)>ub);
        else
            pdp=pdp+1;
            FSn2(t,:)=lb+rand(1,dim).*(ub-lb);
        end
    end
    %% ��ͨ��->����
    for t=1:n3
        r3=rand;
        if r3>Pdp
            dg=0.5+(1.11-0.5)*rand;
            FSn3(t,:)=FSn3(t,:)+dg*Gc*(FSh-FSn3(t,:));
            FSn3(t,:)=FSn3(t,:)+(lb-FSn3(t,:)).*(FSn3(t,:)<lb)+(ub-FSn3(t,:)).*(FSn3(t,:)>ub);
        else
            pdp=pdp+1;
            FSn3(t,:)=lb+rand(1,dim).*(ub-lb);
        end
    end
    %% �жϼ���
    FSn=[FSn2;FSn3];
    Sc=0;
    for z=1:Ns
        for k=1:dim
            Sc=Sc+(FSa(z,k)-FSh(1,k))^2;%���㼾�ڳ���
        end
    end
    Sc=sqrt(Sc);
    Smin=10*exp(-6)/365^(it/(MaxIT/2.5));
    if Sc<Smin
        
        for i=1:N-Ns-1
            FSn(i,:)=FSn(i,:)+levy(1,dim,1.5).*(ub-lb);
            FSn(i,:)=FSn(i,:)+(lb-FSn(i,:)).*(FSn(i,:)<lb)+(ub-FSn(i,:)).*(FSn(i,:)>ub);
        end
    end
    %% ����
    FS=[FSh;FSa;FSn];
    for i=1:N
        sc_count=sc_count+1;
        f(i)=fobj(FS(i,:));
    end
    [F,index]=sort(f);
    if F(1)<G(it-1)
        FSh=FS(index(1),:);
        G(it)=F(1);
    else
        G(it)=G(it-1);
    end
    FSa=FS(index(2:4),:);
    FSn=FS(index(5:end),:);
    Best(it)=G(it);
    gbest=FSh;
%     display(['SSA:The optimum at iteration ' , num2str(it), ' is ', num2str(Best(it))]);
end
gbesst_f=Best(end);
end