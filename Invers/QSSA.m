function [gbest_f,Best,gbest]=QSSA(N,MaxIT,lb,ub,dim,fobj)
Ns=3;     %����
sf=32;
Gc=1.9;
Pdp=0.2;  %�����ߴ��ڸ���
n1=Ns;
sc_count=0;
ps=0.5;
omega1=0.6;
omega2=0.9;
c1=0.6;
c2=0.9;
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
Dic(1) = sum(sum(abs(FS - median(FS))) / N);
%% ��ѭ��
for it=2:MaxIT
    Div(it) = sum(sum(abs(FS - median(FS))) / N);
    omega=omega2+(omega1-omega2)*it/MaxIT;
    %% ����
    for t=1:n1
        if rand>Pdp
            dg=(11*rand+9)/sf;
            FSa(t,:)=FSa(t,:)+dg*Gc*(FSh-FSa(t,:));
        else
            u=rand;
            FSa(t,:)=FSh+omega*abs(FSh--FSa(t,:)).*log(1./u)*(1-2*(rand>=0.5));
        end
        FSa(t,:)=FSa(t,:)+(lb-FSa(t,:)).*(FSa(t,:)<lb)+(ub-FSa(t,:)).*(FSa(t,:)>ub);
    end
    %% ��ͨ��
    for t=1:N-n1-1
        r2=rand;
        r3=rand;
        r4=rand;
        if r2>=Pdp && r4>=ps
            dg=(11*rand+9)/sf;
            t1=randi(3);
            FSn(t,:)=FSn(t,:)+dg*Gc*(FSa(t1,:)-FSn(t,:));
        elseif r4<ps && r3>=Pdp
            dg=(11*rand+9)/sf;
            FSn(t,:)=FSn(t,:)+dg*Gc*(FSh-FSn(t,:));
        else
            rr1=rand;
            rr2=rand;
            phi=(c1*rr1)/(c1*rr1+c2*rr2);
            pnt=phi*sum(FSa)/n1+(1-phi)*FSh;
            u=rand;
            FSn(t,:)=pnt+omega*abs(pnt-FSn(t,:)).*log(1./u)*(1-2*(rand>=0.5));
        end
        FSn(t,:)=FSn(t,:)+(lb-FSn(t,:)).*(FSn(t,:)<lb)+(ub-FSn(t,:)).*(FSn(t,:)>ub);
    end
    %% �жϼ���
    Sc=0;
    for z=1:Ns
        for k=1:dim
            Sc=Sc+(FSa(z,k)-FSh(1,k))^2;%���㼾�ڳ���
        end
    end
    Sc=sqrt(Sc);
    Smin=(10*exp(-6))/365^(it/(MaxIT/2.5));
    if Sc<Smin
%         sc_count=sc_count+1;
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
%     display(['QSSA:The optimum at iteration ' , num2str(it), ' is ', num2str(Best(it))]);
end
gbest_f=Best(end);
end