clc
clear all
%-----------------------------------------
%�������
%-----------------------------------------


for f=[2,20,200]%������� ����������� �������
    T=1/f;%����� �������
    p=6;%ʳ������ ������ � ������
    N=512;%ʳ������ ����� (�����) � ������ �� ʳ������ ����� ��� ������������
    Td=p*T/N;%����� �������������
    fd=1/Td;%������� �������������
    x=Td:Td:p*T;%������ ���� (������ ��������� ��������� ������� � 'N' �����)
    y=3+sin(2*pi*x*f).*exp(-7*x);%������ ������� ������� � 'N' �����
    ly=length(y);%������� �������
    lx=length(x);%������� ������� ����
    %-----------------------------------------
    %���������� ������� �����
    %-----------------------------------------
    i=0;%���������
    div=ly/(2^i);
    while div>1
        div=ly/(2^(i+1));
        if div>=1
            i=i+1;
        end
    end
    %-----------------------------------------
    %���������� ������ �� ������ �� 2^i
    %-----------------------------------------
    if ly-2^i~=0%����� �������� ������
        Info=['������� ������� ������ ' num2str(ly) ' �����'];
        disp('-----------------------------------------------------')
        disp(Info)
        disp('-----------------------------------------------------')
        Info=['��������� ������ �����, i=' num2str(i)];
        disp(Info)
        disp('-----------------------------------------------------')
        Ost=ly-2^i;%������ �� ������ �� 2^i
        Info=['������ �� ������ �� 2^i=' num2str(Ost)];
        disp(Info)
        disp('-----------------------------------------------------')
        %-----------------------------------------
        %������������
        %-----------------------------------------
        na=15;%��������� ������� ����� � ������
        while (rem(lx,na)~=0)
            na=na+1;%���������� �������, ������� �� 15
            %-----------------------------------------
            %��������� ������� ���������� ����� �� ������ �� 32
            %-----------------------------------------
            if na>40
                ji=32*floor(lx/32);
                jj=floor(lx/(lx-ji))-1;
                j=jj;
                while(lx~=ji)
                    x(j)=[];
                    y(j)=[];
                    j=j+jj;
                    lx=length(x);
                    ly=length(y);
                end
                na=30;
            end
            %-----------------------------------------
        end
        inTo=2^(i+1)-ly;%���������� ������� �����, ���������� ��� ������������
        Na=fix(ly/na);%���������� ������� ������ ������������
        vpls=fix(inTo/Na);%ʳ������ ���������� �����
        add=rem(inTo,Na);%������� ������� �����, ��� ��������� � ������
        Lint=vpls+na+1;
        for j=1:Na
            nl=(j-1)*na+1;
            kl=na*j;
            nli=(j-1)*Lint+1;
            kli=Lint*j;
            [xi(nli:kli),L(nli:kli)]=Lagr_rep(x(nl:kl),y(nl:kl),Lint);
            %��������� ������ �������� (����� ��������� �� 115 �����)
            index(nli:kli)=0;
            for ji=nl:kl
                for jj=nli:kli
                    if xi(jj)==x(ji)
                        index(jj)=1;
                    end
                end
            end
            if add~=0
                if j==add
                    xint=xi;
                end
            else
                xint=1;
            end
        end
        ji=Na-add;
        j=1;
        while (ji~=0)
            for jj=(length(xint)+(j-1)*Lint+round(Lint/10)):1:length(xi)%(length(xi)-round(Lint/10)):-j:length(xint)
                if index(jj)==0
                    xi(jj)=66666666666666666.6666666666666;
                    ji=ji-1;
                    j=j+1;
                    break;
                end
            end
        end
        j = 1;
        while j <= length(xi)
            if xi(j)==66666666666666666.6666666666666
                xi(j)=[];
                L(j)=[];
            else
                j=j+1;
            end
        end
        %-----------------------------------------
        %������������
        %-----------------------------------------
        lL=length(L);
        lxi=length(xi);
        if lL<=32
            Na=4;%ʳ������ ������� ������������
        else
            Na=8;
        end
        Pp=5;%������������ ������� �������
        if lL~=lxi
            disp('������� ������� �� ���� ������ ������� ����!')
        else
            na=lxi/Na;%ʳ������ ������� ������������
            for j=1:Na
                R_2_max=0;
                por=0;
                nl=(j-1)*na+1;
                kl=na*j;
                for por_t=1:Pp
                    p=polyfit(xi(nl:kl),L(nl:kl),por_t);
                    a(nl:kl)=polyval(p,xi(nl:kl));
                    Do=sum((L(nl:kl)-a(nl:kl)).^2);
                    yn1=length(L(nl:kl));
                    ysr=sum(L(nl:kl))/yn1;
                    Ds=sum((L(nl:kl)-ysr).^2);
                    if Ds==0
                        disp('ĳ����� �� ����! ����� ����� ���!')
                    else
                        if Do<1e-7
                            Do=0;
                        end
                        R_2=1-Do/Ds;%���������� ����������� �����������
                        if R_2_max<R_2
                            R_2_max=R_2;
                            por=por_t;
                        end
                    end
                end
                p=polyfit(xi(nl:kl),L(nl:kl),por);
                a(nl:kl)=polyval(p,xi(nl:kl));
            end
            %-----------------------------------------
            %������������
            %-----------------------------------------
            lot=lxi/Na;
            olot=round(lot/5);
            for j=1:(Na-1)
                nl=j*lot-olot+1;
                kl=j*lot+olot;
                p=polyfit(xi(nl:kl),a(nl:kl),5);
                a(nl:kl)=polyval(p,xi(nl:kl));
            end
        end
    else
        Ost=0;%������ �� ������ �� 2^i
        Info=['������� ������� ������ ' num2str(ly) ' �����'];
        disp('-----------------------------------------------------')
        disp(Info)
        disp('-----------------------------------------------------')
        disp('������ �� ������ �� 2^i �������')
        disp('-----------------------------------------------------')
        Info=['C����� �����, i=' num2str(i)];
        disp(Info)
        disp('-----------------------------------------------------')
        %-----------------------------------------
        %������������
        %-----------------------------------------
        if ly<=32
            Na=4;%ʳ������ ������� ������������
        else
            Na=8;
        end
        Pp=7;%������������ ������� �������
        if ly~=lx
            disp('������� ������� �� ���� ������ ������� ����!')
        else
            na=lx/Na;%ʳ������ ������� ������������
            for j=1:Na
                R_2_max=0;
                por=0;
                nl=(j-1)*na+1;
                kl=na*j;
                for por_t=1:Pp
                    p=polyfit(x(nl:kl),y(nl:kl),por_t);
                    a(nl:kl)=polyval(p,x(nl:kl));
                    Do=sum((y(nl:kl)-a(nl:kl)).^2);
                    yn1=length(y(nl:kl));
                    ysr=sum(y(nl:kl))/yn1;
                    Ds=sum((y(nl:kl)-ysr).^2);
                    if Ds==0
                        disp('ĳ����� �� ����! ����� ����� ���!')
                    else
                        if Do<1e-7
                            Do=0;
                        end
                        R_2=1-Do/Ds;%���������� ����������� �����������
                        if R_2_max<R_2
                            R_2_max=R_2;
                            por=por_t;
                        end
                    end
                end
                p=polyfit(x(nl:kl),y(nl:kl),por);
                a(nl:kl)=polyval(p,x(nl:kl));
            end
            %-----------------------------------------
            %������������
            %-----------------------------------------
            lot=length(x)/Na;
            olot=round(lot/5);
            for j=1:(Na-1)
                nl=j*lot-olot+1;
                kl=j*lot+olot;
                p=polyfit(x(nl:kl),a(nl:kl),5);
                a(nl:kl)=polyval(p,x(nl:kl));
            end         
    L=y;
    xi=x;
        end
    end
            %-----------------------------------------
            %�������� �������
            %-----------------------------------------
            PohAbs=abs(L-a);
            PohOtn=PohAbs./L*100;
            figure()
            subplot(3,1,1)
            plot(x,y,xi,a,'--'),grid
    %         title('Signal graph')
    %         legend('Initial function','Approximated function')
    %         xlabel('Time, s')
    %         ylabel('Magnitude')
            title('������ �������')
            legend('��������� �������','������������� �������')
            xlabel('���, �')
            ylabel('��������')
            subplot(3,1,2)
            plot(xi,PohAbs),grid
    %         title('Graph of absolute error')
    %         xlabel('Time, s')
    %         ylabel('Magnitude')
            title('������ ��������� �������')
            xlabel('���, �')
            ylabel('��������')
            subplot(3,1,3)
            plot(xi,PohOtn),grid
    %         title('Graph of relative error')
    %         xlabel('Time, s')
    %         ylabel('Relative error, %')
            title('������ ������� �������')
            xlabel('���, �')
            ylabel('³������ �������, %')
end