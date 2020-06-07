function [a,L,xi]=MNK_Lagr(x,y,N)
%x - ���������� ������ ���������
%y - ���������� ������ ������� �������
%N - ʳ������ ����� � ����������� ������
%x� - ��������������� ������ ���������
%a - ������ ������� ������������� �������
%L - ������ ������� �������������� �������
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
a=MNK_Part(xi, L);
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
a=MNK_Part(x, y);    
L=y;
xi=x;
end