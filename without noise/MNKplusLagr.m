clc
clear all
%-----------------------------------------
%Функція
%-----------------------------------------


for f=[2,20,200]%Частота періодичного сигналу
    T=1/f;%Період сигналу
    p=6;%Кількість періодів у сигналі
    N=512;%Кількість відліків (точок) у сигналі та Кількість точок для інтерполяції
    Td=p*T/N;%Період дискретизації
    fd=1/Td;%Частота дискретизації
    x=Td:Td:p*T;%Вектор часу (Вектор аргументів початкової функції з 'N' точок)
    y=3+sin(2*pi*x*f).*exp(-7*x);%Вектор значень функції з 'N' точок
    ly=length(y);%Довжина сигналу
    lx=length(x);%Довжина вектора часу
    %-----------------------------------------
    %Визначення степеня двійки
    %-----------------------------------------
    i=0;%Інкремент
    div=ly/(2^i);
    while div>1
        div=ly/(2^(i+1));
        if div>=1
            i=i+1;
        end
    end
    %-----------------------------------------
    %Визначення остачі від ділення на 2^i
    %-----------------------------------------
    if ly-2^i~=0%Умова наявності остачі
        Info=['Довжина сигналу складає ' num2str(ly) ' відліки'];
        disp('-----------------------------------------------------')
        disp(Info)
        disp('-----------------------------------------------------')
        Info=['Найближча степінь двійки, i=' num2str(i)];
        disp(Info)
        disp('-----------------------------------------------------')
        Ost=ly-2^i;%Остача від ділення на 2^i
        Info=['Остача від ділення на 2^i=' num2str(Ost)];
        disp(Info)
        disp('-----------------------------------------------------')
        %-----------------------------------------
        %Інтерполяція
        %-----------------------------------------
        na=15;%Початкова кількість точок у відрізку
        while (rem(lx,na)~=0)
            na=na+1;%Визначення дільника, більшого за 15
            %-----------------------------------------
            %Зменшення кількості початкових відліків до ділення на 32
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
        inTo=2^(i+1)-ly;%Визначення кількості відліків, необхідних для інтерполяції
        Na=fix(ly/na);%Визначення кількості відрізків інтерполяції
        vpls=fix(inTo/Na);%Кількість додаткових точок
        add=rem(inTo,Na);%Неповна кількість відліків, для додавання у відрізки
        Lint=vpls+na+1;
        for j=1:Na
            nl=(j-1)*na+1;
            kl=na*j;
            nli=(j-1)*Lint+1;
            kli=Lint*j;
            [xi(nli:kli),L(nli:kli)]=Lagr_rep(x(nl:kl),y(nl:kl),Lint);
            %Видалення зайвих елементів (кінець алгоритму на 115 рядку)
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
        %Апроксимація
        %-----------------------------------------
        lL=length(L);
        lxi=length(xi);
        if lL<=32
            Na=4;%Кількість проміжків апроксимації
        else
            Na=8;
        end
        Pp=5;%Максимальний порядок полінома
        if lL~=lxi
            disp('Довжина сигналу не рівна довжині вектора часу!')
        else
            na=lxi/Na;%Кількість проміжків апроксимації
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
                        disp('Ділення на нуль! Змініть вхідні дані!')
                    else
                        if Do<1e-7
                            Do=0;
                        end
                        R_2=1-Do/Ds;%Визначення коефіцієнту детермінації
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
            %Згладжування
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
        Ost=0;%Остача від ділення на 2^i
        Info=['Довжина сигналу складає ' num2str(ly) ' відліки'];
        disp('-----------------------------------------------------')
        disp(Info)
        disp('-----------------------------------------------------')
        disp('Остача від ділення на 2^i відсутня')
        disp('-----------------------------------------------------')
        Info=['Cтепінь двійки, i=' num2str(i)];
        disp(Info)
        disp('-----------------------------------------------------')
        %-----------------------------------------
        %Апроксимація
        %-----------------------------------------
        if ly<=32
            Na=4;%Кількість проміжків апроксимації
        else
            Na=8;
        end
        Pp=7;%Максимальний порядок полінома
        if ly~=lx
            disp('Довжина сигналу не рівна довжині вектора часу!')
        else
            na=lx/Na;%Кількість проміжків апроксимації
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
                        disp('Ділення на нуль! Змініть вхідні дані!')
                    else
                        if Do<1e-7
                            Do=0;
                        end
                        R_2=1-Do/Ds;%Визначення коефіцієнту детермінації
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
            %Згладжування
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
            %Побудова графіків
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
            title('Графік сигналу')
            legend('Початкова функція','Апроксимована функція')
            xlabel('Час, с')
            ylabel('Амплітуда')
            subplot(3,1,2)
            plot(xi,PohAbs),grid
    %         title('Graph of absolute error')
    %         xlabel('Time, s')
    %         ylabel('Magnitude')
            title('Графік абсолютної похибки')
            xlabel('Час, с')
            ylabel('Амплітуда')
            subplot(3,1,3)
            plot(xi,PohOtn),grid
    %         title('Graph of relative error')
    %         xlabel('Time, s')
    %         ylabel('Relative error, %')
            title('Графік відносної похибки')
            xlabel('Час, с')
            ylabel('Відносна похибка, %')
end