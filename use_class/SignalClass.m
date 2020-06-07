classdef SignalClass  < handle
    properties
        N {mustBeNumeric};
        p {mustBeNumeric};
        f {mustBeNumeric};
        ly {mustBeNumeric};
        lx {mustBeNumeric};
        snr=0;
        fd=[];
        x=[];
        y=[];
        change_y=[]
        Y_noise=[];
        name_noise = '';
        PohOtn=[];
        PohAbs=[];
        xi=[];
        a=[];
    end
    methods
        function obj = updateValue(obj)
            T=1/obj.f;
            Td=obj.p*T/obj.N;
            obj.fd=1/Td;
            obj.x=Td:Td:obj.p*T;
            obj.y=3+sin(2*pi*obj.x*obj.f).*exp(-7*obj.x);
            obj.ly=length(obj.y);
            obj.lx=length(obj.x);
            if obj.snr == 0
                obj.change_y = obj.y;
            end
       end
       
       function [koef] = find_koef(obj, signal, noise, snr)
           koef = sqrt( mean(signal.^2) / ( 10^(snr/10)* mean(noise.^2) ) );
       end

       function obj = add_noise(obj, name_noise)
            if obj.snr ~= 0
                if name_noise == "white"
                    obj.Y_noise = function_add_white_noise(1, obj.ly);
                elseif name_noise == "pink"
                    obj.Y_noise = function_add_pink_noise(1, obj.ly);
                elseif name_noise == "violet"
                    obj.Y_noise = function_add_violet_noise(1, obj.ly);
                elseif name_noise == "grey"
                    obj.Y_noise = function_add_grey_noise(1, obj.ly);
                elseif name_noise == "broun"
                    obj.Y_noise = function_add_broun_noise(1, obj.ly);
                elseif name_noise == "blue"
                    obj.Y_noise = function_add_blue_noise(1, obj.ly);
                end 
                if ~isempty(obj.Y_noise)
                    koef = obj.find_koef(obj.y, obj.Y_noise, obj.snr);
                    obj.change_y = obj.y + obj.Y_noise*koef;
                end
            end
       end
       
       function obj =  MNKplusLagr(obj)
           if isempty(obj.Y_noise)
               obj.change_y = obj.y
           end
           i=0;
            div=obj.ly/(2^i);
            while div>1
                div=obj.ly/(2^(i+1));
                if div>=1
                    i=i+1;
                end
            end
            %-----------------------------------------
            %Визначення остачі від ділення на 2^i
            %-----------------------------------------
            if obj.ly-2^i~=0%Умова наявності остачі
                Info=['Довжина сигналу складає ' num2str(obj.ly) ' відліки'];
                disp('-----------------------------------------------------')
                disp(Info)
                disp('-----------------------------------------------------')
                Info=['Найближча степінь двійки, i=' num2str(i)];
                disp(Info)
                disp('-----------------------------------------------------')
                Ost=obj.ly-2^i;%Остача від ділення на 2^i
                Info=['Остача від ділення на 2^i=' num2str(Ost)];
                disp(Info)
                disp('-----------------------------------------------------')
                %-----------------------------------------
                %Інтерполяція
                %-----------------------------------------
                na=15;%Початкова кількість точок у відрізку
                while (rem(obj.lx,na)~=0)
                    na=na+1;%Визначення дільника, більшого за 15
                    %-----------------------------------------
                    %Зменшення кількості початкових відліків до ділення на 32
                    %-----------------------------------------
                    if na>40
                        ji=32*floor(obj.lx/32);
                        jj=floor(obj.lx/(obj.lx-ji))-1;
                        j=jj;
                        while(obj.lx~=ji)
                            obj.x(j)=[];
                            obj.change_y(j)=[];
                            j=j+jj;
                            obj.lx=length(obj.x);
                            obj.ly=length(obj.change_y);
                        end
                        na=30;
                    end
                    %-----------------------------------------
                end
                inTo=2^(i+1)-obj.ly;%Визначення кількості відліків, необхідних для інтерполяції
                Na=fix(obj.ly/na);%Визначення кількості відрізків інтерполяції
                vpls=fix(inTo/Na);%Кількість додаткових точок
                add=rem(inTo,Na);%Неповна кількість відліків, для додавання у відрізки
                Lint=vpls+na+1;
                for j=1:Na
                    nl=(j-1)*na+1;
                    kl=na*j;
                    nli=(j-1)*Lint+1;
                    kli=Lint*j;
                    [obj.xi(nli:kli),L(nli:kli)]=Lagr_rep(obj.x(nl:kl),obj.change_y(nl:kl),Lint);
                    %Видалення зайвих елементів (кінець алгоритму на 115 рядку)
                    index(nli:kli)=0;
                    for ji=nl:kl
                        for jj=nli:kli
                            if obj.xi(jj)==obj.x(ji)
                                index(jj)=1;
                            end
                        end
                    end
                    if add~=0
                        if j==add
                            xint=obj.xi;
                        end
                    else
                        xint=1;
                    end
                end
                ji=Na-add;
                j=1;
                while (ji~=0)
                    for jj=(length(xint)+(j-1)*Lint+round(Lint/10)):1:length(obj.xi)%(length(xi)-round(Lint/10)):-j:length(xint)
                        if index(jj)==0
                            obj.xi(jj)=66666666666666666.6666666666666;
                            ji=ji-1;
                            j=j+1;
                            break;
                        end
                    end
                end
                j = 1;
                while j <= length(obj.xi)
                    if obj.xi(j)==66666666666666666.6666666666666
                        obj.xi(j)=[];
                        L(j)=[];
                    else
                        j=j+1;
                    end
                end
                %-----------------------------------------
                %Апроксимація
                %-----------------------------------------
                lL=length(L);
                lxi=length(obj.xi);
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
                            p=polyfit(obj.xi(nl:kl),L(nl:kl),por_t);
                            obj.a(nl:kl)=polyval(p,obj.xi(nl:kl));
                            Do=sum((L(nl:kl)-obj.a(nl:kl)).^2);
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
                        p=polyfit(obj.xi(nl:kl),L(nl:kl),por);
                        obj.a(nl:kl)=polyval(p,obj.xi(nl:kl));
                    end
                    %-----------------------------------------
                    %Згладжування
                    %-----------------------------------------
                    lot=lxi/Na;
                    olot=round(lot/5);
                    for j=1:(Na-1)
                        nl=j*lot-olot+1;
                        kl=j*lot+olot;
                        p=polyfit(obj.xi(nl:kl),obj.a(nl:kl),5);
                        obj.a(nl:kl)=polyval(p,obj.xi(nl:kl));
                    end
                end
            else
                Ost=0;%Остача від ділення на 2^i
                Info=['Довжина сигналу складає ' num2str(obj.ly) ' відліки'];
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
                if obj.ly<=32
                    Na=4;%Кількість проміжків апроксимації
                else
                    Na=8;
                end
                Pp=7;%Максимальний порядок полінома
                if obj.ly~=obj.lx
                    disp('Довжина сигналу не рівна довжині вектора часу!')
                else
                    na=obj.lx/Na;%Кількість проміжків апроксимації
                    for j=1:Na
                        R_2_max=0;
                        por=0;
                        nl=(j-1)*na+1;
                        kl=na*j;
                        for por_t=1:Pp
                            p=polyfit(obj.x(nl:kl),obj.change_y(nl:kl),por_t);
                            obj.a(nl:kl)=polyval(p,obj.x(nl:kl));
                            Do=sum((obj.change_y(nl:kl)-obj.a(nl:kl)).^2);
                            yn1=length(obj.change_y(nl:kl));
                            ysr=sum(obj.change_y(nl:kl))/yn1;
                            Ds=sum((obj.change_y(nl:kl)-ysr).^2);
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
                        p=polyfit(obj.x(nl:kl),obj.change_y(nl:kl),por);
                        obj.a(nl:kl)=polyval(p,obj.x(nl:kl));
                    end
                    %-----------------------------------------
                    %Згладжування
                    %-----------------------------------------
                    lot=length(obj.x)/Na;
                    olot=round(lot/5);
                    for j=1:(Na-1)
                        nl=j*lot-olot+1;
                        kl=j*lot+olot;
                        p=polyfit(obj.x(nl:kl),obj.a(nl:kl),5);
                        obj.a(nl:kl)=polyval(p,obj.x(nl:kl));
                    end         
            L=obj.change_y;
            obj.xi=obj.x;
                end
            end
        %-----------------------------------------
        %Побудова графіків
        %-----------------------------------------
        obj.PohAbs=abs(L-obj.a);
        obj.PohOtn=obj.PohAbs./L*100;
       end
       
       
       function obj =  plot_signal(obj)

           figure()
            subplot(3,1,1)
            plot(obj.x,obj.change_y,obj.xi,obj.a,'--'),grid
             title('Signal graph')
             legend('Initial function','Approximated function')
             xlabel('Time, s')
            ylabel('Magnitude')
            
            
            subplot(3,1,2)
            plot(obj.xi,obj.PohAbs),grid
             title('Graph of absolute error')
             xlabel('Time, s')
             ylabel('Magnitude')
             
             
            subplot(3,1,3)
            plot(obj.xi,obj.PohOtn),grid
             title('Graph of relative error')
             xlabel('Time, s')
             ylabel('Relative error, %')
       end
    end
end