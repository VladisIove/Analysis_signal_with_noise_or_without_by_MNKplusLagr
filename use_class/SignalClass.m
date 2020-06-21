classdef SignalClass  < handle
    properties
        % ------------------------------
        %    Gloabal properties  object
        % ------------------------------
        N {mustBeNumeric};
        p {mustBeNumeric};
        f {mustBeNumeric};
        ly {mustBeNumeric};
        lx {mustBeNumeric};
        na {mustBeNumeric};
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
           % ------------------------------
           %    Update value it need if we change input value 
           %    and need update object
           % ------------------------------
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
           % ------------------------------
           %    FIND NEED KOEFISIENT BY SIGNAL TO NOISE RATIO
           % ------------------------------
           koef = sqrt( mean(signal.^2) / ( 10^(snr/10)* mean(noise.^2) ) );
       end

       function obj = add_noise(obj, name_noise)
           % ------------------------------
           %    ADD NOISE TO SIGNAL
           % ------------------------------
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
       
       function obj = write_message(obj, i)
          Info=['Dovgina signaly skladae ' num2str(obj.ly) ' vidliki'];
          disp('-----------------------------------------------------')
          disp(Info)
          disp('-----------------------------------------------------')
          Info=['Naibligha stepin dviuki, i=' num2str(i)];
          disp(Info)
          disp('-----------------------------------------------------')
          Ost=obj.ly-2^i;
          Info=['Ostacha dilena na 2^i=' num2str(Ost)];
          disp(Info)
          disp('-----------------------------------------------------')
       end
       
       function obj = interpoletion(obj)
            while (rem(obj.lx,obj.na)~=0)
                obj.na=obj.na+1;%Визначення дільника, більшого за 15
                %-----------------------------------------
                %Зменшення кількості початкових відліків до ділення на 32
                %-----------------------------------------
                if obj.na>40
                    ji=32*floor(obj.lx/32);
                    jj=floor(obj.lx/(obj.lx-ji))-1;
                    while(obj.lx~=ji)
                        obj.x(jj)=[];
                        obj.change_y(jj)=[];
                        obj.lx=length(obj.x);
                        obj.ly=length(obj.change_y);
                    end
                    obj.na=30;
                end
            end
       end
       
       function [por] = find_poradok_polinoma(obj, xi, polinom,L, nl, kl)
        R_2_max=0;
        por=0;
        for por_t=1:polinom
            p=polyfit(xi(nl:kl),L(nl:kl),por_t);
            obj.a(nl:kl)=polyval(p,xi(nl:kl));
            Do=sum((L(nl:kl)-obj.a(nl:kl)).^2);
            Ds=sum((L(nl:kl)-mean(L(nl:kl))).^2); % dispercia sluchanoi velichini !mat ochikuvanie == 1?
            if Ds==0
                disp('Dilena no nol! Zminit vhidni danni!')
            else
                if Do<1e-7
                    Do=0;
                end
                R_2=1-Do/Ds;%Viznechena determinatcii
                if R_2_max<R_2
                    R_2_max=R_2;
                    por=por_t;
                end
            end
        end
       end
       
       function obj = zglagivanie(obj, Na,lxi, xi)
            lot=lxi/Na;
            olot=round(lot/5); 
            for j=1:(Na-1)
                nl=j*lot-olot+1;
                kl=j*lot+olot;
                p=polyfit(xi(nl:kl),obj.a(nl:kl),5);
                obj.a(nl:kl)=polyval(p,xi(nl:kl));
            end
       end
       function [L] = aproximate(obj, polinom, L, xi)
           lL=length(L);
           lxi=length(xi);
           if lL<=32
            Na=4;%Kilkist promogitkiv aproksimatcii
           else
            Na=8;
           end
           if lL~=lxi
            disp('Dovzina signaly ne rivna dovgini vectory chacy!')
           else
            obj.na=lxi/Na;% Kilkist promigkiv aproksimatcii
            for j=1:Na
                nl=(j-1)*obj.na+1;
                kl=obj.na*j;
                por = obj.find_poradok_polinoma(xi, polinom,L, nl, kl)
                p=polyfit(xi(nl:kl),L(nl:kl),por);
                obj.a(nl:kl)=polyval(p,xi(nl:kl));
            end
            obj.zglagivanie(Na,lxi, xi); %Zglagyvanie           
          end
       end
       
       function obj =  MNKplusLagr(obj)
           if isempty(obj.Y_noise)
               obj.change_y = obj.y
           end
           i = floor(log2(obj.ly)); % power of 2
            if obj.ly-2^i~=0%ymova nauavnosti ostachi
                obj.write_message(i)
                obj.interpoletion(); % Interpoletion
                
                inTo=2^(i+1)-obj.ly;% Viznachenie kilkosti vidlikiv dla interpolatcii
                Na=fix(obj.ly/obj.na);% Viznachenie kilkosti vidrizkiv dla inerpolatcii
                vpls=fix(inTo/Na);% Kilkist dodatkovix tochok
                add=rem(inTo,Na);% Nepovna kilkist vidlikiv, dla dodavanna y vidrizky
                Lint=vpls+obj.na+1;
                
                for j=1:Na
                    nl=(j-1)*obj.na+1;
                    kl=obj.na*j;
                    nli=(j-1)*Lint+1;
                    kli=Lint*j;
                    [obj.xi(nli:kli),L(nli:kli)]=Lagr_rep(obj.x(nl:kl),obj.change_y(nl:kl),Lint);
                    % Vidalena zaivih elementiv(Kinetc 115 radok)
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
                L = obj.aproximate(5, L, obj.xi); 
            else
                obj.write_message(i);
                L = obj.aproximate(7,obj.change_y, obj.x);
                obj.xi = obj.x;
            end
        %-----------------------------------------
        %Find absolute and otnositelnie pohibki
        %-----------------------------------------
        obj.PohAbs=abs(L-obj.a);
        obj.PohOtn=obj.PohAbs./abs(L)*100;
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