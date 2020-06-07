function a=MNK_Part(x, y)
%y - Функція для апроксимації
%a - апроксимована функція
%t - вектор часу
%N - кількість відрізків апроксимації
%Pp_max - максимальний порядок полінома апроксимації
%n - кількість точок на один відрізок
%R_2 - коефіцієнт детермінації
%PohAbs - абсолютна похибка апроксимації
%PohOtn - відносна похибка апроксимації
%     %-----------------------------------------
%     %Апроксимація
%     %-----------------------------------------
ly=length(y);
lx=length(x);
if ly<=32
    Na=4;%Кількість проміжків апроксимації
else
    Na=8;
end
Pp=5;%Максимальний порядок полінома
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
    lot=lx/Na;
    olot=round(lot/5);
    for j=1:(Na-1)
        nl=j*lot-olot+1;
        kl=j*lot+olot;
        p=polyfit(x(nl:kl),a(nl:kl),5);
        a(nl:kl)=polyval(p,x(nl:kl));
    end
end