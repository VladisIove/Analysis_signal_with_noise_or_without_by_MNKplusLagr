clc
clear all
%-----------------------------------------
%Функція
%-----------------------------------------
for f=[2, 20, 200]%Частота періодичного сигналу
    T=1/f;%Період сигналу
    p=6;%Кількість періодів у сигналі
    N=512;%Кількість відліків (точок) у сигналі та Кількість точок для інтерполяції
    Td=p*T/N;%Період дискретизації
    fd=1/Td;%Частота дискретизації
    x=Td:Td:p*T;%Вектор часу (Вектор аргументів початкової функції з 'N' точок)
    y=3+sin(2*pi*x*f).*exp(-7*x);%Вектор значень функції з 'N' точок
    ly=length(y);%Довжина сигналу
    lx=length(x);%Довжина вектора часу
    Y_noise_vector = function_add_violet_noise(1,ly);
    koef = find_koef(y,Y_noise_vector, 2, 0.22);
    y = y + Y_noise_vector*koef;
    [a,L,xi]=MNK_Lagr(x,y, N)
            %-----------------------------------------
            %Побудова графіків
            %-----------------------------------------
            PohAbs=abs(L-a);
            PohOtn=PohAbs./L*100;
            figure()
            subplot(2,1,1)
            plot(x,y),grid
    %         title('Signal graph')
    %         legend('Initial function','Approximated function')
    %         xlabel('Time, s')
    %         ylabel('Magnitude')
            title('Графік сигналу')
            legend('Початкова функція')
            xlabel('Час, с')
            ylabel('Амплітуда')
                        subplot(2,1,2)
            plot(xi,a),grid
    %         title('Signal graph')
    %         legend('Initial function','Approximated function')
    %         xlabel('Time, s')
    %         ylabel('Magnitude')
            title('Графік сигналу')
            legend('Апроксимована функція')
            xlabel('Час, с')
            ylabel('Амплітуда')
            
            
            figure()
            subplot(2,1,1)
            plot(xi,PohAbs),grid
    %         title('Graph of absolute error')
    %         xlabel('Time, s')
    %         ylabel('Magnitude')
            title('Графік абсолютної похибки')
            xlabel('Час, с')
            ylabel('Амплітуда')
            subplot(2,1,2)
            plot(xi,PohOtn),grid
    %         title('Graph of relative error')
    %         xlabel('Time, s')
    %         ylabel('Relative error, %')
            title('Графік відносної похибки')
            xlabel('Час, с')
            ylabel('Відносна похибка, %')
  
end