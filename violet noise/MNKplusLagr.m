clc
clear all
%-----------------------------------------
%�������
%-----------------------------------------
for f=[2, 20, 200]%������� ����������� �������
    T=1/f;%����� �������
    p=6;%ʳ������ ������ � ������
    N=512;%ʳ������ ����� (�����) � ������ �� ʳ������ ����� ��� ������������
    Td=p*T/N;%����� �������������
    fd=1/Td;%������� �������������
    x=Td:Td:p*T;%������ ���� (������ ��������� ��������� ������� � 'N' �����)
    y=3+sin(2*pi*x*f).*exp(-7*x);%������ ������� ������� � 'N' �����
    ly=length(y);%������� �������
    lx=length(x);%������� ������� ����
    Y_noise_vector = function_add_violet_noise(1,ly);
    koef = find_koef(y,Y_noise_vector, 2, 0.22);
    y = y + Y_noise_vector*koef;
    [a,L,xi]=MNK_Lagr(x,y, N)
            %-----------------------------------------
            %�������� �������
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
            title('������ �������')
            legend('��������� �������')
            xlabel('���, �')
            ylabel('��������')
                        subplot(2,1,2)
            plot(xi,a),grid
    %         title('Signal graph')
    %         legend('Initial function','Approximated function')
    %         xlabel('Time, s')
    %         ylabel('Magnitude')
            title('������ �������')
            legend('������������� �������')
            xlabel('���, �')
            ylabel('��������')
            
            
            figure()
            subplot(2,1,1)
            plot(xi,PohAbs),grid
    %         title('Graph of absolute error')
    %         xlabel('Time, s')
    %         ylabel('Magnitude')
            title('������ ��������� �������')
            xlabel('���, �')
            ylabel('��������')
            subplot(2,1,2)
            plot(xi,PohOtn),grid
    %         title('Graph of relative error')
    %         xlabel('Time, s')
    %         ylabel('Relative error, %')
            title('������ ������� �������')
            xlabel('���, �')
            ylabel('³������ �������, %')
  
end