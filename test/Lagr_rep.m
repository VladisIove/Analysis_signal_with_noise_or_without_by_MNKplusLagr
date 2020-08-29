function [x,L]=Lagr_rep(xp,yp,N)
%xp - pochatkoviu vektor argymentiv
%yp - pochatkoviu vektor znacheniu fun
%N - kilkist tochok y interpolovanomy vektori
%x - vektor argymentiv dla interpolatcii
%L - vektor znacheniu interpolatcii fun
n=length(xp);% Dovgina pochatkovogo vektora
x=[];% Pochatkove znach vektory
if N<n
    disp('Kilkist tochok interpolatcii ne moge byti menshoy za kilkist ')
    disp('tochok signaly')
else
    dx=(xp(n)-xp(1))/(N-1);
    x=xp(1):dx:xp(n);
    for i=1:n
        for j=2:length(x)-1
            if (x(j+1)>xp(i) && x(j)<xp(i)) 
                k1=x(j+1)-xp(i);
                k2=xp(i)-x(j);
                if k1<k2
                    x(j+1)=xp(i);
                else
                    x(j)=xp(i);
                end
            end
        end
    end
    for l=1:N
        %-----------------------------------------
        % algoritm interpolatcii
        %-----------------------------------------
        for i=1:n
            P(i)=1;
            for j=1:n
                if j~=i
                    P(i)=P(i)*(x(l)-xp(j))/(xp(i)-xp(j));
                    % rozraxynok chiselnika ta znamenika
                end
            end
        end
        L(l)=sum(P.*yp);% Viznachena
        %-----------------------------------------
        % Kinetc algoritmy interpolatcii
        %-----------------------------------------
    end
end