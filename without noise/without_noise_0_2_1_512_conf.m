function [f,p,N,y,Y_noise,x,ly,lx] = without_noise_0_2_1_512_conf()
    N = 512;
    p = 1;
    f = 2; 
    T=1/f;
    Td=p*T/N;
    fd=1/Td;
    x=Td:Td:p*T;
    y=3+sin(2*pi*x*f).*exp(-7*x);
    ly=length(y);
    lx=length(x);
    Y_noise=0;
end

