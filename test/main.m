clc; clear all; clear classes;

p = 3;
f = 20;
N = 1024;
na = 8;
[x_cl,y_cl,xi_cl,a_cl,PohAbs_cl, PohOtn_cl] = getValueClassSignal(p,f, N,  na);

[x_fn,y_fn,xi_fn,a_fn,PohAbs_fn, PohOtn_fn] = getValueMNKPlusLagr(p,f, N,  na);

if x_cl == x_fn & y_cl == y_fn & xi_cl == xi_fn & a_cl == a_fn & PohAbs_cl == PohAbs_fn & PohOtn_cl == PohOtn_fn
    a = 1
end