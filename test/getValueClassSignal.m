function [x,y,xi,a,PohAbs, PohOtn] = getValueClassSignal(p,f, N,  na)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    signal = SignalClass();
    signal.p = p; % count periods
    signal.f = f; % frequency
    signal.N = N; % count dot
    %signal.snr = snr; % signal to noise ratio; optional value
    signal.na = na% interpoletion
    %% STEP 2
    % Create first configuration 
    % Need use if you have reset data
    signal.updateValue();
    %% STEP 3
    % Add noise to signal. Optional step
    %signal.add_noise("white"); % in example use white noise. P.s. check STEP 0
    %% STEP 4
    % Use MNKplusLagr
    signal.MNKplusLagr();
    
    x = signal.x;
    y = signal.change_y;
    xi = signal.xi;
    a = signal.a;
    PohAbs = signal.PohAbs;
    PohOtn = signal.PohOtn;
    
end

