clc; clear all; clear classes;
%% Create signal with noise or without
% STEP 0
% NOISES USE IN PROJECT
% white noise -> "white"
% pink noise -> "pink"
% grey noise -> "grey"
% broun noise -> "broun"
% blue noise -> "blue"
% without noise -> DOSNT use method add_noise
%% STEP 1
% Create object with need starte parameters
signal = SignalClass();
signal.p = 1; % count periods
signal.f = 2; % frequency
signal.N = 256; % count dot
signal.snr = 1/2; % signal to noise ratio; optional value
signal.na = 15% interpoletion
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
%% STEP 5
% Plot
signal.plot_signal();
