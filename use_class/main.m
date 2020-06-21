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

Info= {};
answer = [];
for name = ["","white", "pink", "grey", "broun", "blue", "violet"]
    answer_signal = [];
    Info(end+1)={name};
    for s = [0.5, 1, 2]
        for i = [1,2,3,4,5,6]
            for j = [512,1024,2048,4096,8192]
                for k = [2,20,200]
                    signal = SignalClass();
                    signal.p = i; % count periods
                    signal.f = k; % frequency
                    signal.N = j; % count dot
                    signal.snr = s; % signal to noise ratio; optional value
                    signal.na = 15% interpoletion
                    %% STEP 2
                    % Create first configuration 
                    % Need use if you have reset data
                    signal.updateValue();
                    signal.add_noise(name);
                    signal.MNKplusLagr();
                    
                    m = min(signal.PohOtn);
                    text = ['snr ' num2str(signal.snr) ' p ' num2str(signal.p) ' N ' num2str(signal.N) ' f ' num2str(signal.f) ' min ' num2str(m)];
                    Info(end+1)={text};
                    answer_signal(end+1) = m;
                end
            end
        end
    end
    answer = [answer; answer_signal];
end
    for i = 1:length(Info)
        disp(Info(i));
    end

%{
signal = SignalClass();
signal.p = 1; % count periods
signal.f = 2; % frequency
signal.N = 512; % count dot
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

%}
