[data, fs] = audioread('CartoonShot.wav');

T = 1/fs; %sampling period
L = size(data);
L = L(1); %only care about first column which tells the amt of samples
t = (0:L-1)*T;
fn = fs / 2; %this is the nyquest sampling maximum frequency


plot(fs/L*(0:L-1),20*log(abs(data_fft)),"LineWidth",1)
title("Complex Magnitude of fft Spectrum")
xlabel("f (Hz)")
ylabel("|fft(X)|dB")


%plot(data(:,:))

%Fs = 1000;            % Sampling frequency                    
%T = 1/Fs;             % Sampling period       
%%L = 1500;             % Length of signal
%t = (0:L-1)*T;        % Time vecto

%From Internet
%THIS IS TO SHOW ZACH GIT HUB IS GOATED
%clear sound;
%stops the thing from playing