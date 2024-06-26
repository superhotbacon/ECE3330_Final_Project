clear all; close all;clc;

[data, fs] = audioread('Music\VINSKY V3.wav');
data_fft = abs(fft(data));
data_fft_dB = 20*log10(data_fft);
original = ifft(data_fft);
%sound(original,fs);
filename = 'Music\Output\VNSKYOUT.mp3';

T = 1/fs; %sampling period
L = size(data);
L = L(1); %only care about first column which tells the amt of samples
t = (0:L-1)*T;
fn = fs / 2; %this is the nyquest sampling maximum frequency

fc = 22000;

[b,a] = butter(6,fc/(fs/2), 'low');
dataOut = filter(b,a,data);

dataOut_fft = abs(fft(dataOut));
dataOut_fft_dB = 20*log10(abs(fft(dataOut)));


sound(dataOut,fs);
audiowrite(filename,dataOut,fs);


figure(1)
subplot(2,1,1)
x = fs/L*(0:L-1)/1000;
plot(x,data_fft,"LineWidth",1)
grid on;
xlim([0 6]);
title("Complex Magnitude of input data in Frequency Spectrum")
xlabel("f (kHz)")
ylabel("|fft(X)|")
subplot(2,1,2)

plot(x,dataOut_fft,"LineWidth",1)
grid on;
xlim([0 6]);
title("Complex Magnitude of output data in Frequency Spectrum")
xlabel("f (kHz)")
ylabel("|fft(X)|")

figure(2)
subplot(2,1,1)

x = fs/L*(0:L-1)/1000; %sampling freq / 
plot(x,data_fft_dB,"LineWidth",1)
grid on;
xlim([0 22]);
title("input data in dB")
xlabel("f (kHz)")
ylabel("|fft(X)|dB")

subplot(2,1,2)
x = fs/L*(0:L-1)/1000; %sampling freq / 
plot(x,dataOut_fft_dB,"LineWidth",1)
grid on;
xlim([0 22]);
title("Output data in dB")
xlabel("f (kHz)")
ylabel("|fft(X)|dB")


%the line below only plots half the data, cutting out the higher harmonic 
%that gets repeated because of the nyquest sampling theorem
%plot(x(1:length(x)/2),dataOut_fft_dB(1:length(dataOut_fft_dB)/2),"LineWidth",1)
%plot(fs/L*(0:L-1)/1000,X,"LineWidth",1)
figure(3)
freqz(b,a)
%from 0 to the nyquest frequency
grid on;
