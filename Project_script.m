[data, fs] = audioread('CrackTheSkye.wav');
data_fft = fft(data);

original = ifft(data_fft);
%sound(original,fs);


T = 1/fs; %sampling period
L = size(data);
L = L(1); %only care about first column which tells the amt of samples
t = (0:L-1)*T;
fn = fs / 2; %this is the nyquest sampling maximum frequency


%dB = 20*log10(abs(data_fft));




fc = 1200;

[b,a] = butter(6,fc/(fs/2), 'low');
dataOut = filter(b,a,data);

dataOut_fft = abs(fft(dataOut));
dataOut_fft_dB = 20*log10(abs(fft(dataOut)));


%sound(dataOut,fs);

figure(1)
x = fs/L*(0:L-1)/1000;
plot(x,dataOut_fft,"LineWidth",1)
xlim([0 20]);
title("Complex Magnitude in Frequency Spectrum")
xlabel("f (kHz)")
ylabel("|fft(X)|dB")

figure(2)
x = fs/L*(0:L-1)/1000; %sampling freq / 
plot(x,dataOut_fft_dB,"LineWidth",1)
xlim([0 20]);
title("Output data in dB")
xlabel("f (kHz)")
ylabel("|fft(X)|dB")


%plot(fs/L*(0:L-1)/1000,X,"LineWidth",1)
figure(3)
freqz(b,a)
