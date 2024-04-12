[data, fs] = audioread('CrackTheSkye.wav');
data_fft = fft(data);

original = ifft(data_fft);
%sound(original,fs);


T = 1/fs; %sampling period
L = size(data);
L = L(1); %only care about first column which tells the amt of samples
t = (0:L-1)*T;
fn = fs / 2; %this is the nyquest sampling maximum frequency


dB = 20*log10(abs(data_fft));


%fs = 96000;

fc = 1200;

[b,a] = butter(6,fc/(fs/2), 'low');
dataOut = filter(b,a,data_fft);

X = abs(ifft(dataOut));




sound(X,fs);
figure(1)
plot(fs/L*(0:L-1)/1000,X,"LineWidth",1)
title("Complex Magnitude of fft Spectrum")
xlabel("f (kHz)")
ylabel("|fft(X)|dB")


figure(2)
freqz(b,a)
