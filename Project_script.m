[data, fs] = audioread('MahoganyHallStompLouisArmstrong.wav');
data_fft = fft(data);
plot(abs(data_fft(:,1)))

data_fft = fft(data);
plot(abs(data_fft(:,1)));
title('test-FFT');


%THIS IS TO SHOW ZACH GIT HUB IS GOATED