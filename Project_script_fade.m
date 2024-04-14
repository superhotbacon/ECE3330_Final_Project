clear all; close all;clc;

[data, fs] = audioread('Music\CrackTheSkye.wav');
data_fft = abs(fft(data));
data_fft_dB = 20*log10(data_fft);
original = ifft(data_fft);
%sound(original,fs);


T = 1/fs; %sampling period
L = size(data);
L = L(1); %only care about first column which tells the amt of samples
t = (0:L-1)*T;
fn = fs / 2; %this is the nyquest sampling maximum frequency

t_cutin = 7; %in seconds
samples_per_segment = 100;

num_samples = (t_cutin * fs);

data_cutIn = data(1:num_samples,:); %data for number of 
%samples needing to be filterd
fc_start = 100;
fc_end = 22000;
%FOR LINEAR INCREMENT
%fc_increment = (fc_end - fc_start)/(num_samples/samples_per_segment);

%for linear
%{
%first filter assignment
[b,a] = butter(6,fc_start/(fs/2), 'low');
[data_out_cutIn(1:samples_per_segment,:),zf] =...
    filter(b,a,data_cutIn(1:samples_per_segment,:));

fc_current = fc_start + fc_increment;
for sample_seg = samples_per_segment:1:num_samples/samples_per_segment - 1
    [b,a] = butter(6,(fc_current)/(fs/2), 'low');
    [data_out_cutIn(samples_per_segment*sample_seg:samples_per_segment*(sample_seg+1),:),zf] = filter(b,a,data_cutIn(samples_per_segment*sample_seg:samples_per_segment*(sample_seg+1),:),zf);
    
    fc_current = fc_current + fc_increment;
end
%}

%for quadratic
%quadratic increment
num_segments = (num_samples/samples_per_segment);
%scalar value
k = 1; %implement this
segment_nums(1:(num_segments)) = (1:1:num_segments);
a = (fc_end-fc_start)/((num_segments-1)^2);%a value for quad equation
%quadratic equation of fc_increment
fc_increment(1:(num_segments)) = a.*(segment_nums(1:num_segments)).^2;
    %k*(segment_nums.^2 - segment_nums.*(1 + num_segments) + 1*num_segments); %using num_segments as 'x' value 

%first filter assignment
[b,a] = butter(6,fc_start/(fs/2), 'low');
[data_out_cutIn(1:samples_per_segment,:),zf] =...
    filter(b,a,data_cutIn(1:samples_per_segment,:));

fc_current =  fc_increment(1);
for sample_seg = samples_per_segment:1:num_samples/samples_per_segment - 1
    [b,a] = butter(6,(fc_current)/(fs/2), 'low');
    [data_out_cutIn(samples_per_segment*sample_seg:samples_per_segment*(sample_seg+1),:),zf] = filter(b,a,data_cutIn(samples_per_segment*sample_seg:samples_per_segment*(sample_seg+1),:),zf);
    
    fc_current = fc_increment(sample_seg);
end



dataOut = [data_out_cutIn(:,:); data(t_cutin*fs:L,:)];%cconacntoante it 
%data(num_samples+samples_per_segment*2 is needed for an array offset

dataOut_fft = abs(fft(dataOut));
dataOut_fft_dB = 20*log10(abs(fft(dataOut)));






sound(dataOut,fs);


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


%following plot prints error
%plot(x,dataOut_fft,"LineWidth",1)
%grid on;
%xlim([0 6]);
%title("Complex Magnitude of output data in Frequency Spectrum")
%xlabel("f (kHz)")
%ylabel("|fft(X)|")


%gives error
%figure(2)
%subplot(2,1,1)
%x = fs/L*(0:L-1)/1000; %sampling freq / 
%plot(x,dataOut_fft_dB,"LineWidth",1)
%grid on;
%xlim([0 6]);
%title("Output data in dB")
%xlabel("f (kHz)")
%ylabel("|fft(X)|dB")

subplot(2,1,2)
x = fs/L*(0:L-1)/1000; %sampling freq / 
plot(x,data_fft_dB,"LineWidth",1)
grid on;
xlim([0 6]);
title("input data in dB")
xlabel("f (kHz)")
ylabel("|fft(X)|dB")

%the line below only plots half the data, cutting out the higher harmonic 
%that gets repeated because of the nyquest sampling theorem
%plot(x(1:length(x)/2),dataOut_fft_dB(1:length(dataOut_fft_dB)/2),"LineWidth",1)
%plot(fs/L*(0:L-1)/1000,X,"LineWidth",1)
figure(3)
freqz(b,a)
grid on;
