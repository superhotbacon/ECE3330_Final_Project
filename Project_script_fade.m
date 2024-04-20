clear all; close all;clc;

[data, fs] = audioread('Music\HighRoad.mp3');

T = 1/fs; %sampling period
L = size(data);
L = L(1); %only care about first column which tells the amt of samples
t = (0:L-1)*T;
fn = fs / 2; %this is the nyquest sampling maximum frequency
filename = 'Music\Output\HighRoadOutput.mp3';

t_delta = 5; 
% INPUT PARAMETER
%Desired fade in time
samples_per_segment = 200;
% INPUT PARAMETER
% Samples per segment determines how many audio
% samples are filtered per incremental filter.

num_samples = (t_delta * fs);
%total number of samples to filter

data_cutIn = data(1:num_samples,:); 
%data for number of samples needing to be filterd
data_cutIn_fft = abs(fft(data_cutIn));
data_cutIn_fft_dB = 20*log10(data_cutIn_fft);

fc_start = 200;
fc_end = 22000;
% INPUT PARAMETERS
% starting and ending corner frequency cutoff values

%USING QUADRATIC INCREMENT
%quadratic increment
num_segments = floor(num_samples/samples_per_segment);
%scalar value
segment_nums(1:(num_segments)) = (1:1:num_segments);
a = (fc_end-fc_start)/((num_segments-1)^2);
% 'a' value for quadratic equation

fc_increment(1:(num_segments)) = a.*(segment_nums(1:num_segments)).^2 + fc_start; 
% quadratic equation of fc_increment using input parameters

%first filter assignment
[b,a] = butter(6,fc_start/(fs/2), 'low');
[data_out_cutIn(1:samples_per_segment,:),zf] =...
    filter(b,a,data_cutIn(1:samples_per_segment,:));

fc_current =  fc_increment(1);

%The following for loop iterates through the number of total segments
%needed to reach t_delta. It iterates through and applies the incremental
%low pass filters until it reaches 22000 Hz, which is about the max the 
%human ear can hear.
for sample_seg = 2:1:num_segments -1
    [b,a] = butter(6,(fc_current)/(fs/2), 'low');
    [data_out_cutIn(samples_per_segment*sample_seg+1:samples_per_segment*(sample_seg+1),:),zf]...
    = filter(b,a,data_cutIn(samples_per_segment*sample_seg+1:samples_per_segment*(sample_seg+1),:),zf);
    % samples_per_segment*sample_seg+1 is nessesary to prevent
    % the for loop from filtering the same sample twice, 
    % prevents popping noise
    fc_current = fc_increment(sample_seg);
    % set the next filter frequency based on the matrix defined above
end

dataOut = [data_out_cutIn(:,:); data(num_segments*samples_per_segment:L,:)];%cconacntoante it 
% concanate the filtered data and the rest of the song

dataOut_fft = abs(fft(data_out_cutIn));
dataOut_fft_dB = 20*log10(abs(fft(data_out_cutIn)));

sound(dataOut,fs);
audiowrite(filename,dataOut,fs);

figure(1)
subplot(2,1,1)
x = fs/L*(0:L-1)/1000; %frequency, in kHz
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
plot(x,data_cutIn_fft_dB,"LineWidth",1)
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
grid on;


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
