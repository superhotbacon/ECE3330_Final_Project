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
