function seismo_v_d1 = weight_dataAD(seismo_v,dt,dx,df,offset,sx,gx,min_traces, ...
    fmin,fmax,cr_0,cr_0l,cr_pre_r,cr_pre_l)
%WEIGHT_DATAAD Build the WD adjoint data from physical offset distances.
%
% The receiver selection is shared with RTrADx/RTlADx through
% OFFSET_TRACES.  OFFSET is expressed in the same physical unit as DX.

[nt,ng] = size(seismo_v);
seismo_v_d1 = zeros(nt,ng);

frequency = 2*pi*(fmin:df:fmax)';
first_frequency = round(fmin/df)+1;
npair = numel(frequency);
ntpad = floor(1/dt/df);
win = npair + 2*first_frequency;

[right_data,~,right_mask,right_distance] = ...
    offset_traces(seismo_v,gx,sx,dx,offset,'right');
if sum(right_mask) >= min_traces
    residual = zeros(win,1);
    residual(first_frequency:first_frequency+npair-1) = ...
        frequency .* (cr_0-cr_pre_r) ./ (cr_0.*cr_pre_r);
    residual(isnan(residual)) = 0;
    seismo_v_d1(:,right_mask) = apply_weight( ...
        right_data,right_distance,residual,nt,ntpad);
end

[left_data,~,left_mask,left_distance] = ...
    offset_traces(seismo_v,gx,sx,dx,offset,'left');
if sum(left_mask) >= min_traces
    residual = zeros(win,1);
    residual(first_frequency:first_frequency+npair-1) = ...
        frequency .* (cr_0l-cr_pre_l) ./ (cr_0l.*cr_pre_l);
    residual(isnan(residual)) = 0;
    weighted_left = apply_weight(left_data,left_distance,residual,nt,ntpad);
    seismo_v_d1(:,left_mask) = weighted_left(:,end:-1:1);
end
end

function weighted_data = apply_weight(data,distance,residual,nt,ntpad)
% Apply the frequency-domain wavenumber perturbation to one side.

offset_count = numel(distance);
tempdata = zeros(ntpad,offset_count);
tempdata(1:nt,:) = data;
spectrum = fft(tempdata);
kernel = (-1i/pi/2) * (residual * distance(:)');

for frequency_index = 1:numel(residual)
    positive_index = frequency_index + 1;
    negative_index = ntpad + 1 - frequency_index;
    spectrum(positive_index,:) = spectrum(positive_index,:) .* kernel(frequency_index,:);
    spectrum(negative_index,:) = conj(spectrum(positive_index,:));
end

weighted_data = ifft(spectrum);
weighted_data = weighted_data(1:nt,:);
end
