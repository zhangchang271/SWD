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

sides = {'right','left'};
observed = [cr_0(:),cr_0l(:)];
predicted = [cr_pre_r(:),cr_pre_l(:)];
for side = 1:2
    [data,trace_index,~,distance] = ...
        offset_traces(seismo_v,gx,sx,dx,offset,sides{side});
    if numel(trace_index) < min_traces
        continue
    end

    residual = zeros(win,1);
    residual(first_frequency:first_frequency+npair-1) = ...
        frequency .* (observed(:,side)-predicted(:,side)) ./ ...
        (observed(:,side).*predicted(:,side));
    residual(isnan(residual)) = 0;

    tempdata = zeros(ntpad,numel(trace_index));
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
    seismo_v_d1(:,trace_index) = weighted_data(1:nt,:);
end
end
