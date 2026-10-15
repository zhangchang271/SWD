function [ml,saveForBackward] = RTlADx(seismo_v,df,dt,dx,np,vmin,vmax, ...
    fmin,fmax,offset,gx,sx)
%RTLADX Linear Radon transform for receivers to the left of the shot.

[uxt,~,trace_mask,distance] = offset_traces(seismo_v,gx,sx,dx,offset,'left');
ccn = fix(1/df/dt);
lf = round(fmin/df)+1;
nf = round(fmax/df)+1;
mm = zeros(np,nf);
ll0 = zeros(np,0);

if ~isempty(uxt)
    d = fft(uxt,ccn).';
    pp = 1./(vmin:vmax)';
    ll0 = 1i*2*pi*df*(pp*distance);
    for frequency_index = lf:nf
        mm(:,frequency_index) = exp(ll0*(frequency_index-1))*d(:,frequency_index);
    end
end

ml = abs(mm(:,lf:nf));
for frequency_index = 1:size(ml,2)
    scale = max(ml(:,frequency_index));
    if scale > 0
        ml(:,frequency_index) = ml(:,frequency_index)/scale;
    end
end

saveForBackward = struct('mm',mm,'ll0',ll0,'lf',lf,'nf',nf, ...
    'ccn',ccn,'mll',abs(mm(:,lf:nf)),'uxtposl',trace_mask);
end
