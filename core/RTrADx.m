function [ml,saveForBackward] = RTrADx(seismo_v,df,dt,dx,np,vmin,vmax, ...
    fmin,fmax,offset,gx,sx)
%RTrADX Linear Radon transform for receivers to the right of the shot.

[uxt,~,trace_mask,distance] = offset_traces(seismo_v,gx,sx,dx,offset,'right');
[ml,mm,ll0,lf,nf,ccn] = radon_spectrum(uxt,distance,df,dt,np,vmin,vmax,fmin,fmax);

saveForBackward = struct('mm',mm,'ll0',ll0,'lf',lf,'nf',nf, ...
    'ccn',ccn,'mlr',abs(mm(:,lf:nf)),'uxtposr',trace_mask);
end

function [ml,mm,ll0,lf,nf,ccn] = radon_spectrum(uxt,distance,df,dt,np,vmin,vmax,fmin,fmax)
ccn = fix(1/df/dt);
lf = round(fmin/df)+1;
nf = round(fmax/df)+1;
mm = zeros(np,nf);

if isempty(uxt)
    ml = zeros(np,nf-lf+1);
    ll0 = zeros(np,0);
    return
end

d = fft(uxt,ccn).';
pp = 1./(vmin:vmax)';
ll0 = 1i*2*pi*df*(pp*distance);
for frequency_index = lf:nf
    mm(:,frequency_index) = exp(ll0*(frequency_index-1))*d(:,frequency_index);
end

ml = abs(mm(:,lf:nf));
for frequency_index = 1:size(ml,2)
    scale = max(ml(:,frequency_index));
    if scale > 0
        ml(:,frequency_index) = ml(:,frequency_index)/scale;
    end
end
end
