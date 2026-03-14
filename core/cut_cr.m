function cr_out = cut_cr(cr_0,freq,offset,fillv)
freq = freq(:);
[nf,ncol] = size(cr_0);
if isscalar(offset)
    offset = repmat(offset,1,ncol);
end
cr_out = cr_0;
for i = 1:ncol
    c = cr_0(:,i);
    if all(c==fillv)
        cr_out(:,i) = nan;
        continue
    end
    d = c - freq*offset(i)/1.8;
    k = find(d(1:end-1).*d(2:end)<=0,1,'first');
    if ~isempty(k)
        cr_out(1:k,i) = nan;
    end
end
end