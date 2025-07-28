function [a_data_res] = ADWDgrad3(nt,ng,ns,npair,is,w,m,M,SoftArgNorm,grad_outputr,grad_outputl,space_M,saveForBackwardr,saveForBackwardl)
% back propagation wavelet

a_data_res = zeros(nt,ng);
% ll0: Spatial position matrix, depending only on gx and frequency range
% lf nf: frequency range
% ccn: length FFT
% uxtposr: common shot gather
% offsetr: length of CSG
% mlr: dispersion spectrum
% mmr: dispersion spectrum

if is<=ns-floor(m/M)-round(w/M)
    ll0 = saveForBackwardr.ll0;%Loading Variables from the forward Process
    lf = saveForBackwardr.lf;
    nf = saveForBackwardr.nf;
    ccn = saveForBackwardr.ccn;
    uxtposr = saveForBackwardr.uxtposr;
    offsetr = length(uxtposr);
    mlr = saveForBackwardr.mlr;
    mlr = (mlr).^(1/40);
    cr_r = saveForBackwardr.cr_r;
%     [mlr] = paddingMl_norm(cr_r,mlr);
    mmr = saveForBackwardr.mm;
    ml_Nr = softmax(mlr*SoftArgNorm);
    gradmlr = zeros(size(ml_Nr));
    mmr = mmr(:,lf:nf);
    grad_outputr(grad_outputr==0)=1;
    %$\frac{\partial\kappa}{\partial D(k_i,\omega)}=\beta E_i(i-\sum_{j=1}^njE_j)$
    for i = 1:npair
        test3r = ml_Nr(:,i).*space_M(:,i) - sum((space_M(:,i).*ml_Nr(:,i))'.*(ml_Nr(:,i)),2);
        gradmlr(:,i) = test3r .* 2 .* (grad_outputr(i)).*SoftArgNorm;%*mask_pre(i)*mask_obs(i)
    end
%     gradmlr = norm_trace(gradmlr);
    gradmmr = mmr./abs(mmr);
    gradmmr(isnan(gradmmr)) = 0;
    gradmmr = gradmmr .* gradmlr;
    graddr = zeros(offsetr,ccn);
    for luoj = lf:nf
            temps = exp(ll0*(luoj-1))'*gradmmr(:,luoj-lf+1);
            graddr(:,luoj) = temps;%.*(nf-luoj)./25; 
        
    end
% for i=lf-25:lf
%     temps = graddr(:,lf);
%     graddr(:,i) = temps*((i-lf+25)/25);
% %     graddr(:,i) = temps;
% %     graddr(:,i) = temps.*0.5.*((cr_r(1)-1200)./grad_outputr(1));
% end
% for i=nf:nf+25
%     temps = graddr(:,nf);
%     graddr(:,i) = temps*((nf-i+25)/25);
% end
    
    graduxtr = real(ifft(graddr,ccn,2));
    graduxtr = (graduxtr(:,1:nt));
    a_data_res(:,uxtposr) = graduxtr';
end

if is>=round(w/M)+floor(m/M)+1 
    ll0 = saveForBackwardl.ll0;
    lf = saveForBackwardl.lf;
    nf = saveForBackwardl.nf;
    ccn = saveForBackwardl.ccn;
    uxtposl = saveForBackwardl.uxtposl;
    offsetl = length(uxtposl);
    mll = saveForBackwardl.mll;
    mll = mll.^(1/40);
    cr_l = saveForBackwardl.cr_l;
%     [mll] = paddingMl_norm(cr_l,mll);
    mml = saveForBackwardl.mm;
    ml_Nl = softmax(mll*SoftArgNorm);
    gradmll = zeros(size(ml_Nl));
    mml = mml(:,lf:nf);
    grad_outputl(grad_outputl==0)=1;
    for i = 1:npair
        test3l = ml_Nl(:,i).*space_M(:,i) - sum((space_M(:,i).*ml_Nl(:,i))'.*(ml_Nl(:,i)),2);
        gradmll(:,i) = test3l .* 2 .* (grad_outputl(i)).*SoftArgNorm;%*mask_pre(i)*mask_obs(i)
    end
%     gradmll = norm_trace(gradmll);
    gradmml = mml./abs(mml);
    gradmml(isnan(gradmml)) = 0;
    gradmml = gradmml .* gradmll;
    graddl = zeros(offsetl,ccn);
    for luoj = lf:nf
            temps = exp(ll0*(luoj-1))'*gradmml(:,luoj-lf+1);
            graddl(:,luoj) = temps;%.*(nf-luoj)./25; 
    end
% for i=lf-25:lf
%     temps = graddl(:,lf);
%     graddl(:,i) = temps*((i-lf+25)/25);
% %     graddl(:,i) = temps;
% %     graddl(:,i) = temps.*0.5.*((cr_l(1)-1200)./grad_outputl(1));
% end
% for i=nf:nf+25
%     temps = graddl(:,nf);
%     graddl(:,i) = temps*((nf-i+25)/25);
% end

    graduxtl = real(ifft(graddl,ccn,2));
    graduxtl = (graduxtl(:,1:nt));
    a_data_res(:,uxtposl(1:end-1)) = graduxtl(end:-1:2,:)';
end

a_data_res = a_data_res./max(abs(a_data_res(:)));

end

