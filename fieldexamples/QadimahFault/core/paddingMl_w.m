function [mlpre1,mm1] = paddingMl_w(cr_pre,ml,mm)
[np,npair] = size(ml);
curve_um = zeros(npair,1);
curve_dm = zeros(npair,1);
for pos_f =1:npair
    [~,pfpk1] = findpeaks(-ml(:,pos_f));
    pfpk2 = [1;pfpk1;np];
    % [~,locs_m] = max(mlr(:,pos_f));
    d_loc = find((pfpk2-cr_pre(pos_f))>0,1);
    if isempty(d_loc)
        curve_um(pos_f) = cr_pre(pos_f,1);
        curve_dm(pos_f) = cr_pre(pos_f,1);
    else
        curve_um(pos_f) = pfpk2(d_loc);
        curve_dm(pos_f) = pfpk2(d_loc-1);
    end
end
mlpre1 = zeros(np,npair);
mm1 = zeros(np,npair);
for pos_f =1:npair
    mlpre1(curve_dm(pos_f):curve_um(pos_f),pos_f) = ...
        ml(curve_dm(pos_f):curve_um(pos_f),pos_f);
    mm1(curve_dm(pos_f):curve_um(pos_f),pos_f) = ...
     mm(curve_dm(pos_f):curve_um(pos_f),pos_f);
end
end