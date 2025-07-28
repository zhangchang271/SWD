function [mlpre1,mm1] = paddingMl_w2(cr_pre,ml,mm)
[np,npair] = size(ml);
% [~,npair_sub] = size(cr_pre);
pos_sub1 = cr_pre(2,1) - 1;
curve_um = zeros(npair,1);
curve_dm = zeros(npair,1);
for pos_f = cr_pre(2,:)
    [~,pfpk1] = findpeaks(-ml(:,pos_f));
    pfpk2 = [1;pfpk1;np];
    % [~,locs_m] = max(mlr(:,pos_f));
    d_loc = find((pfpk2-cr_pre(1,pos_f-pos_sub1))>0,1);
    if isempty(d_loc)
        curve_um(pos_f) = cr_pre(1,pos_f-pos_sub1);
        curve_dm(pos_f) = cr_pre(1,pos_f-pos_sub1);
    else
        curve_um(pos_f) = pfpk2(d_loc);
        curve_dm(pos_f) = pfpk2(d_loc-1);
    end


end
mlpre1 = zeros(np,npair);
mm1 = zeros(np,npair);
for pos_f =cr_pre(2,:)
    mlpre1(curve_dm(pos_f):curve_um(pos_f),pos_f) = ...
        ml(curve_dm(pos_f):curve_um(pos_f),pos_f);
    mm1(curve_dm(pos_f):curve_um(pos_f),pos_f) = ...
     mm(curve_dm(pos_f):curve_um(pos_f),pos_f);
end
end