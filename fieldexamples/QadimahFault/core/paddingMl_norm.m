function [mlpre1] = paddingMl_norm(cr_pre,ml)
[np,npair] = size(ml);
mlpre1 = zeros(np,npair);
x = 1:np;
for pos_f =1:npair
mlpre1(:,pos_f) = normpdf(x,cr_pre(pos_f),16);
end

end