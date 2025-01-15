
function [res,cr_pre] = LHDispPick(ml,vmin,cr_0,method,ini)

    if method==1
        [idx,~,~]=pickFDC1(ml,ini,2,4);
        cr_pre = idx' + vmin;
    else
        [~,p_max] = max(ml);
        cr_pre = p_max + vmin;

    end
    res = mean(((cr_pre - cr_0).^2));
    
end