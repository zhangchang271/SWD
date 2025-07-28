
function [res,cr_pre] = LHDispPick(ml,vmin,cr_0,method,saveForBackward)
if isfield(saveForBackward,'uxtposr')
        [~,ini] = max(max(saveForBackward.mlr));
        [~,iniP] = max(saveForBackward.mlr(:,ini));
else
        [~,ini] = max(max(saveForBackward.mll));
        [~,iniP] = max(saveForBackward.mll(:,ini));
end
    if method==1
        if ini==0
            
        end
        [idx,~,~]=pickFDCValue(ml,ini,iniP,2,4);
        cr_pre = idx' + vmin;
    else
        [~,p_max] = max(ml);
        cr_pre = p_max + vmin;

    end
    res = mean(((cr_pre(15:end-15,:) - cr_0(15:end-15,:)).^2));
    
end