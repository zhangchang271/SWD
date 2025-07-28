function [p_max]=pick_disp(ml)
% [~,size_b,size_c] = size(ml);
[~,p_max] = max(ml(:,:));
dline = abs(diff(p_max));
dline(dline<15)=0;
[~,id] = findpeaks(dline);
for j=1:length(id)-1
    if (id(j+1)-id(j))<20
        for k =id(j):id(j+1)
            [~,pfpk] = findpeaks(ml(:,k));
            [~,minp] = min(abs(pfpk-p_max(k-1)));
             if ~isempty(pfpk)
                p_max(k) = pfpk(minp);
             end
        end
    end
end
end