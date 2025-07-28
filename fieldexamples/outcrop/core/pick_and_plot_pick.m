function [line_start_stop,line_part] = pick_and_plot_pick(id2,p_max,size_b,ml)
if id2(1)~=0
    id2 = [0,id2,size_b];
end
line_part1=[0,0];
for num_line=2:length(id2)
l2=p_max(id2(num_line-1)+1:id2(num_line));
if sum(ismember(line_part1,[l2(1),id2(num_line-1)+1],'rows'))==1
    continue
end
if id2(num_line)~=size_b
    for k =id2(num_line)+1:size_b
        [~,pfpk] = findpeaks(ml(:,k));
        if isempty(pfpk)
             [~,pfpk] = max(ml(:,k));
        end
        [~,minp] = min(abs(pfpk-l2(end)));
        %%%Validation of value
        [~,pfpk_val] = findpeaks(ml(:,k-1));
        if isempty(pfpk_val)
            [~,pfpk_val] = max(ml(:,k-1));
        end
        [~,minp_val] = min(abs(pfpk_val-pfpk(minp)));
        if pfpk_val(minp_val)==l2(end) 
            l2(end+1) = pfpk(minp);
        else
    %         start_interp=1;
    %         abandon_range=[abandon_range;k,id(j+1)+1];%%提取区间有问题，直接舍弃
            break
        end

    end
end
if id2(num_line-1)~=0
    for k =id2(num_line-1):-1:1
        [~,pfpk] = findpeaks(ml(:,k));
        if isempty(pfpk)
        [~,pfpk] = max(ml(:,k));
        end
        [~,minp] = min(abs(pfpk-l2(1)));
        %%%Validation of value
        [~,pfpk_val] = findpeaks(ml(:,k+1));
        if isempty(pfpk_val)
            [~,pfpk_val] = max(ml(:,k+1));
        end
        [~,minp_val] = min(abs(pfpk_val-pfpk(minp)));

        if pfpk_val(minp_val)==l2(1)% && pfpk(minp) - l2(end)<20
            l2 = [pfpk(minp), l2];
        else
    %         start_interp=1;
    %         abandon_range=[abandon_range;k,id(j+1)+1];%%提取区间有问题，直接舍弃
            break
        end
    end
    k=k-1;
end
ltfl2 = floor(length(l2)/5);
if  length(l2)<(size_b/5)% || sum(diff(l2(1:ltfl2)))/ltfl2>0.5 
    continue
end
if num_line==2
%     plot(l2,'Color','w','LineStyle','--','LineWidth',2)
    line_start_stop = [1,length(l2)];
    line_part{num_line-1}=[l2;1:length(l2)];
    line_part1=[line_part1;line_part{num_line-1}'];
elseif exist('line_part','var')==0
%     plot(l2,'Color','w','LineStyle','--','LineWidth',2)
    line_start_stop = [k+1,k+length(l2)];
    line_part{num_line-1}=[l2;(k+1):(k+length(l2))];
    line_part1=[line_part1;line_part{num_line-1}'];
else
%     plot((k+1):(k+length(l2)),l2,'Color','w','LineStyle','--','LineWidth',2)
    line_start_stop = [line_start_stop;k+1,k+length(l2)];
    line_part{num_line-1}=[l2;(k+1):(k+length(l2))];
end


end
if exist('line_part','var')==0
        line_start_stop = [1,size_b];
    line_part = {[p_max;1:size_b]};

end
end


% for k =id2(2):-1:1
%     [~,pfpk] = findpeaks(ml(:,k));
%     [~,minp] = min(abs(pfpk-l3(1)));
%     %%%Validation of value
%     [~,pfpk_val] = findpeaks(ml(:,k+1));
%     [~,minp_val] = min(abs(pfpk_val-pfpk(minp)));
%     if pfpk_val(minp_val)==l3(1)
%         l3 = [pfpk(minp), l3];
%     else
% %         start_interp=1;
% %         abandon_range=[abandon_range;k,id(j+1)+1];%%提取区间有问题，直接舍弃
%         break
%     end
% 
% end