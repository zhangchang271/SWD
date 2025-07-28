function [id2,line_start_stop,line_part]=pickdisp_xl(ml)
[size_a,size_b] = size(ml);
% imagesc(ml(:,:))
% hold on
[~,p_max] = max(ml(:,:));
% plot(p_max,'b')
dline = abs(diff(p_max));
dline(dline<10)=0;
[~,id2] = findpeaks(dline);


% if ~isempty(id2)
%     if p_max(id2(1))-p_max(id2(1)+1)>0
%         id2(1)=[];
%     end
% end
if ~isempty(id2)
    [line_start_stop,line_part] = pick_and_plot_pick(id2,p_max,size_b,ml);
%     [line_start_stop,line_part] = CEDC(ml);
else
    line_start_stop = [1,size_b];
    line_part = {[p_max;1:size_b]};
end


line_part(cellfun(@isempty,line_part))=[];
% [line_start_stop,ia,~] = unique(line_start_stop,'rows');
% line_part = line_part(ia);

%%%%%%%%%%%%%%%%%%%%%eneragy check
energavgia=[];
for Tmode=1:length(line_part)
    D = diag(ml(line_part{Tmode}(1,:),line_part{Tmode}(2,:)));
    D1=D;
    D1(D>0.65)=1;
    D1(D<0.65)=0;
    energavg = sum(D1)/(line_start_stop(Tmode,2)-line_start_stop(Tmode,1)+1);
    if energavg > 0.15
        energavgia=[energavgia;Tmode];
    end
end
line_part = line_part(energavgia);
line_start_stop = line_start_stop(energavgia,:);



fundamental_start=p_max(1);
endlen=floor(size_b/5);
CenterPoint=zeros(length(line_part),4);

for Tmode=1:length(line_part)
    lla=(sum(line_part{Tmode},2))/(line_start_stop(Tmode,2)-line_start_stop(Tmode,1)+1);
    CenterPoint(Tmode,1:2)=lla';
    CenterPoint(Tmode,3)=sqrt((lla(1)/size_a)^2+(lla(2)/size_b)^2);
    if fundamental_start == line_part{Tmode}(1,1) && line_part{Tmode}(2,1) < 5
        CenterPoint(Tmode,4)=1;
        fundamental_Tmode = Tmode;
    end
end
if exist('fundamental_Tmode','var')==1
    [~,index] = sortrows(CenterPoint,1);
    index2 = index(find(index==fundamental_Tmode):end);
    CenterPoint = CenterPoint(index2,:);
    line_part = line_part(index2);
    line_start_stop = line_start_stop(index2,:);
end

[~,index3] = sortrows(CenterPoint,3);
line_part = line_part(index3);
line_start_stop = line_start_stop(index3,:);

%%%%%%%%%%%%%%partly identical
linesum=zeros(length(line_part),1);
for Tmode=1:length(line_part)
    linesum(Tmode,1)=sum(line_part{Tmode}(1,end-endlen:end));
end
[~,ib,~] = unique(linesum,'stable');
line_part = line_part(ib);
line_start_stop = line_start_stop(ib,:);



end


% energavgia=[];
% for Tmode=1:length(line_part)
%     D = diag(ml(line_part{Tmode}(1,:),line_part{Tmode}(2,:)));
%     D1=D;
%     D1(D>0.65)=1;
%     D1(D<0.65)=0;
%     energavg = sum(D1)/(line_start_stop(Tmode,2)-line_start_stop(Tmode,1)+1);
%     if energavg > 0.15
%         energavgia=[energavgia;Tmode];
%     end
% end
% line_part = line_part(energavgia);
% line_start_stop = line_start_stop(energavgia,:);

