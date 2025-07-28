function [cr_0,dvmin,dvmax]=pickFDC(ml,lp,fg,lfg)
[a,b]=size(ml);
idx=zeros(b,1);
dvmax=zeros(b,1);
dvmin=zeros(b,1);
qujian=nan(b,1);
chazhi=ones(b,1);
if nargin < 2
   fg=1.5; %高阶与基阶交汇时energy区间变大的倍数
   lp= 10;%起始点（起始点前的曲线pick max energy）
   lfg=fg;
elseif nargin < 3
   fg=1.5; %高阶与基阶交汇时energy区间变大的倍数
   lfg=fg;
elseif nargin < 4
lfg=fg; %高阶与基阶交汇时energy区间变大的倍数left
end
% [~,idx(1:lp)]=max(ml(:,1:lp)); 
[~,idx(lp)]=max(ml(:,lp));
%% left of lp
for i=lp-1:-1:1
    if idx(i+1)==a
        idx(i)=a;
        dvmax(i)=dvmax(i+1);
        dvmin(i)=dvmin(i+1); 
    else
    if idx(i+1)>a-3
     dvmax(i)=a;   
    else
    [~,mlmax]=findpeaks(-ml(idx(i+1):a,i));
       if isempty(mlmax)
         dvmax(i)=a;
       else
         dvmax(i)=idx(i+1)+min(mlmax); 
       end
    end
    if idx(i+1)<3
     dvmin(i)=1;  
    else
     [~,mlmin]=findpeaks(-ml(1:idx(i+1),i));
       if isempty(mlmin)
        dvmin(i)=1;
       else
        dvmin(i)=max(mlmin);   
       end    
    end
    qujian(i)=dvmax(i)-dvmin(i);
    if qujian(i)>qujian(i+1)*lfg
%        if all(diff(ml(dvmin(i-1)+1:dvmax(i-1),i))>=0)|| all(diff(ml(dvmin(i-1)+1:dvmax(i-1),i))<=0)
        
         idx(i)=idx(i+1);
         chazhi(i)=0;
         dvmax(i)=dvmax(i+1);
         dvmin(i)=dvmin(i+1); 
         qujian(i)=dvmax(i)-dvmin(i);
%        else
%          [~,mlmax1]=findpeaks(ml(dvmin(i-1):dvmax(i-1),i));
%          if isempty(mlmax1)
%          idx(i)=idx(i-1);
%          chazhi(i)=0;
%           dvmax(i)=dvmax(i-1);
%           dvmin(i)=dvmin(i-1);
%           qujian(i)=dvmax(i)-dvmin(i);
%          else
%           dvmax(i)=dvmax(i-1);
%         [~,jin]=min(abs(mlmax1-idx(i-1)));
%           dvmin(i)=dvmin(i-1);
%          idx(i)=dvmin(i-1)+mlmax1(jin);
%         qujian(i)=dvmax(i)-dvmin(i);
%          end
%        end
    else
         [~,id]=max(ml(dvmin(i):dvmax(i),i));
         idx(i)=dvmin(i)+id-1;
    end
    end
end
%% right of lp
 for i=lp+1:b
    if idx(i-1)>a-3
     dvmax(i)=a;   
    else
    [~,mlmax]=findpeaks(-ml(idx(i-1):a,i));
       if isempty(mlmax)
         dvmax(i)=a;
       else
         dvmax(i)=idx(i-1)+min(mlmax); 
       end
    end
    if idx(i-1)<3
     dvmin(i)=1;  
    else
     [~,mlmin]=findpeaks(-ml(1:idx(i-1),i));
       if isempty(mlmin)
        dvmin(i)=1;
       else
        dvmin(i)=max(mlmin);   
       end    
    end
    qujian(i)=dvmax(i)-dvmin(i);
    if qujian(i)>fg*qujian(i-1) 
%        if all(diff(ml(dvmin(i-1)+1:dvmax(i-1),i))>=0)|| all(diff(ml(dvmin(i-1)+1:dvmax(i-1),i))<=0)
         idx(i)=idx(i-1);
         chazhi(i)=0;
         dvmax(i)=dvmax(i-1);
         dvmin(i)=dvmin(i-1); 
         qujian(i)=dvmax(i)-dvmin(i);
%        else
%          [~,mlmax1]=findpeaks(ml(dvmin(i-1):dvmax(i-1),i));
%          if isempty(mlmax1)
%          idx(i)=idx(i-1);
%          chazhi(i)=0;
%           dvmax(i)=dvmax(i-1);
%           dvmin(i)=dvmin(i-1);
%           qujian(i)=dvmax(i)-dvmin(i);
%          else
%           dvmax(i)=dvmax(i-1);
%         [~,jin]=min(abs(mlmax1-idx(i-1)));
%           dvmin(i)=dvmin(i-1);
%          idx(i)=dvmin(i-1)+mlmax1(jin);
%         qujian(i)=dvmax(i)-dvmin(i);
%          end
%        end
    else
         [~,id]=max(ml(dvmin(i):dvmax(i),i));
         idx(i)=dvmin(i)+id-1;
    end
end
que=find(chazhi==1);
idxc=interp1(que,idx(que),1:b,'linear');
cr_0=interp1(find(~isnan(idxc)),idxc(~isnan(idxc)),1:b,'nearest','extrap');

