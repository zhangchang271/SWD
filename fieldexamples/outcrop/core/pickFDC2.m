function [cr_0]=pickFDC2(ml,data,dt,ffreq,lp,case1)
[a,b]=size(ml);
% lineall=zeros(a,b);
idx=zeros(1,b);
 chazhi=ones(b,1);
 if nargin < 5
%      win=round(b/20);
%      for i=(1+win):round(2*b/(1+sqrt(5)))
%         eg(i)=sum(sum(ml(:,i-win:i+win)));
%      end
% [~,maxi]=max(ml(:,round(b/5):round(b*3/5)));
%      [~,lp1]= min(maxi);%起始点
%     lp=lp1+round(b/5)-1;
[nt,bb]=size(data);
for i=1:bb
    data(:,i)=data(:,i)/max(abs(data(:,i)));
end
f=fft(data);
amp=abs(f(1:floor(nt/2),:));
ff=1/dt/nt*(0:floor(nt/2)-1);
[~,inda]=max(mean(amp(:,:),2));
[~, lp]=min(abs(ffreq-ff(inda)));
    case1=2;
    elseif nargin < 6
    case1=2; 

 end
[~,idx(lp)]=max(ml(:,lp));
%% left of lp
if (case1 == 1)
[~,lomax1]=findpeaks(ml(:,lp));
[~,lomin1]=findpeaks(-ml(:,lp));
for ixml=lp-1:-1:1
[~,lomax2]=findpeaks(ml(:,ixml));
[~,lomin2]=findpeaks(-ml(:,ixml));
mmax1=length(lomax1);
nmin1=length(lomin1);
length1=mmax1+nmin1;
line1=zeros(length1,1);

if mmax1<1
    line1=lomin1;
elseif nmin1<1
    line1=lomax1;
else
if lomax1(1)>lomin1(1)
    a1=1;
    line1(1:2:length1)=lomin1;
    line1(2:2:length1)=lomax1;
else
    a1=0;
    line1(1:2:length1)=lomax1;
    line1(2:2:length1)=lomin1;
end
end
mmax2=length(lomax2);
nmin2=length(lomin2);
length2=mmax2+nmin2;
line2=zeros(length2,1);
if mmax2<1
    line2=lomin2;
elseif nmin2<1
    line2=lomax2;
else
if lomax2(1)>lomin2(1)
    a2=1;
    line2(1:2:length2)=lomin2;
    line2(2:2:length2)=lomax2;
else
    a2=0;
    line2(1:2:length2)=lomax2;
    line2(2:2:length2)=lomin2;
end
end
m=length(line1);
n=length(line2);
cost=nan(m,n);
if a1==a2
for i=1:m
    for j=1:n
        if mod((i-j),2)==1
        else
        cost(i,j)=abs(line1(i)-line2(j));
        end
    end
end
else
  for i=1:m
    for j=1:n
        if mod((i-j),2)==1
            cost(i,j)=abs(line1(i)-line2(j));
        end
    end
  end  
end
 lj1=zeros(1,m);
 lj2=zeros(1,n);
for i=1:m
[~,lj1(i)]=min(cost(i,:));    
end
for j=1:n
[~,lj2(j)]=min(cost(:,j));    
end
 ljj1=zeros(1,m);
 ljj2=zeros(1,n);
for i=1:m
if lj2(lj1(i))==i
    ljj1(i)=lj1(i);
end
end
for j=1:n
if lj1(lj2(j))==j
    ljj2(j)=lj2(j);
end
end

ljj1(ljj1==0)=[];
ljj2(ljj2==0)=[];
% dian1max=intersect(line1(ljj2),lomax1);
% dian1min=intersect(line1(ljj2),lomin1);
% dian2max=intersect(line2(ljj1),lomax2);
% dian2min=intersect(line2(ljj1),lomin2);

ind1=find(line1(ljj2)==idx(ixml+1));
if ~isempty(ind1) && idx(ixml+1)<a-10
idx(ixml)=line2(ljj1(ind1));
% lineall(dian1max,ixml-1)=1;
% lineall(dian1min,ixml-1)=-1;
% lineall(dian2max,ixml)=1;
% lineall(dian2min,ixml)=-1;

% lineall(line1(ljj2),ixml-1)=1;
% lineall(line2(ljj1),ixml)=1;
lomax1=lomax2;
lomin1=lomin2;
else
    idx(ixml)=idx(ixml+1);
    chazhi(i)=0;
end
end
else
    if idx(lp)>a-3
     dvmax(lp)=a;   
    else
    [~,mlmax]=findpeaks(-ml(idx(lp):a,lp));
       if isempty(mlmax)
         dvmax(lp)=a;
       else
         dvmax(lp)=idx(lp)+min(mlmax); 
       end
    end
    if idx(lp)<3
     dvmin(lp)=1;  
    else
     [~,mlmin]=findpeaks(-ml(1:idx(lp),lp));
       if isempty(mlmin)
        dvmin(lp)=1;
       else
        dvmin(lp)=max(mlmin);   
       end    
    end 
for ixml=lp-1:-1:1
    if idx(ixml+1)==a
        idx(ixml)=a;
        dvmax(ixml)=dvmax(ixml+1);
        dvmin(ixml)=dvmin(ixml+1); 
    else
    if idx(ixml+1)>a-3
     dvmax(ixml)=a;   
    else
    [~,mlmax]=findpeaks(-ml(idx(ixml+1):a,ixml));
       if isempty(mlmax)
         dvmax(ixml)=a;
       else
         dvmax(ixml)=idx(ixml+1)+min(mlmax); 
       end
    end
    if idx(ixml+1)<3
     dvmin(ixml)=1;  
    else
     [~,mlmin]=findpeaks(-ml(1:idx(ixml+1),ixml));
       if isempty(mlmin)
        dvmin(ixml)=1;
       else
        dvmin(ixml)=max(mlmin);   
       end    
    end
%     qujian(i)=dvmax(i)-dvmin(i);
%     if qujian(i)>qujian(i+1)*lfg
% %        if all(diff(ml(dvmin(i-1)+1:dvmax(i-1),i))>=0)|| all(diff(ml(dvmin(i-1)+1:dvmax(i-1),i))<=0)
%         
%          idx(i)=idx(i+1);
%          chazhi(i)=0;
%          dvmax(i)=dvmax(i+1);
%          dvmin(i)=dvmin(i+1); 
%          qujian(i)=dvmax(i)-dvmin(i);
%     else
         [~,id]=max(ml(dvmin(ixml):dvmax(ixml),ixml));
         idx(ixml)=dvmin(ixml)+id-1;
%     end
    end
end
end
%% right of lp
[~,lomax1]=findpeaks(ml(:,lp));
[~,lomin1]=findpeaks(-ml(:,lp));
for ixml=lp+1:b
[~,lomax2]=findpeaks(ml(:,ixml));
[~,lomin2]=findpeaks(-ml(:,ixml));
mmax1=length(lomax1);
nmin1=length(lomin1);
length1=mmax1+nmin1;
line1=zeros(length1,1);
if mmax1<1
    line1=lomin1;
elseif nmin1<1
    line1=lomax1;
else
if lomax1(1)>lomin1(1)
    a1=1;
    line1(1:2:length1)=lomin1;
    line1(2:2:length1)=lomax1;
else
    a1=0;
    line1(1:2:length1)=lomax1;
    line1(2:2:length1)=lomin1;
end
end
mmax2=length(lomax2);
nmin2=length(lomin2);
length2=mmax2+nmin2;
line2=zeros(length2,1);
if mmax2<1
    line2=lomin2;
elseif nmin2<1
    line2=lomax2;
else
if lomax2(1)>lomin2(1)
    a2=1;
    line2(1:2:length2)=lomin2;
    line2(2:2:length2)=lomax2;
else
    a2=0;
    line2(1:2:length2)=lomax2;
    line2(2:2:length2)=lomin2;
end
end
m=length(line1);
n=length(line2);
cost=nan(m,n);
if a1==a2
for i=1:m
    for j=1:n
        if mod((i-j),2)==1
        else
        cost(i,j)=abs(line1(i)-line2(j));
        end
    end
end
else
  for i=1:m
    for j=1:n
        if mod((i-j),2)==1
            cost(i,j)=abs(line1(i)-line2(j));
        end
    end
  end  
end
 lj1=zeros(1,m);
 lj2=zeros(1,n);
for i=1:m
[~,lj1(i)]=min(cost(i,:));    
end
for j=1:n
[~,lj2(j)]=min(cost(:,j));    
end
 ljj1=zeros(1,m);
 ljj2=zeros(1,n);
for i=1:m
if lj2(lj1(i))==i
    ljj1(i)=lj1(i);
end
end
for j=1:n
if lj1(lj2(j))==j
    ljj2(j)=lj2(j);
end
end

ljj1(ljj1==0)=[];
ljj2(ljj2==0)=[];
% dian1max=intersect(line1(ljj2),lomax1);
% dian1min=intersect(line1(ljj2),lomin1);
% dian2max=intersect(line2(ljj1),lomax2);
% dian2min=intersect(line2(ljj1),lomin2);

ind1=find(line1(ljj2)==idx(ixml-1));
if ind1
idx(ixml)=line2(ljj1(ind1));
% lineall(dian1max,ixml-1)=1;
% lineall(dian1min,ixml-1)=-1;
% lineall(dian2max,ixml)=1;
% lineall(dian2min,ixml)=-1;

% lineall(line1(ljj2),ixml-1)=1;
% lineall(line2(ljj1),ixml)=1;
lomax1=lomax2;
lomin1=lomin2;
else
    idx(ixml)=idx(ixml-1);
    chazhi(ixml)=0;
end
end
que=find(chazhi==1);
idxc=interp1(que,idx(que),1:b,'linear');
cr_0=interp1(find(~isnan(idxc)),idxc(~isnan(idxc)),1:b,'nearest','extrap');
end