%
 function [ml,dataLen]=RTl(seismo_v1,is,df,dt,np,vmin,vmax,fmin,fmax,a,b,ds,offset,m,FK,M)
%#######################################################################!
%if is<=b-offset
  nend=M*is-M+1-offset-m;
  if nend<1
      nend=0;
  end
  uxt=fliplr(seismo_v1(:,nend+1:round((M*is-M+1-m))));
  [~,dataLen] = size(uxt);
  offset=M*is-M+1-m-nend;
  x=(m+linspace(0,offset-1,offset))*ds;
  
 %uxt=seismo_v1(:,round((M*is-M+1+m)):offset+m+(M*is-M+1));
 if FK==1% normalization
  [uxt] = fk_filter(uxt, dt, ds);
 end
%% 
ccn=fix(1./df./dt);
d=fft(uxt,ccn);
d=d.';
lf=round(fmin./df)+1;
nf=round(fmax./df)+1;
pp=1./(vmin:vmax)';
ll0=1i.*2.*pi.*df.*(pp*x);
mm=zeros(np,nf);
for luoj=lf:nf
l=exp(ll0.*(luoj-1));
mm(:,luoj)=(l*d(:,luoj));
end
ml=abs(mm);
for luoi=lf:nf
 ml(:,luoi)=ml(:,luoi)./max(ml(:,luoi));
end
 ml=ml(:,lf:nf);
% [maxpowr,idx]=max(ml(:,:));
