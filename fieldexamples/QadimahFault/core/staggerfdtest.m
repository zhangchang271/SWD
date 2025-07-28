
function [seismo_u,seismo_w,seismo_ux,wavefield_gradient]=staggerfdtest(is,nbc,nt,dtx,dx,dt,sx,sz,gx,gz,s,vp,vs,isfs,fsz,fd_order,source_type,parameter_type,dt_wf,nt_wf)

% INPUT
% is:shot number
% sx,sz:shot coordinate
% gx,gx:geophone coordinate
% nbc: padding number (used in pad and ABC)
% nt: time length
% dtx: dt/dx,dt is the time interval, dx is the space interval
% s:source wavelet
% vp: p wave velocity
% vs: s wave velocity
% den: density
% isfs: 0:no free surface; 1: free surface
% fsz: free surface layer
% fd_order:finite difference accuracy order (24,26,28)
% source_type: different ways to add sources:p/s/w/tau_zz
% p: Initiate Strong P wave, weak S wave
% tau_zz/w: Both P and S wave
% s: Strong S wave, weak P wave
% When using free surface boundary conditio, source_type=p is
% recommended.

% OUTPUT
% seismo_u: horizontal displacement component recorded data
% seismo_w: vertical displacement component recorded data

%add wavefield recorded by Zongcai Feng  'wavefield'  2015-03-24
%add parameter_type by Zongcai Feng 'parameter_type'2015-04-07
%if parameter_type=0,input velocity; if parameter_type=1,lam and mu;

% Calculate lambda, mu based on density, p/s wave veolcity
% ca: lambda+2*mu; cl:lambda; cm: mu
%[ca,cm,cl]=calparam(vp,vs,den);

if (nargin)<=18
    dt_wf=dt;
    nt_wf=nt;
end
[nz,nx]=size(vp);
den=ones(nz,nx);
if parameter_type==0
    [ca,cm,cl]=calparam(vp,vs,den);
end
if parameter_type==1
    ca=cl+2*cm;
end

[nz,nx]=size(vp);
%if (nargout)==3
    wavefield_gradient.fux=zeros(nz,nx,nt_wf);
    wavefield_gradient.fuz=zeros(nz,nx,nt_wf);
    wavefield_gradient.bwx=zeros(nz,nx,nt_wf);
    wavefield_gradient.bwz=zeros(nz,nx,nt_wf);
%end

% staggered grid finite difference coeffcients
S21=1.0;
S41=9.0/8.0;S42=-1.0/24.0;
S61=1.17187;S62=-6.51042e-2;S63=4.68750e-3;
S81=1.19629;S82=-7.97526e-2;S83=9.57031e-3;S84=-6.97545e-4;

tic;

ng=numel(gx);
% pad means to expand the model space in order to add the absorbing
% boundary condition
if (isfs)
    pad_top=(fd_order-20)/2+1;
else
    pad_top=nbc;
end

cm=pad(cm,nbc,isfs,pad_top);
cl=pad(cl,nbc,isfs,pad_top);
ca=pad(ca,nbc,isfs,pad_top);
den=pad(den,nbc,isfs,pad_top);

% Adjust source and receiver position because of free surface
% if (isfs)
%     if (sz==fsz)
%         sz=sz;
%     end
%     if (gz==fsz)
%         gz=gz;
%     end
% end

% change source/geophone position because of pad
sx=sx+nbc;sz=sz+pad_top;
gx=gx+nbc;gz=gz+pad_top;

% if (isfs)
    fsz=fsz+pad_top;
% end

[nzbc,nxbc]=size(cm);

% allocate the size of the output
seismo_u=zeros(nt,ng);
seismo_w=zeros(nt,ng);
seismo_ux=zeros(nt,ng);

uu=zeros(nzbc,nxbc); % horizontal displacement wavefield
ww=zeros(nzbc,nxbc); % vertical displacement wavefield
xx=zeros(nzbc,nxbc); % tau_xx wavefield
zz=zeros(nzbc,nxbc); % tau_zz wavefield
xz=zeros(nzbc,nxbc); % tau_xz wavefield

fux=zeros(nzbc,nxbc);
fuz=zeros(nzbc,nxbc);
bwx=zeros(nzbc,nxbc);
bwz=zeros(nzbc,nxbc);
fux1=zeros(nzbc,nxbc);
fuz1=zeros(nzbc,nxbc);
bwx1=zeros(nzbc,nxbc);
bwz1=zeros(nzbc,nxbc);
% calculate dt/dx/dens
b=dtx./den;
% damp is used in the absorbing boundary condition
vmin=min(min(vp(:),vs(:)));
damp=damp_circle(vmin,nzbc,nxbc,nbc,dx,isfs,pad_top);
temp=1-damp*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Such scheme is abandoned.
% Adjust material parameter according to
% 'Free-surface boundary conditions for elastic staggered-grid modeling
% schemes' by Rune Mittet
if (isfs)
   cm(fsz,:)=1.0*cm(fsz,:);
   ca(fsz,:)=2*cm(fsz,:);
   cl(fsz,:)=0.0;
   b(fsz,:)=2*b(fsz,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % k=(pad_top:nzbc-pad_top);i=(pad_top:nxbc-pad_top);
% if (fd_order==22)
     k=(2:nzbc-1);i=(2:nxbc-1);
% elseif (fd_order==24)
%     k=(3:nzbc-2);i=(3:nxbc-2);
% elseif (fd_order==26)
%     k=(4:nzbc-3);i=(4:nxbc-3);
% elseif (fd_order==28)
%     k=(5:nzbc-4);i=(5:nxbc-4);
% end




for it=1:nt
    
    %display(['Modeling it=',num2str(it),' nt=',num2str(nt)]);
    % add source
  
%  if (strcmp(source_type,'p'))
%  
%         zz(sz,sx)=s(it);
%     elseif (strcmp(source_type,'s'))
%          uu(sz,sx)=uu(sz,sx)+s(it);
% %         ww(sz,sx)=s(it);
% %         %  ww(sz,sx)=ww(sz,sx)-s(it);% modified by zongcai
% %         ww(sz,sx+1)=ww(sz,sx+1)+s(it);
%     elseif (strcmp(source_type,'s11'))
%          uu(sz,sx)=uu(sz,sx)+s(it);
%              uu(sz+1,sx)=uu(sz+1,sx)-s(it);
%           ww(sz,sx)=ww(sz,sx)+s(it);
%             ww(sz,sx+1)=ww(sz,sx+1)+s(it);
%       elseif (strcmp(source_type,'s1'))
%           uu(sz,sx)=uu(sz,sx)+s(it);
%           %   uu(sz+1,sx)=uu(sz+1,sx)-s(it);
%           ww(sz,sx)=ww(sz,sx)+s(it);
%            % ww(sz,sx+1)=ww(sz,sx+1)+s(it);
%     elseif (strcmp(source_type,'tau_zz'))
%         zz(sz,sx)=zz(sz,sx)+s(it);
%     elseif (strcmp(source_type,'w'))
%         ww(sz,sx)=ww(sz,sx)+s(it);
%         %  ww(sz,sx)=s(it);
% %         ww(sz,sx+1)=ww(sz,sx+1)+s(it);
% %         ww(sz,sx-1)=ww(sz,sx-1)+s(it);
%         %     elseif (strcmp(source_type,'-p'))
%         %         xx(sz,sx)=xx(sz,sx)+s(it);
%         %         zz(sz,sx)=zz(sz,sx)-s(it);
%         %    elseif (strcmp(source_type,'tau_xz'))
%         %         xz(sz,sx)=xz(sz,sx)+s(it);
%     end

    
    % Calculate particle velocity
    if (fd_order==22)
     %  zz(fsz,:)=0.0;
        uu(k,i)=temp(k,i).*uu(k,i)+b(k,i).*(S21*(xx(k,i+1)-xx(k,i)+...
            xz(k,i)-xz(k-1,i)));
        ww(k,i)=temp(k,i).*ww(k,i)+b(k,i).*(S21*(xz(k,i)-xz(k,i-1)+...
            zz(k+1,i)-zz(k,i)));
    elseif (fd_order==24)
        %zz(fsz,:)=0.0;
        uu(k,i)=temp(k,i).*uu(k,i)+b(k,i).*(S41*(xx(k,i+1)-xx(k,i))+S42*(xx(k,i+2)-xx(k,i-1))+...
            S41*(xz(k,i)-xz(k-1,i))+S42*(xz(k+1,i)-xz(k-2,i)));
        
        ww(k,i)=temp(k,i).*ww(k,i)+b(k,i).*(S41*(xz(k,i)-xz(k,i-1))+S42*(xz(k,i+1)-xz(k,i-2))+...
            S41*(zz(k+1,i)-zz(k,i))+S42*(zz(k+2,i)-zz(k-1,i)));
    elseif (fd_order==26)
        uu(k,i)=temp(k,i).*uu(k,i)+b(k,i).*(S61*(xx(k,i)-xx(k,i-1))+S62*(xx(k,i+1)-xx(k,i-2))+...
            +S63*(xx(k,i+2)-xx(k,i-3))+S61*(xz(k,i)-xz(k-1,i))+S62*...
            (xz(k+1,i)-xz(k-2,i))+S63*(xz(k+2,i)-xz(k-3,i)));
        
        ww(k,i)=temp(k,i).*ww(k,i)+b(k,i).*(S61*(xz(k,i+1)-xz(k,i))+S62*(xz(k,i+2)-xz(k,i-1))+...
            +S63*(xz(k,i+3)-xz(k,i-2))+S61*(zz(k+1,i)-zz(k,i))+S62*...
            (zz(k+2,i)-zz(k-1,i))+S63*(zz(k+3,i)-zz(k-2,i)));
    elseif (fd_order==28)
        uu(k,i)=temp(k,i).*uu(k,i)+b(k,i).*(S81*(xx(k,i)-xx(k,i-1))+S82*(xx(k,i+1)-xx(k,i-2))+...
            +S83*(xx(k,i+2)-xx(k,i-3))+S84*(xx(k,i+3)-xx(k,i-4))+...
            S81*(xz(k,i)-xz(k-1,i))+S82*(xz(k+1,i)-xz(k-2,i))+S83*(xz(k+2,i)-xz(k-3,i))+...
            S84*(xz(k+3,i)-xz(k-4,i)));
        
        ww(k,i)=temp(k,i).*ww(k,i)+b(k,i).*(S81*(xz(k,i+1)-xz(k,i))+S82*(xz(k,i+2)-xz(k,i-1))+...
            S83*(xz(k,i+3)-xz(k,i-2))+S84*(xz(k,i+4)-xz(k,i-3))+...
            S81*(zz(k+1,i)-zz(k,i))+S82*(zz(k+2,i)-zz(k-1,i))+S83*(zz(k+3,i)-zz(k-2,i))+...
            S84*(zz(k+4,i)-zz(k-3,i)));
    end
    
    % Free surface boundary condition for velocity
%    if (isfs)
%      ww(fsz-1,i)=ww(fsz,i)+cl(fsz,i)/ca(fsz,i)*(uu(fsz,i+1)-uu(fsz,i));
% %     
% %      xx(fsz-1,i)=xx(fsz,i)+(xx(fsz+1,i)-xx(fsz,i)+ww(fsz,i+1)-ww(fsz,i))+(ww(fsz-1,i+1)-ww(fsz-1,i));
%   end
 
    
if (strcmp(source_type,'p'))
          xx(sz,sx)=xx(sz,sx)+s(it);
           xz(sz,sx)=xz(sz,sx)+s(it);
        zz(sz,sx)=zz(sz,sx)+s(it);
    elseif (strcmp(source_type,'s'))
         uu(sz,sx)=uu(sz,sx)+s(it);
%         ww(sz,sx)=s(it);
%         %  ww(sz,sx)=ww(sz,sx)-s(it);% modified by zongcai
%         ww(sz,sx+1)=ww(sz,sx+1)+s(it);
    elseif (strcmp(source_type,'s11'))
         uu(sz,sx)=uu(sz,sx)+s(it);
             uu(sz+1,sx)=uu(sz+1,sx)-s(it);
          ww(sz,sx)=ww(sz,sx)+s(it);
            ww(sz,sx+1)=ww(sz,sx+1)+s(it);
      elseif (strcmp(source_type,'s1'))
          uu(sz,sx)=uu(sz,sx)+s(it);
          %   uu(sz+1,sx)=uu(sz+1,sx)-s(it);
          ww(sz,sx)=ww(sz,sx)+s(it);
           % ww(sz,sx+1)=ww(sz,sx+1)+s(it);
    elseif (strcmp(source_type,'tau_zz'))
        zz(sz,sx)=zz(sz,sx)+s(it);
    elseif (strcmp(source_type,'w'))
        ww(sz,sx)=ww(sz,sx)+s(it);
        %  ww(sz,sx)=s(it);
%         ww(sz,sx+1)=ww(sz,sx+1)+s(it);
%         ww(sz,sx-1)=ww(sz,sx-1)+s(it);
        %     elseif (strcmp(source_type,'-p'))
        %         xx(sz,sx)=xx(sz,sx)+s(it);
        %         zz(sz,sx)=zz(sz,sx)-s(it);
        %    elseif (strcmp(source_type,'tau_xz'))
        %         xz(sz,sx)=xz(sz,sx)+s(it);
    end

    % update stress
    if  (fd_order==22)
        fux(k,i)=S21*(uu(k,i)-uu(k,i-1));
        fuz(k,i)=S21*(uu(k+1,i)-uu(k,i));
        bwx(k,i)=S21*(ww(k,i+1)-ww(k,i));
        bwz(k,i)=S21*(ww(k,i)-ww(k-1,i));
%         fux1(k,i)=S21*(fux(k,i)-fux(k,i-1));
%         fuz1(k,i)=S21*(fuz(k,i)-fuz(k,i-1));
%         bwx1(k,i)=S21*(bwx(k,i)-bwx(k,i-1));
%         bwz1(k,i)=S21*(bwz(k,i)-bwz(k,i-1));
    elseif (fd_order==24)
        fux(k,i)=S41*(uu(k,i)-uu(k,i-1))+S42*(uu(k,i+1)-uu(k,i-2));
        fuz(k,i)=S41*(uu(k+1,i)-uu(k,i))+S42*(uu(k+2,i)-uu(k-1,i));
        bwx(k,i)=S41*(ww(k,i+1)-ww(k,i))+S42*(ww(k,i+2)-ww(k,i-1));
        bwz(k,i)=S41*(ww(k,i)-ww(k-1,i))+S42*(ww(k+1,i)-ww(k-2,i));
    elseif (fd_order==26)
        fux(k,i)=S61*(uu(k,i+1)-uu(k,i))+S62*(uu(k,i+2)-uu(k,i-1))+...
            S63*(uu(k,i+3)-uu(k,i-2));
        fuz(k,i)=S61*(uu(k+1,i)-uu(k,i))+S62*(uu(k+2,i)-uu(k-1,i))+...
            S63*(uu(k+3,i)-uu(k-2,i));
        bwx(k,i)=S61*(ww(k,i)-ww(k,i-1))+S62*(ww(k,i+1)-ww(k,i-2))+...
            S63*(ww(k,i+2)-ww(k,i-3));
        bwz(k,i)=S61*(ww(k,i)-ww(k-1,i))+S62*(ww(k+1,i)-ww(k-2,i))+...
            S63*(ww(k+2,i)-ww(k-3,i));
    elseif (fd_order==28)
        fux(k,i)=S81*(uu(k,i+1)-uu(k,i))+S82*(uu(k,i+2)-uu(k,i-1))+...
            S83*(uu(k,i+3)-uu(k,i-2))+S84*(uu(k,i+4)-uu(k,i-3));
        fuz(k,i)=S81*(uu(k+1,i)-uu(k,i))+S82*(uu(k+2,i)-uu(k-1,i))+...
            S83*(uu(k+3,i)-uu(k-2,i))+S84*(uu(k+4,i)-uu(k-3,i));
        bwx(k,i)=S81*(ww(k,i)-ww(k,i-1))+S82*(ww(k,i+1)-ww(k,i-2))+...
            S83*(ww(k,i+2)-ww(k,i-3))+S84*(ww(k,i+3)-ww(k,i-4));
        bwz(k,i)=S81*(ww(k,i)-ww(k-1,i))+S82*(ww(k+1,i)-ww(k-2,i))+...
            S83*(ww(k+2,i)-ww(k-3,i))+S84*(ww(k+3,i)-ww(k-4,i));
    end
    xx=temp.*xx+(ca.*fux+cl.*bwz)*dtx;
    zz=temp.*zz+(ca.*bwz+cl.*fux)*dtx;
    xz=temp.*xz+cm.*(fuz+bwx)*dtx;
    
%     % free surface boundary condition
%            if (isfs)
         zz(fsz,:)=0.0;
%     end
    
    % record data
    for ig=1:ng
        seismo_w(it,ig)=ww(gz(ig),gx(ig));
        seismo_u(it,ig)=uu(gz(ig),gx(ig));
        seismo_ux(it,ig)= fux(gz(ig),gx(ig))/dx;
    end
    
    if (nargout)==4
        in_wf=dt_wf/dt;
        if mod(it-1,in_wf)==0
            it_record=1+floor((it-1)/in_wf);
            wavefield_gradient.fux(:,:,it_record)=fux(pad_top+1:end-nbc,nbc+1:end-nbc)/dx;
            wavefield_gradient.fuz(:,:,it_record)=fuz(pad_top+1:end-nbc,nbc+1:end-nbc)/dx;
            wavefield_gradient.bwx(:,:,it_record)=bwx(pad_top+1:end-nbc,nbc+1:end-nbc)/dx;
            wavefield_gradient.bwz(:,:,it_record)=bwz(pad_top+1:end-nbc,nbc+1:end-nbc)/dx;
%             wavefield_gradient.fux(:,1:end-1,it_record)=(fux(pad_top+1:end-nbc,nbc+2:end-nbc)/dx-fux(pad_top+1:end-nbc,nbc+1:end-nbc-1)/dx)/dx;
%             wavefield_gradient.fuz(:,1:end-1,it_record)=(fuz(pad_top+1:end-nbc,nbc+2:end-nbc)/dx-fuz(pad_top+1:end-nbc,nbc+1:end-nbc-1)/dx)/dx;
%             wavefield_gradient.bwx(:,1:end-1,it_record)=(bwx(pad_top+1:end-nbc,nbc+2:end-nbc)/dx-bwx(pad_top+1:end-nbc,nbc+1:end-nbc-1)/dx)/dx;
%             wavefield_gradient.bwz(:,1:end-1,it_record)=(bwz(pad_top+1:end-nbc,nbc+2:end-nbc)/dx-bwz(pad_top+1:end-nbc,nbc+1:end-nbc-1)/dx)/dx;
% wavefield_gradient.fux(:,end,it_record)=wavefield_gradient.fux(:,end-1,it_record);
% wavefield_gradient.fuz(:,end,it_record)=wavefield_gradient.fuz(:,end-1,it_record);
% wavefield_gradient.bwx(:,end,it_record)=wavefield_gradient.bwx(:,end-1,it_record);
% wavefield_gradient.bwz(:,end,it_record)=wavefield_gradient.bwz(:,end-1,it_record);
            %             wavefield_gradient.uu(:,:,it_record)=uu(pad_top+1:end-nbc,nbc+1:end-nbc)/dx;
%             wavefield_gradient.ww(:,:,it_record)=ww(pad_top+1:end-nbc,nbc+1:end-nbc)/dx;
        end
    end
    
    %     % plot the snapshot of the wavefield
    %     if it/30 == round(it/30);
    %         plotit(dx,ww,nxbc,nzbc,it,dt);
    %         %   figure(8);imagesc(wavefield_gradient.bwx(:,:,it));
    %     end
    %
    
end

%display(num2str(is),'th shot');
end
