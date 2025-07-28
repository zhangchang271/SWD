
function [seismo_u,seismo_w,www]=staggerfd8(is,nbc,nt,dtx,dx,dt,sx,sz,gx,gz,s,vp,vs,isfs,fsz,fd_order,source_type,parameter_type,dt_wf,nt_wf)

% INPUT
% is:shot number
% sx,sz:shot coordinate
% gx,gx:geophone coordinate
% nbc: padding number (used in pad and ABC)
% nt: time length
% dtx: dt/dx,dt is the time interval, dx is the space interval
% s:source waveletfsz
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
[nz,nx]=size(vs);
den=single(ones(nz,nx));
if parameter_type==0
    [ca,cm,cl]=calparam(vp,vs,den);
% [cl,cm] = calepara(vs,den,vp);
end
if parameter_type==1
    ca=cl+2*cm;
end

[nz,nx]=size(vp);
%if (nargout)==3
%     wavefield_gradient.fux=zeros(nz,nx,nt_wf);
%     wavefield_gradient.fuz=zeros(nz,nx,nt_wf);
%     wavefield_gradient.bwx=zeros(nz,nx,nt_wf);
%     wavefield_gradient.bwz=zeros(nz,nx,nt_wf);
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

den=pad(den,nbc,isfs,pad_top);

% Adjust source and receiver position because of free surface


% change source/geophone position because of pad
sx=sx+nbc;sz=sz+pad_top;
gx=gx+nbc;gz=gz+pad_top;


[nzbc,nxbc]=size(cm);
nzbc=single(nzbc);nxbc=single(nxbc);
% allocate the size of the output
seismo_u=zeros(nt,ng);
seismo_w=zeros(nt,ng);

uu=zeros(nzbc,nxbc); % horizontal displacement wavefield
ww=zeros(nzbc,nxbc); % vertical displacement wavefield
xx=zeros(nzbc,nxbc); % tau_xx wavefield
zz=zeros(nzbc,nxbc); % tau_zz wavefield
xz=zeros(nzbc,nxbc); % tau_xz wavefield

fux=zeros(nzbc,nxbc);
fuz=zeros(nzbc,nxbc);
bwx=zeros(nzbc,nxbc);
bwz=zeros(nzbc,nxbc);

% calculate dt/dx/dens

% damp is used in the absorbing boundary condition
vmin=min(min(vp(:),vs(:)));
damp=damp_circle(vmin,nzbc,nxbc,nbc,dx,isfs,pad_top);
temp=single(1-damp*dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Such scheme is abandoned.
% Adjust material parameter according to
% 'Free-surface boundary conditions for elastic staggered-grid modeling
% schemes' by Rune Mittet
if (isfs)
    ca=pad(ca,nbc,isfs,pad_top);
   cm(pad_top,:)=1.0*cm(pad_top,:);
   ca(pad_top,:)=2*cm(pad_top,:);
   cl(pad_top,:)=0.0;
   b=dtx./den;
   b(pad_top,:)=2*b(pad_top,:);
      b1=b;
   cm1=cm;
end


% if (isfs)
%    den(1:pad_top,:) = 0.5*den(1:pad_top,:);
%    den1=den;
%    den1(pad_top,:)=2*den1(pad_top,:);
%    cm(1:pad_top,:)=0.5*cm(1:pad_top,:);
%    cm1=cm;
%    cl(1:pad_top,:)=0.0;
%    cm1(pad_top,:)=2*cm1(pad_top,:);
%    ca=cl + 2*cm;
%    b=dtx./den;
%    b1=b;
%    b1=dtx./den1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input_vector = [nt,nzbc,nxbc,dtx,ng,source_type,fd_order];
    if (strcmp(source_type,'p'))
        source_type_num=1;
    elseif (strcmp(source_type,'s'))
        source_type_num=2;
    elseif (strcmp(source_type,'s1'))
        source_type_num=3;
    elseif (strcmp(source_type,'z'))
        source_type_num=4;
    elseif (strcmp(source_type,'w'))
        source_type_num=5;
    end
fd_order_num=fd_order;in_wf=dt_wf/dt;format_num=3;%format_num: dim of wavefield.fux...
input_vector = [nt,nzbc,nxbc,dtx,ng,sz,sx,gz(1),gx(1),gx(2)-gx(1),source_type_num,fd_order_num,in_wf,nz,nx,format_num];
[seismo_u,seismo_w,www]= eigen_staggerfd_single(input_vector,temp,ca,cl,cm,cm1,b,b1,s);
% res2=seismo_w;
%display(num2str(is),'th shot');
end
