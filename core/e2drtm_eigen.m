% 2D elastic modeling by staggered grid method
% @version 2 2014-09-26
% @author Bowen Guo

%Born modeling by Zongcai Feng 2015-03-24
%add parameter: vp_refl,vs_refl

function [cl_img,cm_img,illum_div]=e2drtm_eigen(wavefield_gradient,seismo_w,is,nbc,nt,dtx,dx,dt,gx,gz,s,vp,vs,isfs,fd_order,parameter_type,in_wf)



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
% vp_refl: P wave reflectivity, vp=vp0(1+vp_refl)
% vs_refl: s wave reflectivity, vs=vs0(1+vs_refl)
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


% Calculate lambda, mu based on density, p/s wave veolcity
% ca: lambda+2*mu; cl:lambda; cm: mu
%[ca,cm,cl]=calparam(vp,vs,den);

[nz,nx]=size(vp);
den=ones(nz,nx,'single');
if parameter_type==0
    [ca,cm,cl]=calparam(vp,vs,den);
end
if parameter_type==1
    ca=cl+2*cm;
end



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



% change source/geophone position because of pad
gx=gx+nbc;gz=gz+pad_top;


[nzbc,nxbc]=size(cm);



% calculate dt/dx/dens
b=dtx./den;
% damp is used in the absorbing boundary condition
vmin=min(min(vp(:),vs(:)));
damp=damp_circle(vmin,nzbc,nxbc,nbc,dx,isfs,pad_top);
temp=single(1-damp*dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Such scheme is abandoned.
% Adjust material parameter according tosave WD_test_iluu.mat vs_d vs vs1 cr_1 cr_0 residual residual_m dk_vs vp_d
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



%Back propogate modeling  by Zongcai
fd_order_num=fd_order;
input_vector = [nt,nzbc,nxbc,dtx,ng,nbc,gz(1),gx(1),gx(2)-gx(1),fd_order_num,in_wf,nz,nx,dt];
[cl_img,cm_img,illum_div] = eigen_e2drtm_single(input_vector,temp,ca,cl,cm,cm1,b,b1,seismo_w,wavefield_gradient.fux,wavefield_gradient.fuz,wavefield_gradient.bwx,wavefield_gradient.bwz);
end
