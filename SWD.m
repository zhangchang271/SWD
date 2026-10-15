% Skeletonized Wave-Equation Dispersion Spectrum Inversion (SWD)
% Run this script for the physical-distance model example.
% Select WD or SWD with METHOD below. Parameters and workflow are here.
%
% Developers: Zhang Chang and Jing Li, Jilin University.
%
% License:
%
% GNU GENERAL PUBLIC LICENSE
% Version 3, 29 June 2007
%
% Copyright (C) 2024 Zhang Chang, Li jing and CSIM group
%
% This file is part of SWD.
%
% SWD is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%% Model and parameters
clc;
clear
close all;
method = 'WD';  % 'WD' or 'SWD'
root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir,'core'));
load(fullfile(root_dir,'model','model8_2.mat'));
vs_d = model8_2;
vsmin=min(vs_d(:));vsmax=max(vs_d(:));
vpmax=vsmax*1.732;
[nz,nx]=size(vs_d);
vs_d=single(vs_d);
vs = zeros(nz,nx); % Initial model
for i=1:nx
    vs(:,i)=linspace(vsmin,vsmax,nz);
end
vp=vs*1.732;  % vp is constant
fr=30;
dx=1;  % grid spacing in metres
dt=dx/vpmax*0.5;
dtx=dt/dx;
pickMethod=1;  %1==FDC 2==argmax
nt_wf=5;  % Storage wave field interval
nt=floor((nx*dx/vsmin/dt+1000)/10)*10-500;  % time step; Must be divisible by nt_wf
[s,~]=ricker(fr,dt,nt); s =single(s); % source wavelet
nbc=40;   % boundary layer
% define acquisition geometry
source_step=6;  % grid cells between shots
receiver_step=2;  % grid cells between receivers
sx=single(1:source_step:nx); sz=zeros(size(sx),'single')+1; [~,ns]=size(sx);
gx=single(1:receiver_step:nx); gz=zeros(size(gx),'single')+1; ng=numel(gx);
% Inversion settings
pur=-0.05;
offset=90;  % maximum source-receiver distance in metres
minimum_offset=6;  % stop when the offset window reaches this distance
offset_step=6;  % decrease the distance window by this amount
parameter_type=0;
fd_order=22; source_type='w';
isfs=1;
vmin=floor(vsmin/2);  % min phase-velocity of RT
vmax=floor(vsmax*1.2); % max phase-velocity of RT
np=vmax-vmin+1;
df=0.1;
fmin=10;  % min freq of the data
fmax=80;  % max freq of the data
ini=10;   % initial point to extract dispersion
minimum_traces=3;  % minimum traces required on one side
err=0.01;  % relative improvement threshold for reducing offset
SoftArgNorm = 1e+4;
iteration=190;
smoothZ = 1;smoothX=2;  % change if the gradient is incorrect

%% Observed gathers and inversion storage

parallel_init(ns);
vs_all=zeros(nz,nx,iteration); % Save the model of the iteration
dk_vs_all=zeros(nz,nx,iteration); % Save the gradient of the iteration
freq=(fmin:df:fmax);
npair=length(freq);

offsets = zeros(iteration,1);
residual = zeros(iteration,1);
space_M = repmat((1:np)',1,npair);

parfor is=1:ns
    [~,seismo_v_d(:,:,is)]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs_d,isfs,fd_order,source_type,parameter_type,nt_wf);
end
imagesc(seismo_v_d(:,:,1))

k =1;

%% Inversion
while (k<=iteration)
    tic;
    display(['Elastic_LSM, k=',num2str(k),' iteration=',num2str(iteration)]);

    ml = zeros(np,npair,ns);
    ml1 = zeros(np,npair,ns);

    cr_0 = 1.*ones(npair,ns);
    cr_0l = 1.*ones(npair,ns);

    parfor is=1:ns
        [ml(:,:,is),~]=RTrADx(seismo_v_d(:,:,is),df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
        [~,cr_0(:,is)] = LHDispPick(ml(:,:,is),vmin,cr_0(:,is),pickMethod,ini);
        [ml1(:,:,is),~]=RTlADx(seismo_v_d(:,:,is),df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
        [~,cr_0l(:,is)] = LHDispPick(ml1(:,:,is),vmin,cr_0l(:,is),pickMethod,ini);
    end

    cr_pre_r = 1.*ones(npair,ns);
    cr_pre_l = 1.*ones(npair,ns);
    g_cl=zeros(nz,nx);
    g_cm=zeros(nz,nx);
    g_illum=zeros(nz,nx);

    parfor is=1:ns
        [~,seismo_v,wavefield_gradient]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs,isfs,fd_order,source_type,parameter_type,nt_wf);
        [mlr,saveForBackwardr]=RTrADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
        [res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,vmin,cr_0(:,is),pickMethod,ini);

        [mll,saveForBackwardl]=RTlADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
        [res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,vmin,cr_0l(:,is),pickMethod,ini);
        switch method
            case 'WD'
                seismo_v_d1 = weight_dataAD(seismo_v,dt,dx,df,offset,sx(is),gx,minimum_traces,fmin,fmax, ...
                    cr_0(:,is),cr_0l(:,is),cr_pre_r(:,is),cr_pre_l(:,is));
            case 'SWD'
                grad_outputr = cr_pre_r(:,is)-cr_0(:,is);
                grad_outputl = cr_pre_l(:,is)-cr_0l(:,is);
                seismo_v_d1 = ADWDgrad_w(nt,ng,npair,minimum_traces,SoftArgNorm, ...
                    grad_outputr,grad_outputl,space_M,saveForBackwardr,saveForBackwardl);
            otherwise
                error('SWD:InvalidMethod','method must be ''WD'' or ''SWD''.');
        end
        [cl_img,cm_img,illum_div]=e2drtm_eigen(wavefield_gradient,single(seismo_v_d1),is,nbc,nt,dtx,dx,dt,gx,gz,s,vp,vs,isfs,fd_order,parameter_type,nt_wf);
        g_cl = g_cl+cl_img;g_cm = g_cm+cm_img;g_illum = g_illum+illum_div;
    end

    residual(k)=mean(mean(res_r));
    display(['residual = ',num2str( residual(k) ),' k=',num2str(k)]);
    res0=residual(k);
    g_cl=g_cl./g_illum;g_cm=g_cm./g_illum;
    dk_vs = -4*vs.*g_cl+2*vs.*g_cm;
    dk_vs=single(smooth2a(double(dk_vs),smoothZ,smoothX));  % Smooth the Vs gradient

    if k==1
        f1=0.5;
    end
    v_mean=(sum(vs(:).*vs(:)))^0.5;
    g_mean=(sum(dk_vs(:).*dk_vs(:)))^0.5;
    alpha=v_mean/g_mean*pur;
    display(['v_mean=',num2str(v_mean),' g_mean=',num2str(g_mean),' alpha=',num2str(alpha)]);
    vs1=vs+alpha*f1*dk_vs;
    vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax;

    parfor is=1:ns
        [~,seismo_v,~]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs1,isfs,fd_order,source_type,parameter_type,nt_wf);
        [mlr,~]=RTrADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
        [res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,vmin,cr_0(:,is),pickMethod,ini);
        [mll,~]=RTlADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
        [res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,vmin,cr_0l(:,is),pickMethod,ini);
    end
    res1=mean(mean(res_r));
    display(['f1= ',num2str(f1),' res1= ',num2str(res1)]);

    if res1>res0
        while res1>res0 && f1>0.0001
            f2=f1; res2=res1;
            f1=f1*0.5;
            vs1=vs+alpha*f1*dk_vs;
            vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax;

            parfor is=1:ns
                [~,seismo_v,~]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs1,isfs,fd_order,source_type,parameter_type,nt_wf);
                [mlr,~]=RTrADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
                [res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,vmin,cr_0(:,is),pickMethod,ini);
                [mll,~]=RTlADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
                [res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,vmin,cr_0l(:,is),pickMethod,ini);
            end
            res1=mean(mean(res_r));
            display(['f1= ',num2str(f1),' res1= ',num2str(res1)]);
        end
    else
        f2=f1*2;
        vs1=vs+alpha*f2*dk_vs;
        vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax;

        parfor is=1:ns
            [~,seismo_v,~]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs1,isfs,fd_order,source_type,parameter_type,nt_wf);
            [mlr,~]=RTrADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
            [res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,vmin,cr_0(:,is),pickMethod,ini);
            [mll,~]=RTlADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
            [res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,vmin,cr_0l(:,is),pickMethod,ini);
        end
        res2=mean(mean(res_r));
        display(['f2= ',num2str(f2),' res2= ',num2str(res2)]);
    end
    gama=(f1^2*(res0-res2)+f2^2*(res1-res0))/(2*res0*(f1-f2)+2*res1*f2-2*res2*f1);

    if isinf(gama)
        gama=0;
    end
    display(['gama= ',num2str(gama),' numerical step_length= ',num2str(gama*alpha)]);
    vs1=vs+alpha*gama*dk_vs;
    vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax;

    parfor is=1:ns
        [~,seismo_v,~]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx(is),sz(is),gx,gz,s,vp,vs1,isfs,fd_order,source_type,parameter_type,nt_wf);
        [mlr,~]=RTrADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
        [res_r(is),cr_pre_r(:,is)] = LHDispPick(mlr,vmin,cr_0(:,is),pickMethod,ini);
        [mll,~]=RTlADx(seismo_v,df,dt,dx,np,vmin,vmax,fmin,fmax,offset,gx,sx(is));
        [res_l(is),cr_pre_l(:,is)] = LHDispPick(mll,vmin,cr_0l(:,is),pickMethod,ini);
    end

    res3=mean(mean(res_r));
    display(['res3= ',num2str(res3)]);
    if (res3>res1 || res3>res2)
        if res1>res2
            res0=res2;
            gama=f2;
        else
            res0=res1;
            gama=f1;
        end
        vs1=vs+alpha*gama*dk_vs;
        vs1(vs1<vsmin)=vsmin;vs1(vs1>vsmax)=vsmax;

    else
        res0=res3;
    end
    vs=vs1;
    vs_all(:,:,k) = vs1;
    dk_vs_all(:,:,k) = dk_vs;
    offsets(k) = offset;

    if (k>1) && ((residual(k) - res0)/residual(k)) < err
        offset = max(minimum_offset,offset-offset_step);
        f1=0.5;
    end
    k = k +1;

    if offset<=minimum_offset
        break
    end
end

toc;

dfs = mean(mean((vs-vs_d).^2))
