function test_fd
% Numerical checks, separate from the single SWD.m application entry point.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'core'));
build_mex;
f=load(fullfile(root,'tests','fixtures','flat_reference.mat'));
[nz,nx]=size(f.vs); nbc=f.nbc; dx=f.dx; nt=f.nt;
vp=f.vp; vs=f.vs; dt=f.dt; s=f.s;
[u,w,tape]=staggerfd_eigen(1,nbc,nt,dt/dx,dx,dt,f.sx,f.sz,f.gx,f.gz,s,vp,vs,1,22,'w',0,5);
path=tape.boundary_file; assert(isfile(path));
[gl,gm,il]=e2drtm_eigen(tape,f.res,1,nbc,nt,dt/dx,dx,dt,f.gx,f.gz,s,vp,vs,1,22,0,5);
assert(~isfile(path),'Boundary file was not removed.');
assert(relative(u,f.u)<1e-6 && relative(w,f.w)<1e-6);
assert(relative(gl,f.cl)<2e-4 && relative(gm,f.cm)<2e-4 && relative(il,f.illum)<1e-5);
fprintf('Original flat-grid data and gradients: passed.\n');

% Unordered, fractional, repeated receivers, including receivers below the surface.
terrain=1.2*sin((0:nx-1)*pi/40);
gx=single([70.2 4.1 12.3 28.2 35.7 44.4 62.9]);
gz=1+interp1(1:nx,terrain,double(gx),'pchip');
gz(3)=gz(3)+1.3;
sx=single(20.2); sz=1+interp1(1:nx,terrain,double(sx),'pchip');
grid=fd_grid(terrain,gx,nz,nx,nbc,dx,22);
assert(grid.mapped && grid.dt_scale>0 && grid.dt_scale<=1);
N=size(grid.G,2)/2;
assert(norm(grid.G*[ones(N,1);zeros(N,1)])<1e-10);
assert(norm(grid.G*[zeros(N,1);ones(N,1)])<1e-10);
dt=double(grid.dt_scale)*.5/double(max(vp(:))); s=single(ricker(30,dt,nt));
[~,w,tape]=staggerfd_eigen(1,nbc,nt,dt/dx,dx,dt,sx,sz,gx,gz,s,vp,vs,1,22,'w',0,1,grid);
rng(11); residual=single(randn(size(w)));
[gl,gm,~]=e2drtm_eigen(tape,residual,1,nbc,nt,dt/dx,dx,dt,gx,gz,s,vp,vs,1,22,0,1,grid);
permutation=[7 3 1 6 2 5 4 3];
other=fd_grid(terrain,gx(permutation),nz,nx,nbc,dx,22);
[~,wr]=staggerfd_eigen(1,nbc,nt,dt/dx,dx,dt,sx,sz,gx(permutation),gz(permutation),s,vp,vs,1,22,'w',0,1,other);
assert(relative(wr,w(:,permutation))<1e-6);
[~,one]=staggerfd_eigen(1,nbc,nt,dt/dx,dx,dt,sx,sz,gx(3),gz(3),s,vp,vs,1,22,'w',0,1,grid);
assert(relative(one,w(:,3))<1e-6);

% Model perturbation is internal: keep the padded material and sponge fixed.
dv=single(smooth2a(randn(size(vs)),2,2));dv=dv/max(abs(dv(:)));
dv(1:3,:)=0;dv(end-2:end,:)=0;dv(:,1:3)=0;dv(:,end-2:end)=0;
analytic=sum(double((-4*vs.*gl+2*vs.*gm).*dv),'all');
[~,wp]=staggerfd_eigen(1,nbc,nt,dt/dx,dx,dt,sx,sz,gx,gz,s,vp,vs+dv,1,22,'w',0,1,grid);
[~,wm]=staggerfd_eigen(1,nbc,nt,dt/dx,dx,dt,sx,sz,gx,gz,s,vp,vs-dv,1,22,'w',0,1,grid);
numeric=dt*sum(double(wp(:)-wm(:)).*double(residual(:)))/2;
error=abs(numeric-analytic)/max(abs(analytic),1e-30);
assert(error<2e-3,'Mapped-grid model derivative failed.');
fprintf('Flexible positions and mapped-grid derivative: passed (%.3g).\n',error);

% The test kernel stores a full trajectory solely to check reconstruction.
folder=tempname; mkdir(folder); addpath(folder);
cleanup=onCleanup(@()remove_test_kernel(folder));
compiler={};
if isunix && ~ismac && isfile('/usr/bin/g++-9'), compiler={'CXX=/usr/bin/g++-9'}; end
include='/usr/include/eigen3';
if ~isfile(fullfile(include,'Eigen','Core')), include=fullfile(root,'core','eigen3'); end
mex('-R2018a',['-I',include],'-DNDEBUG',compiler{:}, ...
    fullfile(root,'tests','fd_reconstruction_check.cpp'),'-outdir',folder);
types={'p','s','s1','z','w'};
max_error=0;
for mapped=[false true]
    for order=[22 24 26 28]
        if mapped && order~=22, continue; end
        if mapped
            geom=grid; x=gx; z=gz; px=sx; pz=sz; step=dt;
        else
            x=f.gx;z=f.gz;px=single(1);pz=single(1);step=f.dt;
            geom=fd_grid(zeros(1,nx),x,nz,nx,nbc,dx,order);
        end
        for type=1:5
            if mapped && type~=5, continue; end
            [temp,ca,cl,cm,b]=fd_model(vp,vs,nbc,dx,step,1,order,0);
            [rx,rz]=fd_position(x,z,geom); [qx,qz]=fd_position(px,pz,geom);
            input=[nt,size(cm,1),size(cm,2),step/dx,numel(x),qz,qx,step,0,0,type,order,5,nz,nx,3];
            source=single(ricker(30,step,nt)); r=single(randn(nt,numel(x)));
            [errors,cl0,cm0,il0]=fd_reconstruction_check(input,temp,ca,cl,cm,cm,b,b,source,[rz;rx],geom,r);
            [~,~,tape]=staggerfd_eigen(1,nbc,nt,step/dx,dx,step,px,pz,x,z,source,vp,vs,1,order,types{type},0,5,geom);
            [cl1,cm1,il1]=e2drtm_eigen(tape,r,1,nbc,nt,step/dx,dx,step,x,z,source,vp,vs,1,order,0,5,geom);
            assert(max(errors)<1e-4 && relative(cl1,cl0)<1e-3 && relative(cm1,cm0)<1e-3 && relative(il1,il0)<1e-4);
            max_error=max(max_error,max(errors));
        end
    end
end
fprintf('Every reconstructed time step and full-storage gradients: passed (max %.3g).\n',max_error);
rayleigh_check;
fprintf('All finite-difference checks passed.\n');
end

function error=relative(a,b)
error=norm(double(a(:))-double(b(:)))/max(norm(double(b(:))),1e-30);
end

function remove_test_kernel(folder)
clear fd_reconstruction_check;
rmpath(folder);rmdir(folder,'s');
end

function rayleigh_check
% A rotated homogeneous half-space has the same Rayleigh speed along its surface.
nx=241;nz=81;nbc=24;vs=single(ones(nz,nx)*300);vp=single(vs*sqrt(3));
gx=single(61:2:201);sx=single(21);nt=2400;
ratio=fzero(@(q)q^6-8*q^4+(24-16/3)*q^2-16*(1-1/3),[.8 .99]);
for slope=[0 .3]
    terrain=1+slope*((0:nx-1)-(nx-1)/2);
    grid=fd_grid(terrain,gx,nz,nx,nbc,1,22);
    dt=double(.5/max(vp(:)))*grid.dt_scale;s=single(ricker(15,dt,nt));
    gz=1+interp1(1:nx,terrain,double(gx),'pchip');sz=1+interp1(1:nx,terrain,double(sx),'pchip');
    [~,w]=staggerfd_eigen(1,nbc,nt,dt,1,dt,sx,sz,gx,gz,s,vp,vs,1,22,'w',0,5,grid);
    expected=300*ratio/sqrt(1+slope^2);
    distance=double(gx-sx);velocity=(210:340)';t=(0:nt-1)'*dt;
    for j=1:numel(distance)
        center=distance(j)/expected+.075;
        w(:,j)=w(:,j).*single(exp(-((t-center)/.12).^8));
    end
    frequencies=12:20;energy=zeros(numel(velocity),numel(frequencies));
    for j=1:numel(frequencies)
        f=frequencies(j);spectrum=exp(-2i*pi*f*t')*double(w);
        energy(:,j)=abs(exp(2i*pi*f*(1./velocity)*distance)*spectrum.');
    end
    [~,pick]=max(energy); measured=median(velocity(pick));
    assert(abs(measured-expected)/expected<.02,'Rayleigh-wave speed check failed.');
    fprintf('Surface slope %.1f: Rayleigh speed %.1f, analytical %.1f m/s.\n',slope,measured,expected);
end
end
