function [cl_img,cm_img,illum_div]=e2drtm_eigen(wavefield,seismo_w,is,nbc,nt,dtx,dx,dt,gx,gz,s,vp,vs,isfs,fd_order,parameter_type,in_wf,grid)
% Reconstruct the source wavefield while propagating the receiver adjoint.
[nz,nx]=size(vs);
if nargin<18, grid=fd_grid(zeros(1,nx),gx,nz,nx,nbc,dx,fd_order); end
[rx,rz]=fd_position(gx,gz,grid);
[temp,ca,cl,cm,b]=fd_model(vp,vs,nbc,dx,dt,isfs,fd_order,parameter_type);
input=[nt,size(cm,1),size(cm,2),dtx,numel(gx),nbc,0,0,0,fd_order,in_wf,nz,nx,dt];
cleanup=onCleanup(@()delete(wavefield.boundary_file));
[cl_img,cm_img,illum_div]=eigen_e2drtm_single(input,temp,ca,cl,cm,cm,b,b,seismo_w,wavefield,[rz;rx],s,grid);
end
