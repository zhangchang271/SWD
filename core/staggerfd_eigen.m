function [seismo_u,seismo_w,wavefield]=staggerfd_eigen(is,nbc,nt,dtx,dx,dt,sx,sz,gx,gz,s,vp,vs,isfs,fd_order,source_type,parameter_type,in_wf,grid)
% Elastic forward propagation. Positions retain grid units, with origin 1.
% The third output contains a terminal state and boundary history, not volume snapshots.
[nz,nx]=size(vs);
if nargin<19, grid=fd_grid(zeros(1,nx),gx,nz,nx,nbc,dx,fd_order); end
assert(~grid.mapped||isfs,'SWD:Surface','The terrain grid requires isfs=1.');
[rx,rz]=fd_position(gx,gz,grid);[sx,sz]=fd_position(sx,sz,grid);
[temp,ca,cl,cm,b]=fd_model(vp,vs,nbc,dx,dt,isfs,fd_order,parameter_type);
types={'p','s','s1','z','w'};type=find(strcmp(source_type,types));
assert(~isempty(type),'SWD:Source','Unknown source_type.');
assert(~grid.mapped||type==5,'SWD:Source','The terrain grid currently uses the original w source.');
input=[nt,size(cm,1),size(cm,2),dtx,numel(gx),sz,sx,0,0,0,type,fd_order,in_wf,nz,nx,3];
if nargout>2
    path=[tempname,'.swd'];
    try
        [seismo_u,seismo_w,wavefield]=eigen_staggerfd_single(input,temp,ca,cl,cm,cm,b,b,s,[rz;rx],grid,path);
    catch exception
        if isfile(path), delete(path); end
        rethrow(exception);
    end
else
    [seismo_u,seismo_w]=eigen_staggerfd_single(input,temp,ca,cl,cm,cm,b,b,s,[rz;rx],grid,'');
end
end
