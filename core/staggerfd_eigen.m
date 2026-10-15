function [seismo_u,seismo_w,forward]=staggerfd_eigen(vp,vs,s,sx,sz,grid)
% Elastic forward propagation. The third output is passed directly to e2drtm_eigen.
% Positions retain grid units, with origin 1. Ordinary forward calls need two outputs.
[nz,nx]=size(vs);nt=numel(s);ng=size(grid.receivers,2);
assert(isequal([nz,nx],[grid.nz,grid.nx]),'SWD:ModelSize','Model and grid dimensions must agree.');
[sx,sz]=fd_position(sx,sz,grid);
[temp,ca,cl,cm,b]=fd_model(vp,vs,grid.nbc,grid.dx,grid.dt,grid.isfs,grid.fd_order);
types={'p','s','s1','z','w'};type=find(strcmp(grid.source_type,types));
assert(~isempty(type),'SWD:Source','Unknown source_type.');
assert(~grid.mapped||type==5,'SWD:Source','The terrain grid currently uses the original w source.');
dtx=grid.dt/grid.dx;
input=[nt,size(cm,1),size(cm,2),dtx,ng,sz,sx,0,0,0,type,grid.fd_order,grid.nt_wf,nz,nx,3];
if nargout>2
    path=[tempname,'.swd'];
    try
        [seismo_u,seismo_w,forward]=eigen_staggerfd_single(input,temp,ca,cl,cm,cm,b,b,s,grid.receivers,grid,path);
    catch exception
        if isfile(path), delete(path); end
        rethrow(exception);
    end
    % Retain exactly the model and settings that generated this shot.
    forward.input=[nt,size(cm,1),size(cm,2),dtx,ng,grid.nbc,0,0,0,grid.fd_order,grid.nt_wf,nz,nx,grid.dt];
    forward.temp=temp;forward.ca=ca;forward.cl=cl;forward.cm=cm;forward.b=b;
    forward.s=s;forward.grid=grid;
else
    [seismo_u,seismo_w]=eigen_staggerfd_single(input,temp,ca,cl,cm,cm,b,b,s,grid.receivers,grid,'');
end
end
