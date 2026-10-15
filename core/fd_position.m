function [x,z]=fd_position(x,z,grid)
% Physical grid coordinates (origin 1) to the fixed computational grid.
validateattributes(x,{'numeric'},{'real','finite','vector','>=',1,'<=',grid.nx});
validateattributes(z,{'numeric'},{'real','finite','vector'});
x=double(x(:).');z=double(z(:).');
if isscalar(z), z=repmat(z,size(x)); end
assert(numel(x)==numel(z),'SWD:PositionSize','Each position needs an x and a z.');
if grid.mapped
    h=interp1(1:grid.nx,double(grid.terrain)/grid.dx,x,'pchip');
    eta=(z-1-h)./(1-h/(grid.nz-1));
    tol=8*double(eps('single'))*max(1,max(abs([z,h])));
    assert(all(eta>=-tol & eta<=grid.nz-1+tol), ...
        'SWD:PositionDepth','Source or receiver lies outside the terrain-following model.');
    z=min(max(eta,0),grid.nz-1)+1;
    lo=ones(size(x));hi=lo*grid.nx;
    for j=1:40
        mid=(lo+hi)/2;
        value=interp1((1:numel(grid.x))-grid.nbc,grid.x,mid,'pchip');
        left=value<x;lo(left)=mid(left);hi(~left)=mid(~left);
    end
    x=(lo+hi)/2;
else
    assert(all(z>=1 & z<=grid.nz),'SWD:PositionDepth','Position is outside the model.');
end
x=x+grid.nbc;z=z+grid.top-1;
end
