function grid=fd_grid(terrain,gx,gz,model_size,nbc,dx,dt,isfs,fd_order,source_type,nt_wf)
%FD_GRID Prepare the fixed grid, receiver sampling, and propagation settings once.
% terrain: physical z at x=(0:nx-1)*dx, positive down. Zero is a flat surface.
% Receiver coordinates may be fractional and unordered. The map is fixed during inversion.
nz=model_size(1);nx=model_size(2);
if isempty(terrain), terrain=zeros(1,nx); end
validateattributes(terrain,{'numeric'},{'real','finite','vector','numel',nx});
validateattributes(gx,{'numeric'},{'real','finite','vector','nonempty','>=',1,'<=',nx});
terrain=double(terrain(:).')/dx;
r=(fd_order-20)/2; rows=nz+nbc+r+1; cols=nx+2*nbc; top=r+2;
% Align separated receivers with nearby nodes; close pairs still use interpolation.
p=unique(double(gx(:).')); node=round(p); keep=[true,diff(node)>0] & [diff(node)>0,true];
keep=keep & node>1 & node<nx;
anchor=[1,node(keep),nx]; position=[1,p(keep),nx];
% Closely spaced receivers are sampled directly instead of squeezing a cell.
while any(diff(position)./diff(anchor)<.5)
    close=diff(position)./diff(anchor)<.5;
    remove=[false,close]|[close,false];remove([1,end])=false;
    anchor(remove)=[];position(remove)=[];
end
xi=(1:cols)-nbc; eta=((1:rows)'-top);
X=interp1(anchor,position,min(max(xi,1),nx),'pchip');
X(xi<1)=xi(xi<1); X(xi>nx)=xi(xi>nx);
grid.mapped=any(terrain~=0)||any(abs(X-xi)>1e-12);
grid.terrain=terrain*dx; grid.x=X; grid.top=top; grid.nbc=nbc;
grid.nz=nz;grid.nx=nx;grid.dx=dx;grid.dt_scale=1;
grid.isfs=isfs;grid.fd_order=fd_order;grid.source_type=source_type;grid.nt_wf=nt_wf;
assert(~grid.mapped||isfs,'SWD:Surface','The terrain grid requires isfs=1.');
if grid.mapped
    assert(fd_order==22,'SWD:GridOrder','The terrain grid uses fd_order=22.');
    assert(nz>=4 && nbc>=3,'SWD:GridSize','The terrain grid needs nz>=4 and nbc>=3.');
    map=@(q) interp1(xi,X,q,'pchip','extrap');
    height=@(q) interp1(1:nx,terrain,min(max(map(q),1),nx),'pchip');
    H=nz-1;
    [mN,aN,cN,JN,slope]=metrics(xi,eta,map,height,H);
    [mS,aS,cS,JS]=metrics(xi+.5,eta+.5,map,height,H);
    [~,~,~,JU]=metrics(xi+.5,eta,map,height,H);
    [~,~,~,JW]=metrics(xi,eta+.5,map,height,H);
    assert(min([JN(:);JS(:);JU(:);JW(:)])>.2, ...
        'SWD:FoldedGrid','Terrain or receiver shifts make the grid too compressed.');
    invU=1./JU; invW=1./JW;
    JN(top,:)=JN(top,:)/2; invU(top,:)=2*invU(top,:);
    invU(1:top-1,:)=0;invW(1:top-1,:)=0;
    invU([1 end],:)=0;invU(:,[1 end])=0;
    invW([1 end],:)=0;invW(:,[1 end])=0;
    JN(1:top-1,:)=0;JS(1:top-1,:)=0;
    N=rows*cols; index=reshape(1:N,rows,cols); terms={};
    % Physical strain G on staggered nodes; its scatter transpose supplies the force.
    k=top+1:rows-1; i=2:cols-1; P=index(k,i);
    terms{end+1}=entries(P,P,aN(k,i));terms{end+1}=entries(P,P-rows,-aN(k,i));
    terms{end+1}=entries(P,P+1,-cN(k,i)/4);terms{end+1}=entries(P,P+1-rows,-cN(k,i)/4);
    terms{end+1}=entries(P,P-1,cN(k,i)/4);terms{end+1}=entries(P,P-1-rows,cN(k,i)/4);
    terms{end+1}=entries(P+N,P+N,mN(k,i));terms{end+1}=entries(P+N,P-1+N,-mN(k,i));
    k=top:rows-1;P=index(k,i);
    terms{end+1}=entries(P+2*N,P+1,mS(k,i));terms{end+1}=entries(P+2*N,P,-mS(k,i));
    terms{end+1}=entries(P+2*N,P+rows+N,aS(k,i));terms{end+1}=entries(P+2*N,P+N,-aS(k,i));
    % Central mixed derivative in the bulk, one-sided at the free surface.
    k=top+1:rows-1;P=index(k,i);
    terms{end+1}=entries(P+2*N,P+1+N,-cS(k,i)/4);terms{end+1}=entries(P+2*N,P+1+rows+N,-cS(k,i)/4);
    terms{end+1}=entries(P+2*N,P-1+N,cS(k,i)/4);terms{end+1}=entries(P+2*N,P-1+rows+N,cS(k,i)/4);
    P=index(top,i);a=cS(top,i)/2;
    terms{end+1}=entries(P+2*N,P+1+N,-a);terms{end+1}=entries(P+2*N,P+1+rows+N,-a);
    terms{end+1}=entries(P+2*N,P+N,a);terms{end+1}=entries(P+2*N,P+rows+N,a);
    % At the surface only tangential strain remains after eliminating zero traction.
    % The surface normal stiffness is ca-cl^2/ca, evaluated in the C++ kernel.
    a=aN(top,i)./(1+slope(top,i).^2); t=slope(top,i);
    terms{end+1}=entries(P,P,a);terms{end+1}=entries(P,P-rows,-a);
    terms{end+1}=entries(P,P+rows+N,.75*a.*t);terms{end+1}=entries(P,P-rows+N,-.75*a.*t);
    terms{end+1}=entries(P,P+1+rows+N,-.25*a.*t);terms{end+1}=entries(P,P+1-rows+N,.25*a.*t);
    entries_all=vertcat(terms{:});
    grid.G=sparse(entries_all(:,1),entries_all(:,2),entries_all(:,3),3*N,2*N);
    grid.JN=single(JN);grid.JS=single(JS);grid.invU=single(invU);grid.invW=single(invW);
    % A conservative stiffness bound gives the leapfrog time step without a user option.
    A=abs(grid.G); mass=sqrt([invU(:);invW(:)]);
    bound=max(mass.*(A'*([2*JN(:);2*JN(:);JS(:)].*(A*mass))));
    grid.dt_scale=min(1,3.6/sqrt(bound));
end
grid.dt=dt*grid.dt_scale;
[rx,rz]=fd_position(gx,gz,grid);
grid.receivers=[rz;rx];
end

function [m,a,c,J,s]=metrics(q,z,map,height,H)
h=height(q); hz=1-h/H;
xq=map(q+.5)-map(q-.5); hq=height(q+.5)-height(q-.5);
m=ones(size(z))*(1./hz); a=ones(size(z))*(1./xq);
c=(1-z/H)*(hq./(xq.*hz)); J=ones(size(z))*(xq.*hz);
s=ones(size(z))*(hq./xq);
end

function term=entries(row,col,value)
term=[row(:),col(:),value(:)];
end
