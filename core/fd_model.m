function [temp,ca,cl,cm,b]=fd_model(vp,vs,nbc,dx,dt,isfs,fd_order)
% Material padding and sponge shared by forward and adjoint propagation.
% Original elastic modeling: Bowen Guo (2014), Zongcai Feng (2015).
validateattributes(vs,{'numeric'},{'real','finite','positive','2d'});
validateattributes(vp,{'numeric'},{'real','finite','positive','size',size(vs)});
[ca,cm,cl]=calparam(vp,vs,ones(size(vs),'single'));
if isfs, pad_top=(fd_order-20)/2+1; else, pad_top=nbc; end
ca=single(pad(ca,nbc,isfs,pad_top));
cl=single(pad(cl,nbc,isfs,pad_top));
cm=single(pad(cm,nbc,isfs,pad_top));
b=repmat(single(dt/dx),size(cm));
damp=damp_circle(min([vp(:);vs(:)]),size(cm,1),size(cm,2),nbc,dx,isfs,pad_top);
temp=single(1-damp*dt);
if isfs
    ca(pad_top,:)=2*cm(pad_top,:);cl(pad_top,:)=0;
    b(pad_top,:)=2*b(pad_top,:);
end
end
