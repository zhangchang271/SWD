function [ca,cm,cl]=calparam(vp,vs,den)

ca=(vp.^2).*den;
cm=(vs.^2).*den;
cl=ca-2*cm;

end