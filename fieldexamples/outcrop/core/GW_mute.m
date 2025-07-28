function [seismo_w_mute,zoffset_mask]=GW_mute(seismo_w,T1,gx,sx,T0,vg,dt,ztr)
zoffset_mask=mute_near_offset(seismo_w,gx-sx,ztr); 
tau_o=T0+ceil(abs(gx-sx)/vg/dt); 
[nt,ng] = size(seismo_w);
tapering = ones(size(seismo_w));  
for k = 1:ng
  %  b2 = tau_o(k) + 0.1*T1; 
      b2 = tau_o(k);
    sig = 0.05* T1;   
  for j = 1:nt
      if j> b2      
          tapering(j,k) = exp(-(j-b2)^2/(2*sig^2)) ;
      end
  end
end    
%tapering=tapering.*zoffset_mask; 
seismo_w_mute = seismo_w.*(1-tapering);    
end
%% 
function [output]=mute_near_offset(input,h,lamba)
output=ones(size(input)); 
dh=h(2)-h(1);
idh=abs(h/dh);
[~,ix2]=find(idh==lamba);
if ~isempty(ix2)
    for ix=1:length(h)
        if (idh(ix)<=lamba)
           output(:,ix)=0; 
        end 
    end 
end 
end 

