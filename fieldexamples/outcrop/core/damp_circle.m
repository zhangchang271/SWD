% To form a circle for abc damping
% @version 1 2014/09/28
% @author Bowen Guo

function [damp]=damp_circle(vmin,nzbc,nxbc,nbc,dx,isfs,pad_top)


nz=nzbc-nbc-pad_top;
nx=nxbc-2*nbc;

a=(nbc-1)*dx;
kappa = 3.0 * vmin * log(10000000.0) / (2.0 * a);

% setup 1D BC damping array
damp1d=kappa*((0:(nbc-1))*dx/a).^2;

damp=zeros(nzbc,nxbc);


% Take care of the 4 boundaries
   % Left and right
for iz=1+pad_top:nz+pad_top
    damp(iz,1:nbc)=damp1d(nbc:-1:1);
    damp(iz,nx+nbc+1:nx+2*nbc)=damp1d(:);
end

for ix=1+nbc:nx+nbc
    if (isfs)
       damp(1:pad_top,ix)=0.0;
    else
       damp(1:pad_top,ix)=damp1d(nbc:-1:1);
    end
    damp(nzbc-nbc+1:nzbc,ix)=damp1d(:);
end

% Take care of the 4 corners
    % Upper left
if (~isfs)
   for iz=1:pad_top
      for ix=1:nbc
         dist=sqrt((ix-nbc-1)^2+(iz-pad_top-1)^2);
         damp(iz,ix)=kappa*(dist/nbc).^2;
      end
   end
     % Upper right
   for iz=1:pad_top
      for ix=nx+nbc+1:nxbc
         dist=sqrt((ix-nx-nbc)^2+(iz-pad_top-1)^2);
         damp(iz,ix)=kappa*(dist/nbc).^2;
       end
   end
end
   
   % Lower left
for iz=nz+pad_top+1:nzbc
    for ix=1:nbc
       dist=sqrt((ix-nbc-1)^2+(iz-nz-pad_top)^2);
       damp(iz,ix)=kappa*(dist/nbc).^2;
    end
end
   % Lower right
for iz=nz+pad_top+1:nzbc
    for ix=nx+nbc+1:nxbc
       dist=sqrt((ix-nbc-nx)^2+(iz-nz-pad_top)^2);
       damp(iz,ix)=kappa*(dist/nbc).^2;
    end
end



end


  
    
    
    
  