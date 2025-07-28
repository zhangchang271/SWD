function [w,tw] = ricker(f,dt,nt)
%RICKER: Ricker wavelet of central frequency f.
%
%  [w,tw] = ricker(f,dt);
%
%  IN   f : central freq. in Hz (f <<1/(2dt) )
%       dt: sampling interval in sec
%       nt: number of sample
%
%  OUT  w:  the Ricker wavelet
%       tw: axis
%
%  Example
%
%    [w,tw] = ricker(10,0.004,101);
%    plot(tw,w);

nw=2.2/f/dt;
nw=2*floor(nw/2)+1;
nc=floor(nw/2);
%w0 = zeros(nw,1);

k=(1:1:nw)';

alpha = (nc-k+1).*f*dt*pi;
beta=alpha.^2;
w0 = (1.-beta.*2).*exp(-beta);

if nargin==3
    if nt<length(w0)
        error('nt is smaller than condition!');
    else
        w=zeros(nt,1);
        w(1:length(w0))=w0;
    end
else
    w=w0;
end

nw=length(w);

if nargout>1;
    tw=(0:nw-1)*dt;
    %     tw = -(nc+1-(1:1:nw))*dt;
end
end