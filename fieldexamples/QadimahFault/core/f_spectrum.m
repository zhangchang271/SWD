% To display the frequency spectrum of a time series vector
% @version 1 2014-10-02
% @author Bowen Guo

function f_spectrum(v,dt)

[nt,ng]=size(v);

f=fft(v);

ff=1/dt/nt*(0:floor(nt/10)-1);


% if (ng>1)
%    figure;plot(ff,sum(abs(f(1:floor(nt/2),:)),2));xlabel('frequency','fontsize',14);   ylabel('Amplitude','fontsize',14);
%    figure;imagesc((1:ng),ff,(abs(hilbert(f(1:floor(nt/2),:)))));
%    xlabel('trace No.','fontsize',14);
%    ylabel('frequency','fontsize',14);
% else
   figure;plot(ff,abs(f(1:floor(nt/10))));
   xlabel('frequency','fontsize',14);
   ylabel('Amplitude','fontsize',14);
% end


end



