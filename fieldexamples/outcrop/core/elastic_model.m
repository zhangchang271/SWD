function [nz,nx,dx,dt,dtx,vp,vp_ss,vs,vs_ss,den] = elastic_model(vs,fr)

% if model_type==82
    %amp=0.5;
   [nz,nx]=size(vs);
   vp=vs*1.732;
%     load C_model.mat
%     vp=vp(8:32,11:180);
%     vs=vs(8:32,11:180);
    %nx=240*amp;
     dx=(min(vs(:))/fr/12);
 if dx<1
     dx=(min(vs(:))/fr/12);
 else
     dx=floor(min(vs(:))/fr/12);
 end
%      dx=0.25;
    %nz=52*amp;
  
    %dt=0.00025;
  
    x=(0:nx-1)*dx;z=(0:nz-1)*dx;
%  vp=ones(nz,nx)*1200;vs=ones(nz,nx)*700;
   den=ones(nz,nx);
   %vp_ss=vp; vs_ss=vs;
% nt=1000;
% generate source
  
%           for i=1:nz
%             for j=1:nx
%                if vs(i,j)==1000
%                   vs(i,j)=600;
%                end
%            end
%         end
%     [vp,~]=vel_smooth(vp,3,3,1);
%      [vs,~]=vel_smooth(vs,3,3,1);
  vp=floor(vp);
  vs=floor(vs);
  vmin=min(vs(:));
   vmax=max(vs(:));
%        for i=1:nx
%             vp(:,i)=linspace(vmax+700,vmax+900,nz);
%        end

% %     [vs_ss,~]=vel_smooth(vs,3,5,150);
    for i=1:nx
        vs_ss(:,i)=linspace(vmin,vmax,nz);
%         vp_ss(:,i)=linspace(500,3000,nz);
        %vs_ss(:,i)=linspace(vmin-200,vmin-200+10*nz,nz);
    end
   %  [vs_ss,~]=vel_smooth(vs,5,5,45);
%      vp=vs*2.4;
     %dt=dx/3000*0.5;
     dt=dx/max(vp(:))*0.5;
     dtx=dt/dx;
 %[vp,~]=vel_smooth(vp,3,3,5);
  vp_ss=vs_ss.*1.732;
%      vp_ss=vp;
%    for i=1:nx
%         vs_ss(:,i)=linspace(mean(vs(:)),mean(vs(:)),nz);
%     end

 
%     figure(1);
%     subplot(221);imagesc(x,z,vp);colorbar;  ylabel('Z (m)','fontsize',14);title('a) True P-wave velocity','fontsize',14);
%     subplot(223);imagesc(x,z,vp_ss);colorbar;xlabel('X (m)','fontsize',14); ylabel('Z (m)','fontsize',14);title('b) Start P-wave velocity','fontsize',14);
%     subplot(222);imagesc(x,z,vs);colorbar;title('c) S-wave Velocity','fontsize',14);
%     subplot(224);imagesc(x,z,vs_ss);colorbar; xlabel('X (m)','fontsize',14); title('d) Start S-wave velocity','fontsize',14);
end
% end



