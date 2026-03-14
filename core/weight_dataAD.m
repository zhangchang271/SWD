function [seismo_v_d1]=weight_dataAD(seismo_v,dt,dx,df,offset,sx,gx,w,fmin,fmax,cr_0,cr_0l,cr_pre_r,cr_pre_l)



[nt,ng]=size(seismo_v);
seismo_v_d1=zeros(nt,ng); % Define the backprogate data

fre=2*pi*linspace(fmin,fmax,(fmax-fmin)/df+1); 
freq=(fmin:df:fmax);   
ind=fmin/df+1;
npair=length(freq);
win=npair+2*ind;
 %% Calcuate the weight source
 ntpad=floor(1/dt/df);
 %spad=zeros(ntpad);  
 % Caluate the weight residual data

xp = (gx-sx)*dx;
rightdata = xp>0;
xp = xp(rightdata);
rangeoffset = xp<=offset;
x = xp(rangeoffset);
y = x(:)';
uxtposr = rightdata&(abs((sx-gx)*dx)<=offset);
offsetr = sum(uxtposr);
if offsetr>=w
    deta_cc2=zeros(win,1);
    deta_cc2(ind:ind+npair-1) = fre'.*(cr_0-cr_pre_r)./(cr_0.*cr_pre_r);
    deta_cc2(isnan(deta_cc2))=0;
    tempdata=zeros(ntpad,offsetr);

    tempdata(1:nt,:)=seismo_v(:,uxtposr);    

    deta_seismo_v1(:,1:offsetr)=fft(tempdata); %TRansform shot gather to FK domain
    y1(:)=deta_cc2(:); % Define the delta k
    
    for ii=1:offsetr

         deta_seismo_v4(:,ii)=(-1i*y1(:)*y(ii))/pi/2;   %% predicted data

    end

    for omega=1:win
        deta_seismo_v1(omega+1,:)=deta_seismo_v1(omega+1,:).*deta_seismo_v4(omega,:);
        deta_seismo_v1(ntpad+1-omega,:)=conj(deta_seismo_v1(omega+1,:));
    end
    deta_seismo_vv=(ifft(deta_seismo_v1(:,:)));
    seismo_v_d1(:,uxtposr)=deta_seismo_vv(1:nt,1:offsetr);
 
end



%% Left side part
xp = (sx-gx)*dx;
leftdata = xp>0;
xp = xp(leftdata);
rangeoffset = xp<=offset;
x1 = xp(rangeoffset);
y = (x1)';
uxtposl = leftdata&(abs((sx-gx)*dx)<=offset);
offsetl = sum(uxtposl);

if offsetl>=w
    deta_cc22 = zeros(win,1);
    deta_cc22(ind:ind+npair-1) = fre'.*(cr_0l-cr_pre_l)./(cr_0l.*cr_pre_l);
    deta_cc22(isnan(deta_cc22))=0;

    tempdata=zeros(ntpad,offsetl);

    tempdata(1:nt,:)=seismo_v(:,uxtposl);  

    deta_seismo_v11(:,1:offsetl)=fft(tempdata); %TRansform shot gather to FK domain
    y1(:)=deta_cc22(:); % Define the delta k
    
    for ii=1:offsetl
        deta_seismo_v44(:,ii)=(-1i*y1(:)*y(ii))/pi/2;   %% predicted data

    end

    for omega=1:win
        deta_seismo_v11(omega+1,:)=deta_seismo_v11(omega+1,:).*deta_seismo_v44(omega,:);
        deta_seismo_v11(ntpad+1-omega,:)=conj(deta_seismo_v11(omega+1,:));
    end
    
   deta_seismo_vv=(ifft(deta_seismo_v11(:,:)));
   seismo_v_d1(:,uxtposl)=deta_seismo_vv(1:nt,1:offsetl);
end
%clear deta_seismo_vv deta_seismo_v1 deta_seismo_v11 deta_seismo_v4 deta_seismo_v44 deta_cc deta_cc1 y y1 tempdata
end