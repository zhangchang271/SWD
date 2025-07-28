function [dfk] = fk_filter(d, dt, dx)

%
%     Usage : [dfk, fktype, vrange] = fk_filter(d, dt, dx) 
%   
%    Inputs : 
%         d : The 2-D input data matrix (GPR section)
%        dt : Sampling rate along columns (time-dimension)
%        dx : Sampling rate along rows (trace spacing in the space
%             dimension) 
%
%    Output : 
%       dfk : The F-K filtered data matrix (GPR section)
%    fktype : Flag indicating the type of filtering performed. 
%    vrange : Two-column matrix containing the vertices of the polygonal
%             pass / stop zones when fktype = 1, 2, or, the velocity pass /
%             stop ranges when fktype = 3, 4. (fan filtering). Returned as
%             an empty variable when fktype = 5, 6 (up-dipping /
%             down-dipping event separation). 
%
%  Requires : setpolygon.m, isinpoly.m
%

% Get environmental variables
% global ENVAR
% vrange = [];
%%%%%    Inquire filter type
%fktype = menu('SELECT FILTER TYPE','Polygonal-Zone Pass',...
%     'Polygonal-Zone Stop','Velocity range pass', 'Velocity range stop', ...
%     'Up-dip', 'Down-dip', 'Cancel');
% if fktype == 7,
%     dfk = [];
%     return
% end
% Determine whether the pass / stop zone coordinates are predefined (import
% from file), or will be designed on-screen
% if fktype == 1 || fktype == 2,
%     inmode = menu('PLEASE DECIDE INPUT METHOD',...
%         'Import Zone Coordinates from File','Design Filter on Screen',...
%         'Cancel');
%     if inmode == 3,
%         dfk = [];
%         return
%     end
% else
%     inmode = 0;
% end
% Pad with lots of zeros to twice the size of the input data array
[ns, ntr] = size(d);
D = [d; zeros(ns, ntr)];
D = [D zeros(2*ns, ntr)];
D = fft2(D);
D = fftshift(D);
%-----------------------------------------------------------------------
%%%%%   Useful parameters
[nf,nk] = size(D);
nfhalf  = nf/2;
fn      = 1./(2*dt);                 % Nyquist frequency
df      = 2*fn/(nf - 1);             % frequency interval
ff      = -fn:df:fn;                 % frequency spectrum
nkhalf  = nk/2;
if dx < 0,
    dx = -dx;
end
kn      = 1./(2*dx);                 % Nyquist wavenumber
dk      = 2*kn/(nk - 1);             % wavenumber interval
kk      = -kn:dk:kn;                 % wavenumber spectrum           

% switch fktype
% case {1, 2, 3, 4}                    % The filter is zone-pass or zone-stop
%     %%%%%%% Use adaptive figure sizing and posistioning  %%%%%%%%%%%%%%%%%%
%     scrsz  = get(0,'screensize');
%     ffkpos = round([440*scrsz(3)/1680 200*scrsz(4)/1050 ...
%         800*scrsz(3)/1680 700*scrsz(4)/1050]);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ffk = figure('Numbertitle','off','tag','fkfiltfig',...
%         'Name','F-K spectrum','menubar','none', ...
%     'position', ffkpos);
% % Display every second element of the F-K spectrum for faster rendering
%     pcolor(kk(1:2:nk),ff(nfhalf+1:2:nf),log10(abs(D(nfhalf+1:2:nf,1:2:nk)))); 
%     shading('flat');
%     set(gca,'ydir','reverse');
%     if ~isempty(ENVAR) && ~isempty(ENVAR.colormap),
%         colormap(ENVAR.colormap);
%     end
% %    imagecolors;
%     ylabel('Frequency ( GHz )');
%     xlabel('Wavenumber (m^-^1)')
%     set(gca,'nextplot','add')

%%% Filter is Zone pass / stop and zone data imported from disk file
%     if inmode == 1 && (fktype == 1 || fktype ==2),
%         [kfname, kfpath]= uigetfile('*.dat;*.txt',...
%             'Please give file with K, F coordinates of Pass/Stop Zone', ...
%             ENVAR.currentworkdir);
%         if kfname == 0,                           % Canceled
%             delete(ffk);
%             dfk = [];
%             return;  
%         end; 
%         poly = importdata([kfpath kfname]);   % get k,f vertices of polygon
%         k = poly(:,1);
%         f = poly(:,2);
%         k       = [ k ; k(1)];
%         f       = [ f ; f(1)];                % close polygon
%         vrange  = poly;
%         clear poly
%         patch3 = patch(k,f,'k');
%         set(patch3,'facealpha',0.1)
% %
%         reply = questdlg('Is this OK ? ');
%         if strcmp(reply,'Cancel') || strcmp(reply,'No'),
%             delete(ffk); 
%             dfk = [];
%             vrange = [];
%             return
%         end
% 
% %%% Filter is Zone pass/stop and zone will be designed on screen
%     elseif inmode == 2 && (fktype == 1 || fktype ==2),
%         % Tell user how to enter the pass/stop zones using mouse
%         txt = cell(1);
%         txt(1) = cellstr('To set vertices of polygonal pass / stop zone,');
%         txt(2) = cellstr('click LEFT mouse button. To correct mistakes ');
%         txt(3) = cellstr('click the MIDDLE button. When finished, click ');
%         txt(4) = cellstr('the RIGHT button.');
%         msg = msgbox(txt);
%         pos=get(msg,'position');
%         set(msg,'pos',[20 20 pos(3) pos(4)]);
%         pause(0.5)
%         reply = 'No';
%         while strcmp(reply,'Yes')~=1,
%             figure(ffk)
%         % k, f are the wavenumber and frequency coordinates of the polygon
%         % vertices respectively
%             [k, f] = setpolygon(' m^(-1)', ' GHz');
%             if isempty(k) || isempty(f),
%                 delete(ffk);
%                 if exist('msg') && ishandle(msg),
%                     delete(msg)
%                 end
%                 dfk = [];
%                 return
%             end
%             k       = [ k ; k(1)];
%             f       = [ f ; f(1)];           % close polygon
%             vrange  = [k(:) f(:)];
%             patch3 = patch(k,f,'k');
%             set(patch3,'facealpha',0.1)
%             %%%%%   Make sure that pass-zone has been set to satisfaction
%             reply = questdlg('Is this OK ? ');
%             if strcmp(reply,'Cancel')==1,
%                 delete(ffk); 
%                 if exist('msg') && ishandle(msg),
%                     delete(msg)
%                 end
%                 dfk = [];
%                 return
%             end
%             if strcmp(reply,'No')==1,
%                 delete(patch3);
%                 vrange = [];
%                 pause(0.5)
%             end
%         end
%     end
            
%%% Filter is velocity-range pas/stop - use dialog boxes only
%     if fktype == 3 || fktype == 4,
%         reply = 'No';
%         while strcmp(reply,'Yes')~=1,
%             vrange = zeros(4,1);
%             clb = cell(4,1);          
%             clb(1) = cellstr('Give UPPER negative velocity in m/ns');
%             clb(2) = cellstr('Give LOWER negative velocity in m/ns');
%             clb(3) = cellstr('Give LOWER positive velocity in m/ns');
%             clb(4) = cellstr('Give UPPER positive velocity in m/ns');
%             for i=1:4
%                 cldef(i) = cellstr(num2str(vrange(i)));
%             end
%             answer = inputdlg(clb,'Please give velocities',1 );  
%             if isempty(answer),
%                 delete(ffk); 
%                 dfk = [];  fktype=[];  vrange=[];
%                 return
%             end
%             lb = char(answer);                                                  
%             for i=1:4
%                 comma = findstr(lb(i,:),',');
%                 if ~isempty(comma), 
%                     for j=1:length(comma); 
%                         lb(i,comma) = '.'; 
%                     end;
%                     clear comma
%                 end
%                 dummy = str2num(lb(i,:));
%                 vrange(i) = dummy(1); 
%             end
% %%%%%   Make absolutely sure no bullshit has been entered
%             v1 = min(abs(vrange(1)),abs(vrange(2)));
%             v2 = max(abs(vrange(1)),abs(vrange(2)));
%             vrange(1) = -1.*v1;
%             vrange(2) = -1.*v2;
%             v1 = min(abs(vrange(3)),abs(vrange(4)));
%             v2 = max(abs(vrange(3)),abs(vrange(4)));
%             vrange(3) = v2;
%             vrange(4) = v1;
% %%%%%   Define lines of constant velocity and associated polygon
%             k = [0  fn/vrange(1)  fn/vrange(2)  0 ...
%                 fn/vrange(3)   fn/vrange(4)   0]';
%             f = -1*[0  kn*vrange(1)  kn*vrange(2)  0 ...
%                 -kn*vrange(3)  -kn*vrange(4)   0]';
%             patch3 = patch(k,f,'k'); 
%             set(patch3,'facealpha',0.1)
%                     
% %%%%%   Make sure that pass-zone has been set to satisfaction
%             reply = questdlg('Is this OK ? ');
%             if strcmp(reply,'Cancel')==1,
%                 delete(ffk); 
%                 dfk = [];    
%                 return
%             end
%             if strcmp(reply,'No')==1,
%                 delete(patch3); 
%                 vrange = [];
%                 pause(0.5)
%             end
%         end                           % while reply   
%     end                            % fktype       
%%% Clear first message to make room for the second!
%     if exist('msg') && ishandle(msg),
%         delete(msg)
%     end
%%% Pass/stop zone defined - proceed 
    dfk     = D; 
%%% Set up wavenumber and frequency grids for "inpolygon" and "isinpoly"
    [waven, freq] = meshgrid(kk, ff);
%%% Tell user to wait - lots of work to do...    
%     msg = msgbox('Working ... Please wait!','FK_FILTER: INFO','help');
% %%% The MATLAB function "inpolygon" may not exist for MATLAB releases
% %%% earlier that 13! Check and if not, switch to Kirill Pankratov's
% %%% "isinpoly".  
%     find_inpolygon = which('inpolygon');
%     if ~isempty(find_inpolygon),
%         if fktype==1 || fktype==3, 
%             IN = inpolygon(waven(nfhalf+1:nf,:),freq(nfhalf+1:nf,:),k,f);
%             dfk(nfhalf+1:nf,:) = dfk(nfhalf+1:nf,:).*IN;
%             IN = inpolygon(waven(1:nfhalf,:),freq(1:nfhalf,:),k,-f);
%             dfk(1:nfhalf,:) = dfk(1:nfhalf,:).*IN;
%         elseif fktype==2 || fktype==4,
%             IN = ~inpolygon(waven(nfhalf+1:nf,:),freq(nfhalf+1:nf,:),k,f);
%             dfk(nfhalf+1:nf,:) = dfk(nfhalf+1:nf,:).*IN;
%             IN = ~inpolygon(waven(1:nfhalf,:),freq(1:nfhalf,:),k,-f);
%             dfk(1:nfhalf,:) = dfk(1:nfhalf,:).*IN;
%         end;
%     elseif isempty(find_inpolygon),
%         disp('F-K FILTER > INPOLYGON not found! Working with ISINPOLY');
%         if fktype==1 || fktype==3, 
%             IN = isinpoly(waven(nfhalf+1:nf,:),freq(nfhalf+1:nf,:),k,f);
%             dfk(nfhalf+1:nf,:) = dfk(nfhalf+1:nf,:).*IN;
%             IN = isinpoly(waven(1:nfhalf,:),freq(1:nfhalf,:),k,-f);
%             dfk(1:nfhalf,:) = dfk(1:nfhalf,:).*IN;
%         elseif fktype==2 || fktype==4,
%             IN = ~isinpoly(waven(nfhalf+1:nf,:),freq(nfhalf+1:nf,:),k,f);
%             dfk(nfhalf+1:nf,:) = dfk(nfhalf+1:nf,:).*IN;
%             IN = ~isinpoly(waven(1:nfhalf,:),freq(1:nfhalf,:),k,-f);
%             dfk(1:nfhalf,:) = dfk(1:nfhalf,:).*IN;
%         end;
%     end
%     clear waven freq IN             % free some memory
%     delete(ffk); clear ffk          % delete the figure
    
% case 5                              % up-dip filter
%     dfk=D; 
%     for i = 1:nfhalf
%         for j=1:nkhalf
%             dfk(i,nkhalf+j) = 0.0;
%             dfk(i+nfhalf,j) = 0.0;
%         end
%     end
    
% case 6                              % down-dip filter
    dfk=D; 
    for i = 1:nfhalf
        for j=1:nkhalf
            dfk(i,j) = 0.0;
            dfk(i+nfhalf,j+nkhalf) = 0.0;
        end
    end
%end                                 % end switch loop

%%%%% Back to time - space domain 
dfk     = fftshift(dfk);
dfk     = real(ifft2(dfk));
%%%%% Recover filtered array from zero-padded array
dfk = dfk(1:ns,1:ntr);
% -------------------------------------------------
%%%%%   Remove the message box
% if exist('msg') && ishandle(msg),
%     delete(msg)
% end
return

