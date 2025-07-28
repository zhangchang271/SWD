function np=parallel_init(nn)
%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.
%
%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Sep 26, 2012
%
%  PARALLEL_INIT: initialize the parallel mode for matlab 
%
%  np = parallel_init
%
%  OUT  np : number of processors
%
%  See also parall
v = version('-release');

 if strcmp(v(1:4),'2012')
     np=parpool('size');
 else
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        np = 0;
    else
        np = poolobj.NumWorkers;
    end    
 end
 display(['Checking......', num2str(nn),' CPUs are currently using!']);
 
    
if np==0
    if strcmp(v(1:4),'2012')
        matlabpool('local',nn);
    else
        parpool('local',nn);
    end
    
%     np=matlabpool('size');
end
    display([num2str(nn),' CPUs are using now!']);
    

    
    
    
    
    
    

end