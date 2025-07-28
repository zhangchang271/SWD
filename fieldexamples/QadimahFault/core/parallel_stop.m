function np=parallel_stop
%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.
%
%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Sep 26, 2012
%
%  PARALLEL_INIT: stop the parallel mode for matlab 
%
%  np = parallel_init
%
%  OUT  np : number of processors
%
%  See also parall
np=matlabpool('size');
if np>0
    matlabpool close;
    display('Exit Matl');
end

end