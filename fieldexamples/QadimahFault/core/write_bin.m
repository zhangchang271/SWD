%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.

%  file:       write_bin
%  purpose:    write to binary file 
%  assumption: None
%  input:      filename   ------ filename for reading
%              data       ------ data for writing

%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Oct 21, 2011

function write_bin(filename,data)

fid=fopen(filename,'wb');
fwrite(fid,data,'float');fclose(fid);

end
