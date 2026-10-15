% SWD Run the physical-distance SWD workflow.
%
% The implementation lives in SWD_distance.m.  This entry point keeps the
% familiar SWD.m command while making the distance-based workflow the default.

run(fullfile(fileparts(mfilename('fullpath')),'SWD_distance.m'));
