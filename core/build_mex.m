function build_mex
% Compile the two active kernels once before starting parallel workers.
core=fileparts(mfilename('fullpath'));
header=dir(fullfile(core,'elastic_fd.hpp'));
names={'eigen_staggerfd_single','eigen_e2drtm_single'};
include='/usr/include/eigen3';
if ~isfile(fullfile(include,'Eigen','Core'))
    include=fullfile(core,'eigen3');
end
for j=1:numel(names)
    source=fullfile(core,[names{j},'.cpp']); binary=fullfile(core,[names{j},'.',mexext]);
    src=dir(source); bin=dir(binary);
    if isempty(bin)||bin.datenum<max(src.datenum,header.datenum)
        assert(isfile(fullfile(include,'Eigen','Core')),'SWD:Eigen', ...
            'Install Eigen 3 headers or place them in core/eigen3.');
        % R2021a on this host supplies an older C++ runtime than the system g++.
        compiler={};
        if isunix && ~ismac && isfile('/usr/bin/g++-9')
            compiler={'CXX=/usr/bin/g++-9'};
        end
        clear(names{j});
        mex('-R2018a',['-I',include],'-DNDEBUG',compiler{:},source,'-outdir',core);
    end
end
end
