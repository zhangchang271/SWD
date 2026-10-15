function [cl_img,cm_img,illum_div]=e2drtm_eigen(seismo_w,forward)
% Reconstruct the source field using the matching forward model and settings.
cleanup=onCleanup(@()delete(forward.boundary_file));
[cl_img,cm_img,illum_div]=eigen_e2drtm_single(forward.input,forward.temp, ...
    forward.ca,forward.cl,forward.cm,forward.cm,forward.b,forward.b, ...
    seismo_w,forward,forward.grid.receivers,forward.s,forward.grid);
end
