function vmodel = generateCheckerboard(nz, nx, vlayer, hlayer, baseV, floatV,samepad)
    vmodel = baseV * ones(nz, nx);
    
    dz = floor(nz / vlayer);
    dx = floor(nx / hlayer);
    [X, Z] = meshgrid(1:nx, 1:nz);
    checker = mod(floor((X-1)/dx) + floor((Z-1)/dz), 2) == 0;    
    vmodel(checker) = baseV.*(1+floatV);
    vmodel(~checker) = baseV.*(1-floatV);
    if samepad==1
        vmodel(dz*vlayer+1:nz,:) = repmat(vmodel(dz*vlayer,:),length(dz*vlayer+1:nz),1);
    end
end
