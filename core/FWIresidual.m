function [seismo_v_d1,res_r] = FWIresidual(seismo_v,seismo_v_d)
seismo_v_d1 = seismo_v - seismo_v_d;
res_r = seismo_v - seismo_v_d;