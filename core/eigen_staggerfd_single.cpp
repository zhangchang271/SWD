// Elastic forward propagation; terminal state and boundary saving for reverse reconstruction.
// Copyright (C) 2024 Zhang Chang; distributed under GNU GPL v3.
#include "elastic_fd.hpp"
using namespace swd;
void mexFunction(int nlhs,mxArray* out[],int nrhs,const mxArray* in[]) {
    try {
        if(nrhs!=12 || nlhs<2 || nlhs>3) throw std::runtime_error("Use staggerfd_eigen to call this kernel.");
        Elastic e=model(in,true);geometry(e,in[10]);
        auto receivers=points(in[9],e);auto receivers_u=points(in[9],e,false);
        Point source(number(in[0],5)-1-(e.mapped?.5:0),number(in[0],6)-1,e.nzbc,e.nxbc,e.mapped?e.pad_top+1:0);
        int source_type=int(number(in[0],10));
        MatrixXf s=matrix(in[8]);
        if(s.size()<e.nt) throw std::runtime_error("Source wavelet is shorter than nt.");
        State v(e.nzbc,e.nxbc);Strain strain(e.nzbc,e.nxbc);
        std::string path=nlhs==3?filename(in[11]):"";
        Boundary boundary(e,nlhs==3,path);
        MatrixXf before(nlhs==3?e.nt:0,1),u(e.nt,e.ng),w(e.nt,e.ng);
        for(int it=0;it<e.nt;++it) {
            e.velocity(v);
            float old=e.source(v,source,s.data()[it],source_type);
            e.strain(v,strain);e.stress(v,strain);
            for(int j=0;j<e.ng;++j) {u(it,j)=receivers_u[j].sample(v.uu);w(it,j)=receivers[j].sample(v.ww);}
            if(nlhs==3) {boundary.save(v,it);before(it)=old;}
        }
        out[0]=array(u);out[1]=array(w);
        if(nlhs==3) {
            const char* names[]={"last","boundary_file","source_before","source"};
            out[2]=mxCreateStructMatrix(1,1,4,names);
            mxSetField(out[2],0,"last",stateArray(v));
            mxSetField(out[2],0,"boundary_file",mxCreateString(path.c_str()));
            mxSetField(out[2],0,"source_before",array(before));
            MatrixXf p(3,1);p<<number(in[0],5),number(in[0],6),source_type;
            mxSetField(out[2],0,"source",array(p));
        }
    } catch(const std::exception& ex) {mexErrMsgIdAndTxt("SWD:Forward", "%s",ex.what());}
}
