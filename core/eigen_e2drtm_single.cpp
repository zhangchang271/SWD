// Reconstruct the source field while propagating the original adjoint field.
// Copyright (C) 2024 Zhang Chang; distributed under GNU GPL v3.
#include "elastic_fd.hpp"
using namespace swd;
void mexFunction(int nlhs,mxArray* out[],int nrhs,const mxArray* in[]) {
    try {
        if(nrhs!=13 || nlhs!=3) throw std::runtime_error("Use e2drtm_eigen to call this kernel.");
        Elastic e=model(in,false);geometry(e,in[12]);auto receivers=points(in[10],e);
        MatrixXf residual=matrix(in[8]),s=matrix(in[11]);
        if(residual.rows()!=e.nt||residual.cols()!=e.ng||s.size()<e.nt) throw std::runtime_error("Invalid residual or source dimensions.");
        const mxArray* tape=in[9];
        if(!mxIsStruct(tape)) throw std::runtime_error("Invalid reconstruction data.");
        MatrixXf p=matrix(mxGetField(tape,0,"source")),before=matrix(mxGetField(tape,0,"source_before"));
        if(p.size()!=3||before.size()!=e.nt) throw std::runtime_error("Invalid source history.");
        Point source(p(0)-1-(e.mapped?.5:0),p(1)-1,e.nzbc,e.nxbc,e.mapped?e.pad_top+1:0);int source_type=int(p(2));
        State v(e.nzbc,e.nxbc),adj(e.nzbc,e.nxbc);Strain strain(e.nzbc,e.nxbc);
        readState(v,mxGetField(tape,0,"last"));
        Boundary boundary(e,true,filename(mxGetField(tape,0,"boundary_file")),true);
        MatrixXf cl_img=MatrixXf::Zero(e.nz,e.nx),cm_img=cl_img,illum=cl_img;
        const int z=e.pad_top+1,x=e.nbc,nz=e.nz,nx=e.nx;
        for(int it=e.nt-1;it>=0;--it) {
            e.strain(v,strain);
            if(e.mapped) e.mappedAdjoint(adj,strain,receivers,residual,it,source,source_type,cl_img,cm_img,illum);
            else if(it<e.nt-1) {
                e.adjointVelocity(adj);
                for(int j=0;j<e.ng;++j) receivers[j].add(adj.ww,e.dt*residual(it,j));
                e.adjointStress(adj);
                if(it%e.nt_interval==0) {
                    auto fux=strain.fux.block(z,x,nz,nx),bwz=strain.bwz.block(z,x,nz,nx);
                    auto fuz=strain.fuz.block(z,x,nz,nx),bwx=strain.bwx.block(z,x,nz,nx);
                    auto xx=adj.xx.block(z,x,nz,nx),zz=adj.zz.block(z,x,nz,nx),xz=adj.xz.block(z,x,nz,nx);
                    cl_img+=(fux+bwz).cwiseProduct(xx+zz);
                    cm_img+=2*(xx.cwiseProduct(fux)+zz.cwiseProduct(bwz)+0.5*xz.cwiseProduct(fuz+bwx));
                    illum+=(fux+bwz).cwiseProduct(fux+bwz);
                }
            }
            reverse(e,v,strain,boundary,source,s.data()[it],before(it),source_type,it);
        }
        out[0]=array(cl_img);out[1]=array(cm_img);out[2]=array(illum);
    } catch(const std::exception& ex) {mexErrMsgIdAndTxt("SWD:Adjoint", "%s",ex.what());}
}
