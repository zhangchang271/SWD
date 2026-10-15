// Test-only full trajectory: compare every reconstructed field and gradient.
#include "../core/elastic_fd.hpp"
using namespace swd;
void mexFunction(int nlhs,mxArray* out[],int nrhs,const mxArray* in[]) {
    if(nrhs!=12||nlhs!=4) mexErrMsgIdAndTxt("SWD:Test","Expected forward inputs followed by residual.");
    Elastic e=model(in,true);geometry(e,in[10]);
    MatrixXf residual=matrix(in[11]);
    // dt is supplied separately in the unused forward parameter slot 7.
    e.dt=float(number(in[0],7));
    auto receivers=points(in[9],e);
    Point source(number(in[0],5)-1-(e.mapped?.5:0),number(in[0],6)-1,e.nzbc,e.nxbc,e.mapped?e.pad_top+1:0);
    int type=int(number(in[0],10));MatrixXf s=matrix(in[8]);
    State v(e.nzbc,e.nxbc),adj(e.nzbc,e.nxbc);Strain strain(e.nzbc,e.nxbc);
    Boundary boundary(e);std::vector<State> states;std::vector<Strain> strains;VectorXf before(e.nt);
    for(int it=0;it<e.nt;++it) {
        e.velocity(v);before(it)=e.source(v,source,s.data()[it],type);
        e.strain(v,strain);e.stress(v,strain);boundary.save(v,it);
        states.push_back(v);strains.push_back(strain);
    }
    double err=0,norm=0,erre=0,norme=0;
    MatrixXf gl=MatrixXf::Zero(e.nz,e.nx),gm=gl,illum=gl;
    int z=e.pad_top+1,x=e.nbc;
    for(int it=e.nt-1;it>=0;--it) {
        e.strain(v,strain);
        auto f=v.fields(),ref=states[it].fields();
        for(int j=0;j<5;++j) {
            err+=((*f[j]-*ref[j]).block(z,x,e.nz,e.nx)).cast<double>().squaredNorm();
            norm+=ref[j]->block(z,x,e.nz,e.nx).cast<double>().squaredNorm();
        }
        for(auto pair:{std::make_pair(&strain.fux,&strains[it].fux),std::make_pair(&strain.fuz,&strains[it].fuz),std::make_pair(&strain.bwx,&strains[it].bwx),std::make_pair(&strain.bwz,&strains[it].bwz)}) {
            erre+=(*pair.first-*pair.second).block(z,x,e.nz,e.nx).cast<double>().squaredNorm();
            norme+=pair.second->block(z,x,e.nz,e.nx).cast<double>().squaredNorm();
        }
        if(e.mapped) e.mappedAdjoint(adj,strains[it],receivers,residual,it,source,type,gl,gm,illum);
        else if(it<e.nt-1) {
            e.adjointVelocity(adj);for(int j=0;j<e.ng;++j) receivers[j].add(adj.ww,e.dt*residual(it,j));e.adjointStress(adj);
            if(it%e.nt_interval==0) {
                auto ex=strains[it].fux.block(z,x,e.nz,e.nx),ez=strains[it].bwz.block(z,x,e.nz,e.nx);
                auto g=strains[it].fuz.block(z,x,e.nz,e.nx)+strains[it].bwx.block(z,x,e.nz,e.nx);
                auto ax=adj.xx.block(z,x,e.nz,e.nx),az=adj.zz.block(z,x,e.nz,e.nx),at=adj.xz.block(z,x,e.nz,e.nx);
                gl+=(ex+ez).cwiseProduct(ax+az);gm+=2*(ax.cwiseProduct(ex)+az.cwiseProduct(ez)+.5*at.cwiseProduct(g));
                illum+=(ex+ez).cwiseProduct(ex+ez);
            }
        }
        reverse(e,v,strain,boundary,source,s.data()[it],before(it),type,it);
    }
    MatrixXf errors(2,1);errors<<std::sqrt(err/norm),std::sqrt(erre/norme);
    out[0]=array(errors);out[1]=array(gl);out[2]=array(gm);out[3]=array(illum);
}
