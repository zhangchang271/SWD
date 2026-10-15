// Elastic finite differences shared by forward propagation and reconstruction.
// Copyright (C) 2024 Zhang Chang; distributed under GNU GPL v3.
#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <cstdio>
#include <string>
using Eigen::MatrixXf;
using Eigen::VectorXf;

namespace swd {
struct State {
    MatrixXf uu,ww,xx,zz,xz;
    State(int rows,int cols): uu(MatrixXf::Zero(rows,cols)),ww(uu),xx(uu),zz(uu),xz(uu) {}
    std::array<MatrixXf*,5> fields() {return {{&uu,&ww,&xx,&zz,&xz}};}
};
struct Strain {
    MatrixXf fux,fuz,bwx,bwz;
    Strain(int rows,int cols): fux(MatrixXf::Zero(rows,cols)),fuz(fux),bwx(fux),bwz(fux) {}
};
struct Point {
    std::array<int,4> index;
    std::array<float,4> weight;
    Point(double z,double x,int rows,int cols,int first=0) {
        int iz=std::max(first,int(std::floor(z))),ix=int(std::floor(x));
        if (iz<0 || iz+1>=rows || ix<0 || ix+1>=cols) throw std::runtime_error("Position is outside the grid.");
        float dz=float(z-iz),dx=float(x-ix);
        index={{iz+rows*ix,iz+1+rows*ix,iz+rows*(ix+1),iz+1+rows*(ix+1)}};
        weight={{(1-dz)*(1-dx),dz*(1-dx),(1-dz)*dx,dz*dx}};
    }
    float sample(const MatrixXf& a) const {float v=0; for(int j=0;j<4;++j) if(weight[j]!=0) v+=weight[j]*a.data()[index[j]]; return v;}
    void add(MatrixXf& a,float v,int shift=0) const {for(int j=0;j<4;++j) if(weight[j]!=0) a.data()[index[j]+shift]+=weight[j]*v;}
    void prescribe(MatrixXf& a,float v) const {
        for(int j=0;j<4;++j) if(weight[j]==1) {a.data()[index[j]]=v;return;}
        float n=0; for(float w:weight) n+=w*w;
        add(a,(v-sample(a))/n);
    }
};
struct Elastic {
    int nt,nzbc,nxbc,ng,nt_interval,nz,nx,nbc,fd_order_num,pad_top,k,i;
    float dtx,dt;
    MatrixXf temp,ca,cl,cm,cm1,b,b1;
    bool mapped=false;
    Eigen::SparseMatrix<float> G;
    MatrixXf jn,js,invU,invW;
    VectorXf packVelocity(const State& v) const {
        VectorXf q(2*nzbc*nxbc);q<<Eigen::Map<const VectorXf>(v.uu.data(),v.uu.size()),Eigen::Map<const VectorXf>(v.ww.data(),v.ww.size());return q;
    }
    VectorXf packStress(const State& v,bool weighted) const {
        int N=nzbc*nxbc;VectorXf q(3*N);
        for(int n=0;n<N;++n) {
            q(n)=v.xx.data()[n]*(weighted?jn.data()[n]:1);
            q(n+N)=v.zz.data()[n]*(weighted?jn.data()[n]:1);
            q(n+2*N)=v.xz.data()[n]*(weighted?js.data()[n]:1);
        }return q;
    }
    void mappedVelocity(State& v) const {
        VectorXf q=G.transpose()*packStress(v,true);int N=nzbc*nxbc;
        for(int n=0;n<N;++n) {
            v.uu.data()[n]=temp.data()[n]*v.uu.data()[n]-dtx*invU.data()[n]*q(n);
            v.ww.data()[n]=temp.data()[n]*v.ww.data()[n]-dtx*invW.data()[n]*q(n+N);
        }
    }
    void mappedStrain(const State& v,Strain& strain) const {
        VectorXf q=G*packVelocity(v);int N=nzbc*nxbc;
        for(int n=0;n<N;++n) {strain.fux.data()[n]=q(n);strain.bwz.data()[n]=q(n+N);strain.fuz.data()[n]=q(n+2*N);}
    }
    void constitutive(const Strain& e,State& stress) const {
        stress.xx=ca.cwiseProduct(e.fux)+cl.cwiseProduct(e.bwz);
        stress.zz=ca.cwiseProduct(e.bwz)+cl.cwiseProduct(e.fux);
        stress.xz=cm1.cwiseProduct(e.fuz+e.bwx);
        if(mapped) {
            int z=pad_top+1;
            stress.xx.row(z)=(ca.row(z)-cl.row(z).cwiseProduct(cl.row(z)).cwiseQuotient(ca.row(z))).cwiseProduct(e.fux.row(z));
            stress.zz.row(z).setZero();
        }
    }
    void mappedAdjoint(State& adj,const Strain& strain,const std::vector<Point>& receivers,
                       const MatrixXf& residual,int it,const Point& source,int type,
                       MatrixXf& gl,MatrixXf& gm,MatrixXf& illum) const {
        int N=nzbc*nxbc,z0=pad_top+1;
        // Exact transpose of the leapfrog step, with transpose receiver injection.
        if(it%nt_interval==0) for(int x=0;x<nx;++x) for(int z=0;z<nz;++z) {
            int n=z+z0+nzbc*(x+nbc);float ex=strain.fux.data()[n],ez=strain.bwz.data()[n],g=strain.fuz.data()[n];
            float ax=adj.xx.data()[n],az=adj.zz.data()[n],at=adj.xz.data()[n];
            if(z==0) {
                float ratio=cl.data()[n]/ca.data()[n];
                gl(z,x)+=dtx*(1-ratio)*(1-ratio)*ax*ex;
                gm(z,x)+=dtx*(2+2*ratio*ratio)*ax*ex+dtx*at*g;
                illum(z,x)+=(1-ratio)*(1-ratio)*ex*ex;
            } else {
                gl(z,x)+=dtx*(ex+ez)*(ax+az);
                gm(z,x)+=dtx*(2*ax*ex+2*az*ez+at*g);
                illum(z,x)+=(ex+ez)*(ex+ez);
            }
        }
        Strain as(nzbc,nxbc);as.fux=adj.xx;as.bwz=adj.zz;as.fuz=adj.xz;
        State c(nzbc,nxbc);constitutive(as,c);
        VectorXf q=G.transpose()*packStress(c,false);
        for(int n=0;n<N;++n) {adj.uu.data()[n]+=dtx*q(n);adj.ww.data()[n]+=dtx*q(n+N);}
        for(int j=0;j<ng;++j) receivers[j].add(adj.ww,dt*residual(it,j));
        if(type==5) source.prescribe(adj.ww,0);
        State a(nzbc,nxbc);a.uu=adj.uu.cwiseProduct(invU);a.ww=adj.ww.cwiseProduct(invW);
        q=G*packVelocity(a);
        for(int n=0;n<N;++n) {
            adj.xx.data()[n]=temp.data()[n]*adj.xx.data()[n]-dtx*jn.data()[n]*q(n);
            adj.zz.data()[n]=temp.data()[n]*adj.zz.data()[n]-dtx*jn.data()[n]*q(n+N);
            adj.xz.data()[n]=temp.data()[n]*adj.xz.data()[n]-dtx*js.data()[n]*q(n+2*N);
        }
        adj.uu=temp.cwiseProduct(adj.uu);adj.ww=temp.cwiseProduct(adj.ww);
    }

    const float S21=1, S41=1.1250f,S42=-0.0416666667f;
    const float S61=1.17187f,S62=-6.51042E-2f,S63=4.68750E-3f;
    const float S81=1.19629f,S82=-7.97526E-2f,S83=9.57031E-3f,S84=-6.97545E-4f;
    void dimensions() {
        pad_top=(fd_order_num-20)/2;
        if(pad_top<1 || pad_top>4 || fd_order_num%2) throw std::runtime_error("Unsupported finite-difference order.");
        k=nzbc-2*pad_top; i=nxbc-2*pad_top;
    }
    void velocity(State& v) const {
        if(mapped) {mappedVelocity(v);return;}
        auto &uu=v.uu,&ww=v.ww,&xx=v.xx,&zz=v.zz,&xz=v.xz;
       if (fd_order_num == 22){
         uu.block(1,1,k,i) = temp.block(1,1,k,i).cwiseProduct(uu.block(1,1,k,i)) + b.block(1,1,k,i).cwiseProduct(
                 xx.block(1,1+1,k,i) - xx.block(1,1,k,i) + xz.block(1,1,k,i) - xz.block(1-1,1,k,i));
         ww.block(1,1,k,i) = temp.block(1,1,k,i).cwiseProduct(ww.block(1,1,k,i)) + b1.block(1,1,k,i).cwiseProduct(
                 xz.block(1,1,k,i) - xz.block(1,1-1,k,i) + zz.block(1+1,1,k,i) - zz.block(1,1,k,i));}
       else if(fd_order_num == 24){
         uu.block(2,2,k,i) = temp.block(2,2,k,i).cwiseProduct(uu.block(2,2,k,i)) + b.block(2,2,k,i).cwiseProduct(
                 S41*(xx.block(2,2,k,i) - xx.block(2,2-1,k,i)) + S42*(xx.block(2,2+1,k,i) - xx.block(2,2-2,k,i)) +
                 S41*(xz.block(2,2,k,i) - xz.block(2-1,2,k,i)) + S42*(xz.block(2+1,2,k,i) - xz.block(2-2,2,k,i)));
         ww.block(2,2,k,i) = temp.block(2,2,k,i).cwiseProduct(ww.block(2,2,k,i)) + b1.block(2,2,k,i).cwiseProduct(
                 S41*(xz.block(2,2+1,k,i) - xz.block(2,2,k,i)) + S42*(xz.block(2,2+2,k,i) - xz.block(2,2-1,k,i)) +
                 S41*(zz.block(2+1,2,k,i) - zz.block(2,2,k,i)) + S42*(zz.block(2+2,2,k,i) - zz.block(2-1,2,k,i)));}
       else if(fd_order_num == 26){
         uu.block(3,3,k,i) = temp.block(3,3,k,i).cwiseProduct(uu.block(3,3,k,i)) + b.block(3,3,k,i).cwiseProduct(
                 S61*(xx.block(3,3,k,i) - xx.block(3,3-1,k,i)) + S62*(xx.block(3,3+1,k,i) - xx.block(3,3-2,k,i)) +
                 S63*(xx.block(3,3+2,k,i) - xx.block(3,3-3,k,i)) + S61*(xz.block(3,3,k,i) - xz.block(3-1,3,k,i)) +
                 S62*(xz.block(3+1,3,k,i) - xz.block(3-2,3,k,i)) + S63*(xz.block(3+2,3,k,i) - xz.block(3-3,3,k,i)));
         ww.block(3,3,k,i) = temp.block(3,3,k,i).cwiseProduct(ww.block(3,3,k,i)) + b1.block(3,3,k,i).cwiseProduct(
                 S61*(xz.block(3,3+1,k,i) - xz.block(3,3,k,i)) + S62*(xz.block(3,3+2,k,i) - xz.block(3,3-1,k,i)) +
                 S63*(xz.block(3,3+3,k,i) - xz.block(3,3-2,k,i)) + S61*(zz.block(3+1,3,k,i) - zz.block(3,3,k,i)) +
                 S62*(zz.block(3+2,3,k,i) - zz.block(3-1,3,k,i)) + S63*(zz.block(3+3,3,k,i) - zz.block(3-2,3,k,i)));}
       else if(fd_order_num == 28){
         uu.block(4,4,k,i) = temp.block(4,4,k,i).cwiseProduct(uu.block(4,4,k,i)) + b.block(4,4,k,i).cwiseProduct(
                 S81*(xx.block(4,4,k,i) - xx.block(4,4-1,k,i)) + S82*(xx.block(4,4+1,k,i) - xx.block(4,4-2,k,i)) +
                 S83*(xx.block(4,4+2,k,i) - xx.block(4,4-3,k,i)) + S84*(xx.block(4,4+3,k,i) - xx.block(4,4-4,k,i)) +
                 S81*(xz.block(4,4,k,i) - xz.block(4-1,4,k,i)) + S82*(xz.block(4+1,4,k,i) - xz.block(4-2,4,k,i)) +
                 S83*(xz.block(4+2,4,k,i) - xz.block(4-3,4,k,i)) + S84*(xz.block(4+3,4,k,i) - xz.block(4-4,4,k,i)));
         ww.block(4,4,k,i) = temp.block(4,4,k,i).cwiseProduct(ww.block(4,4,k,i)) + b1.block(4,4,k,i).cwiseProduct(
                 S81*(xz.block(4,4+1,k,i) - xz.block(4,4,k,i)) + S82*(xz.block(4,4+2,k,i) - xz.block(4,4-1,k,i)) +
                 S83*(xz.block(4,4+3,k,i) - xz.block(4,4-2,k,i)) + S84*(xz.block(4,4+4,k,i) - xz.block(4,4-3,k,i)) +
                 S81*(zz.block(4+1,4,k,i) - zz.block(4,4,k,i)) + S82*(zz.block(4+2,4,k,i) - zz.block(4-1,4,k,i)) +
                 S83*(zz.block(4+3,4,k,i) - zz.block(4-2,4,k,i)) + S84*(zz.block(4+4,4,k,i) - zz.block(4-3,4,k,i)));}
    }
    void strain(const State& v,Strain& e) const {
        if(mapped) {mappedStrain(v,e);return;}
        const auto &uu=v.uu,&ww=v.ww;
        auto &fux=e.fux,&fuz=e.fuz,&bwx=e.bwx,&bwz=e.bwz;
       if(fd_order_num == 22){
            fux.block(1,1,k,i) = uu.block(1,1,k,i) - uu.block(1,1-1,k,i);
            fuz.block(1,1,k,i) = uu.block(1+1,1,k,i) - uu.block(1,1,k,i);
            bwx.block(1,1,k,i) = ww.block(1,1+1,k,i) - ww.block(1,1,k,i);
            bwz.block(1,1,k,i) = ww.block(1,1,k,i) - ww.block(1-1,1,k,i);}
       else if(fd_order_num == 24){
            fux.block(2,2,k,i) = S41*(uu.block(2,2+1,k,i) - uu.block(2,2,k,i)) + S42*(uu.block(2,2+2,k,i) - uu.block(2,2-1,k,i));
            fuz.block(2,2,k,i) = S41*(uu.block(2+1,2,k,i) - uu.block(2,2,k,i)) + S42*(uu.block(2+2,2,k,i) - uu.block(2-1,2,k,i));
            bwx.block(2,2,k,i) = S41*(ww.block(2,2,k,i) - ww.block(2,2-1,k,i)) + S42*(ww.block(2,2+1,k,i) - ww.block(2,2-2,k,i));
            bwz.block(2,2,k,i) = S41*(ww.block(2,2,k,i) - ww.block(2-1,2,k,i)) + S42*(ww.block(2+1,2,k,i) - ww.block(2-2,2,k,i));}
       else if(fd_order_num == 26){
        fux.block(3,3,k,i) = S61*(uu.block(3,3+1,k,i) - uu.block(3,3,k,i)) + S62*(uu.block(3,3+2,k,i) - uu.block(3,3-1,k,i))+
            S63*(uu.block(3,3+3,k,i) - uu.block(3,3-2,k,i));
        fuz.block(3,3,k,i) = S61*(uu.block(3+1,3,k,i) - uu.block(3,3,k,i)) + S62*(uu.block(3+2,3,k,i) - uu.block(3-1,3,k,i))+
            S63*(uu.block(3+3,3,k,i) - uu.block(3-2,3,k,i));
        bwx.block(3,3,k,i) = S61*(ww.block(3,3,k,i) - ww.block(3,3-1,k,i)) + S62*(ww.block(3,3+1,k,i) - ww.block(3,3-2,k,i))+
            S63*(ww.block(3,3+2,k,i) - ww.block(3,3-3,k,i));
        bwz.block(3,3,k,i) = S61*(ww.block(3,3,k,i) - ww.block(3-1,3,k,i)) + S62*(ww.block(3+1,3,k,i) - ww.block(3-2,3,k,i))+
            S63*(ww.block(3+2,3,k,i) - ww.block(3-3,3,k,i));}
      else if(fd_order_num == 28){
        fux.block(4,4,k,i) = S81*(uu.block(4,4+1,k,i) - uu.block(4,4,k,i)) + S82*(uu.block(4,4+2,k,i) - uu.block(4,4-1,k,i))+
            S83*(uu.block(4,4+3,k,i) - uu.block(4,4-2,k,i)) + S84*(uu.block(4,4+4,k,i) - uu.block(4,4-3,k,i));
        fuz.block(4,4,k,i) = S81*(uu.block(4+1,4,k,i) - uu.block(4,4,k,i)) + S82*(uu.block(4+2,4,k,i) - uu.block(4-1,4,k,i))+
            S83*(uu.block(4+3,4,k,i) - uu.block(4-2,4,k,i)) + S84*(uu.block(4+4,4,k,i) - uu.block(4-3,4,k,i));
        bwx.block(4,4,k,i) = S81*(ww.block(4,4,k,i) - ww.block(4,4-1,k,i)) + S82*(ww.block(4,4+1,k,i) - ww.block(4,4-2,k,i))+
            S83*(ww.block(4,4+2,k,i) - ww.block(4,4-3,k,i)) + S84*(ww.block(4,4+3,k,i) - ww.block(4,4-4,k,i));
        bwz.block(4,4,k,i) = S81*(ww.block(4,4,k,i) - ww.block(4-1,4,k,i)) + S82*(ww.block(4+1,4,k,i) - ww.block(4-2,4,k,i))+
            S83*(ww.block(4+2,4,k,i) - ww.block(4-3,4,k,i)) + S84*(ww.block(4+3,4,k,i) - ww.block(4-4,4,k,i));}
    }
    void stress(State& v,const Strain& e) const {
        if(mapped) {
            State q(nzbc,nxbc);constitutive(e,q);
            v.xx=temp.cwiseProduct(v.xx)+q.xx*dtx;
            v.zz=temp.cwiseProduct(v.zz)+q.zz*dtx;
            v.xz=temp.cwiseProduct(v.xz)+q.xz*dtx;return;
        }
        v.xx=temp.cwiseProduct(v.xx)+(ca.cwiseProduct(e.fux)+cl.cwiseProduct(e.bwz))*dtx;
        v.zz=temp.cwiseProduct(v.zz)+(ca.cwiseProduct(e.bwz)+cl.cwiseProduct(e.fux))*dtx;
        v.xz=temp.cwiseProduct(v.xz)+cm1.cwiseProduct(e.fuz+e.bwx)*dtx;
        v.zz.row(pad_top).setZero();
    }
    float source(State& v,const Point& p,float s,int type) const {
        float old=p.sample(v.ww);
        if(type==1) {p.add(v.xx,s);p.add(v.zz,s);}
        else if(type==2) {p.add(v.uu,s);p.add(v.uu,-s,-1);p.add(v.ww,s);p.add(v.ww,s,nzbc);}
        else if(type==3) {p.add(v.uu,s);p.add(v.ww,s);}
        else if(type==4) p.add(v.zz,s);
        else if(type==5) p.prescribe(v.ww,s);
        return old;
    }
    void adjointVelocity(State& v) const {
        auto &uu_b=v.uu,&ww_b=v.ww,&xx_b=v.xx,&zz_b=v.zz,&xz_b=v.xz;
    if(fd_order_num==22){
        uu_b.block(1,1,k,i) = temp.block(1,1,k,i).cwiseProduct(uu_b.block(1,1,k,i)) - b.block(1,1,k,i).cwiseProduct(S21*(ca.block(1,1,k,i).cwiseProduct(xx_b.block(1,1,k,i)) - ca.block(1,1-1,k,i).cwiseProduct(xx_b.block(1,1-1,k,i))) +
            S21*(cl.block(1,1,k,i).cwiseProduct(zz_b.block(1,1,k,i))-cl.block(1,1-1,k,i).cwiseProduct(zz_b.block(1,1-1,k,i))) + S21*(cm1.block(1,1,k,i).cwiseProduct(xz_b.block(1,1,k,i))-cm1.block(1-1,1,k,i).cwiseProduct(xz_b.block(1-1,1,k,i))));

        ww_b.block(1,1,k,i) = temp.block(1,1,k,i).cwiseProduct(ww_b.block(1,1,k,i)) - b1.block(1,1,k,i).cwiseProduct( S21*(cl.block(1+1,1,k,i).cwiseProduct(xx_b.block(1+1,1,k,i)) - cl.block(1,1,k,i).cwiseProduct(xx_b.block(1,1,k,i))) +
            S21*(ca.block(1+1,1,k,i).cwiseProduct(zz_b.block(1+1,1,k,i)) - ca.block(1,1,k,i).cwiseProduct(zz_b.block(1,1,k,i))) + S21*(cm1.block(1,1+1,k,i).cwiseProduct(xz_b.block(1,1+1,k,i)) - cm1.block(1,1,k,i).cwiseProduct(xz_b.block(1,1,k,i))));}
    else if(fd_order_num == 24){
    uu_b.block(2,2,k,i) = temp.block(2,2,k,i).cwiseProduct(uu_b.block(2,2,k,i)) - b.block(2,2,k,i).cwiseProduct(S41*(ca.block(2,2,k,i).cwiseProduct(xx_b.block(2,2,k,i)) - ca.block(2,2-1,k,i).cwiseProduct(xx_b.block(2,2-1,k,i))) +
       S42*(ca.block(2,2+1,k,i).cwiseProduct(xx_b.block(2,2+1,k,i)) - ca.block(2,2-2,k,i).cwiseProduct(xx_b.block(2,2-2,k,i))) +
       S41*(cl.block(2,2,k,i).cwiseProduct(zz_b.block(2,2,k,i)) - cl.block(2,2-1,k,i).cwiseProduct(zz_b.block(2,2-1,k,i))) + S42*(cl.block(2,2+1,k,i).cwiseProduct(zz_b.block(2,2+1,k,i)) - cl.block(2,2-2,k,i).cwiseProduct(zz_b.block(2,2-2,k,i))) +
       S41*(cm1.block(2,2,k,i).cwiseProduct(xz_b.block(2,2,k,i)) - cm1.block(2-1,2,k,i).cwiseProduct(xz_b.block(2-1,2,k,i))) + S42*(cm1.block(2+1,2,k,i).cwiseProduct(xz_b.block(2+1,2,k,i)) - cm1.block(2-2,2,k,i).cwiseProduct(xz_b.block(2-2,2,k,i))));

    ww_b.block(2,2,k,i) = temp.block(2,2,k,i).cwiseProduct(ww_b.block(2,2,k,i)) - b1.block(2,2,k,i).cwiseProduct(S41*(cl.block(2+1,2,k,i).cwiseProduct(xx_b.block(2+1,2,k,i)) - cl.block(2,2,k,i).cwiseProduct(xx_b.block(2,2,k,i))) +
       S42*(cl.block(2+2,2,k,i).cwiseProduct(xx_b.block(2+2,2,k,i)) - cl.block(2-1,2,k,i).cwiseProduct(xx_b.block(2-1,2,k,i))) +
       S41*(ca.block(2+1,2,k,i).cwiseProduct(zz_b.block(2+1,2,k,i)) - ca.block(2,2,k,i).cwiseProduct(zz_b.block(2,2,k,i))) + S42*(ca.block(2+2,2,k,i).cwiseProduct(zz_b.block(2+2,2,k,i)) - ca.block(2-1,2,k,i).cwiseProduct(zz_b.block(2-1,2,k,i))) +
       S41*(cm1.block(2,2+1,k,i).cwiseProduct(xz_b.block(2,2+1,k,i)) - cm1.block(2,2,k,i).cwiseProduct(xz_b.block(2,2,k,i))) + S42*(cm1.block(2,2+2,k,i).cwiseProduct(xz_b.block(2,2+2,k,i)) - cm1.block(2,2-1,k,i).cwiseProduct(xz_b.block(2,2-1,k,i))));
    }
    else if(fd_order_num == 26){
        uu_b.block(3,3,k,i) = temp.block(3,3,k,i).cwiseProduct(uu_b.block(3,3,k,i)) - b.block(3,3,k,i).cwiseProduct(S61*(ca.block(3,3,k,i).cwiseProduct(xx_b.block(3,3,k,i)) - ca.block(3,3-1,k,i).cwiseProduct(xx_b.block(3,3-1,k,i))) +
                S62*(ca.block(3,3+1,k,i).cwiseProduct(xx_b.block(3,3+1,k,i)) - ca.block(3,3-2,k,i).cwiseProduct(xx_b.block(3,3-2,k,i))) +
            S63*(ca.block(3,3+2,k,i).cwiseProduct(xx_b.block(3,3+2,k,i)) - ca.block(3,3-3,k,i).cwiseProduct(xx_b.block(3,3-3,k,i))) +
            S61*(cl.block(3,3,k,i).cwiseProduct(zz_b.block(3,3,k,i)) - cl.block(3,3-1,k,i).cwiseProduct(zz_b.block(3,3-1,k,i))) + S62*(cl.block(3,3+1,k,i).cwiseProduct(zz_b.block(3,3+1,k,i)) - cl.block(3,3-2,k,i).cwiseProduct(zz_b.block(3,3-2,k,i))) +
                S63*(cl.block(3,3+2,k,i).cwiseProduct(zz_b.block(3,3+2,k,i)) - cl.block(3,3-3,k,i).cwiseProduct(zz_b.block(3,3-3,k,i))) +
            S61*( cm1.block(3,3,k,i).cwiseProduct(xz_b.block(3,3,k,i)) - cm1.block(3-1,3,k,i).cwiseProduct(xz_b.block(3-1,3,k,i))) + S62*( cm1.block(3+1,3,k,i).cwiseProduct(xz_b.block(3+1,3,k,i)) - cm1.block(3-2,3,k,i).cwiseProduct(xz_b.block(3-2,3,k,i))) +
                S63*( cm1.block(3+2,3,k,i).cwiseProduct(xz_b.block(3+2,3,k,i)) - cm1.block(3-3,3,k,i).cwiseProduct(xz_b.block(3-3,3,k,i))));

        ww_b.block(3,3,k,i) = temp.block(3,3,k,i).cwiseProduct(ww_b.block(3,3,k,i)) - b1.block(3,3,k,i).cwiseProduct(S61*(cl.block(3+1,3,k,i).cwiseProduct(xx_b.block(3+1,3,k,i)) - cl.block(3,3,k,i).cwiseProduct(xx_b.block(3,3,k,i))) +
                S62*(cl.block(3+2,3,k,i).cwiseProduct(xx_b.block(3+2,3,k,i)) - cl.block(3-1,3,k,i).cwiseProduct(xx_b.block(3-1,3,k,i))) +
            S63*(cl.block(3+3,3,k,i).cwiseProduct(xx_b.block(3+3,3,k,i)) - cl.block(3-2,3,k,i).cwiseProduct(xx_b.block(3-2,3,k,i))) +
            S61*(ca.block(3+1,3,k,i).cwiseProduct(zz_b.block(3+1,3,k,i)) - ca.block(3,3,k,i).cwiseProduct(zz_b.block(3,3,k,i))) + S62*(ca.block(3+2,3,k,i).cwiseProduct(zz_b.block(3+2,3,k,i)) - ca.block(3-1,3,k,i).cwiseProduct(zz_b.block(3-1,3,k,i))) +
                S63*(ca.block(3+3,3,k,i).cwiseProduct(zz_b.block(3+3,3,k,i)) - ca.block(3-2,3,k,i).cwiseProduct(zz_b.block(3-2,3,k,i))) +
            S61*(cm1.block(3,3+1,k,i).cwiseProduct(xz_b.block(3,3+1,k,i)) - cm1.block(3,3,k,i).cwiseProduct(xz_b.block(3,3,k,i))) + S62*(cm1.block(3,3+2,k,i).cwiseProduct(xz_b.block(3,3+2,k,i)) - cm1.block(3,3-1,k,i).cwiseProduct(xz_b.block(3,3-1,k,i))) +
                S63*(cm1.block(3,3+3,k,i).cwiseProduct(xz_b.block(3,3+3,k,i)) - cm1.block(3,3-2,k,i).cwiseProduct(xz_b.block(3,3-2,k,i))));
    }
        else if(fd_order_num == 28){
        uu_b.block(4,4,k,i) = temp.block(4,4,k,i).cwiseProduct(uu_b.block(4,4,k,i)) - b.block(4,4,k,i).cwiseProduct(S81*(ca.block(4,4,k,i).cwiseProduct(xx_b.block(4,4,k,i)) - ca.block(4,4-1,k,i).cwiseProduct(xx_b.block(4,4-1,k,i))) +
                S82*(ca.block(4,4+1,k,i).cwiseProduct(xx_b.block(4,4+1,k,i)) - ca.block(4,4-2,k,i).cwiseProduct(xx_b.block(4,4-2,k,i))) +
            +S83*(ca.block(4,4+2,k,i).cwiseProduct(xx_b.block(4,4+2,k,i)) - ca.block(4,4-3,k,i).cwiseProduct(xx_b.block(4,4-3,k,i))) + S84*(ca.block(4,4+3,k,i).cwiseProduct(xx_b.block(4,4+3,k,i)) - ca.block(4,4-4,k,i).cwiseProduct(xx_b.block(4,4-4,k,i))) +
            S81*(cl.block(4,4,k,i).cwiseProduct(zz_b.block(4,4,k,i)) - cl.block(4,4-1,k,i).cwiseProduct(zz_b.block(4,4-1,k,i))) + S82*(cl.block(4,4+1,k,i).cwiseProduct(zz_b.block(4,4+1,k,i)) - cl.block(4,4-2,k,i).cwiseProduct(zz_b.block(4,4-2,k,i))) +
                S83*(cl.block(4,4+2,k,i).cwiseProduct(zz_b.block(4,4+2,k,i)) - cl.block(4,4-3,k,i).cwiseProduct(zz_b.block(4,4-3,k,i))) +
            S84*(cl.block(4,4+3,k,i).cwiseProduct(zz_b.block(4,4+3,k,i)) - cl.block(4,4-4,k,i).cwiseProduct(zz_b.block(4,4-4,k,i))) +
            S81*(cm1.block(4,4,k,i).cwiseProduct(xz_b.block(4,4,k,i)) - cm1.block(4-1,4,k,i).cwiseProduct(xz_b.block(4-1,4,k,i))) + S82*(cm1.block(4+1,4,k,i).cwiseProduct(xz_b.block(4+1,4,k,i)) - cm1.block(4-2,4,k,i).cwiseProduct(xz_b.block(4-2,4,k,i))) +
                S83*(cm1.block(4+2,4,k,i).cwiseProduct(xz_b.block(4+2,4,k,i)) - cm1.block(4-3,4,k,i).cwiseProduct(xz_b.block(4-3,4,k,i))) +
            S84*(cm1.block(4+3,4,k,i).cwiseProduct(xz_b.block(4+3,4,k,i)) - cm1.block(4-4,4,k,i).cwiseProduct(xz_b.block(4-4,4,k,i))));

        ww_b.block(4,4,k,i) = temp.block(4,4,k,i).cwiseProduct(ww_b.block(4,4,k,i)) - b1.block(4,4,k,i).cwiseProduct(S81*(cl.block(4+1,4,k,i).cwiseProduct(xx_b.block(4+1,4,k,i)) - cl.block(4,4,k,i).cwiseProduct(xx_b.block(4,4,k,i))) +
                S82*(cl.block(4+2,4,k,i).cwiseProduct(xx_b.block(4+2,4,k,i)) - cl.block(4-1,4,k,i).cwiseProduct(xx_b.block(4-1,4,k,i))) +
            +S83*(cl.block(4+3,4,k,i).cwiseProduct(xx_b.block(4+3,4,k,i)) - cl.block(4-2,4,k,i).cwiseProduct(xx_b.block(4-2,4,k,i))) + S84*(cl.block(4+4,4,k,i).cwiseProduct(xx_b.block(4+4,4,k,i)) - cl.block(4-3,4,k,i).cwiseProduct(xx_b.block(4-3,4,k,i))) +
            S81*(ca.block(4+1,4,k,i).cwiseProduct(zz_b.block(4+1,4,k,i)) - ca.block(4,4,k,i).cwiseProduct(zz_b.block(4,4,k,i))) + S82*(ca.block(4+2,4,k,i).cwiseProduct(zz_b.block(4+2,4,k,i)) - ca.block(4-1,4,k,i).cwiseProduct(zz_b.block(4-1,4,k,i))) +
                S83*(ca.block(4+3,4,k,i).cwiseProduct(zz_b.block(4+3,4,k,i)) - ca.block(4-2,4,k,i).cwiseProduct(zz_b.block(4-2,4,k,i))) +
            S84*(ca.block(4+4,4,k,i).cwiseProduct(zz_b.block(4+4,4,k,i)) - ca.block(4-3,4,k,i).cwiseProduct(zz_b.block(4-3,4,k,i))) +
            S81*(cm1.block(4,4+1,k,i).cwiseProduct(xz_b.block(4,4+1,k,i)) - cm1.block(4,4,k,i).cwiseProduct(xz_b.block(4,4,k,i))) + S82*(cm1.block(4,4+2,k,i).cwiseProduct(xz_b.block(4,4+2,k,i)) - cm1.block(4,4-1,k,i).cwiseProduct(xz_b.block(4,4-1,k,i))) +
                S83*(cm1.block(4,4+3,k,i).cwiseProduct(xz_b.block(4,4+3,k,i)) - cm1.block(4,4-2,k,i).cwiseProduct(xz_b.block(4,4-2,k,i))) +
            S84*(cm1.block(4,4+4,k,i).cwiseProduct(xz_b.block(4,4+4,k,i)) - cm1.block(4,4-3,k,i).cwiseProduct(xz_b.block(4,4-3,k,i))));}


    }
    void adjointStress(State& v) const {
        auto &uu_b=v.uu,&ww_b=v.ww,&xx_b=v.xx,&zz_b=v.zz,&xz_b=v.xz;
    if(fd_order_num == 22){
        xx_b.block(1,1,k,i)=temp.block(1,1,k,i).cwiseProduct(xx_b.block(1,1,k,i))-dtx*S21*(uu_b.block(1,1+1,k,i)-uu_b.block(1,1,k,i));
        zz_b.block(1,1,k,i)=temp.block(1,1,k,i).cwiseProduct(zz_b.block(1,1,k,i))-dtx*S21*(ww_b.block(1,1,k,i)-ww_b.block(1-1,1,k,i));
        xz_b.block(1,1,k,i)=temp.block(1,1,k,i).cwiseProduct(xz_b.block(1,1,k,i))-dtx*S21*(uu_b.block(1+1,1,k,i)-uu_b.block(1,1,k,i))-dtx*S21*(ww_b.block(1,1,k,i)-ww_b.block(1,1-1,k,i));}
    else if(fd_order_num == 24){
        xx_b.block(2,2,k,i)=temp.block(2,2,k,i).cwiseProduct(xx_b.block(2,2,k,i)) - dtx*(S41*(uu_b.block(2,2+1,k,i)-uu_b.block(2,2,k,i))+S42*(uu_b.block(2,2+2,k,i) - uu_b.block(2,2-1,k,i)));
        zz_b.block(2,2,k,i)=temp.block(2,2,k,i).cwiseProduct(zz_b.block(2,2,k,i)) - dtx*(S41*(ww_b.block(2,2,k,i)-ww_b.block(2-1,2,k,i))+S42*(ww_b.block(2+1,2,k,i) - ww_b.block(2-2,2,k,i)));
        xz_b.block(2,2,k,i)=temp.block(2,2,k,i).cwiseProduct(xz_b.block(2,2,k,i)) -
            dtx*(S41*(uu_b.block(2+1,2,k,i)-uu_b.block(2,2,k,i))+S42*(uu_b.block(2+2,2,k,i)-uu_b.block(2-1,2,k,i))) -
            dtx*(S41*(ww_b.block(2,2,k,i)-ww_b.block(2,2-1,k,i))+S42*(ww_b.block(2,2+1,k,i)-ww_b.block(2,2-2,k,i)));}
    else if(fd_order_num == 26){
        xx_b.block(3,3,k,i)=temp.block(3,3,k,i).cwiseProduct(xx_b.block(3,3,k,i))-dtx*(S61*(uu_b.block(3,3+1,k,i)-uu_b.block(3,3,k,i))+S62*(uu_b.block(3,3+2,k,i)-uu_b.block(3,3-1,k,i))+S63*(uu_b.block(3,3+3,k,i)-uu_b.block(3,3-2,k,i)));
        zz_b.block(3,3,k,i)=temp.block(3,3,k,i).cwiseProduct(zz_b.block(3,3,k,i))-dtx*(S61*(ww_b.block(3,3,k,i)-ww_b.block(3-1,3,k,i))+S62*(ww_b.block(3+1,3,k,i)-ww_b.block(3-2,3,k,i))+S63*(ww_b.block(3+2,3,k,i)-ww_b.block(3-3,3,k,i)));
        xz_b.block(3,3,k,i)=temp.block(3,3,k,i).cwiseProduct(xz_b.block(3,3,k,i)) -
            dtx*(S61*(uu_b.block(3+1,3,k,i)-uu_b.block(3,3,k,i))+S62*(uu_b.block(3+2,3,k,i)-uu_b.block(3-1,3,k,i))+S63*(uu_b.block(3+3,3,k,i)-uu_b.block(3-2,3,k,i))) -
            dtx*(S61*(ww_b.block(3,3,k,i)-ww_b.block(3,3-1,k,i))+S62*(ww_b.block(3,3+1,k,i)-ww_b.block(3,3-2,k,i))+S63*(ww_b.block(3,3+2,k,i)-ww_b.block(3,3-3,k,i)));}
    else if(fd_order_num == 28){
        xx_b.block(4,4,k,i)=temp.block(4,4,k,i).cwiseProduct(xx_b.block(4,4,k,i))-dtx*(S81*(uu_b.block(4,4+1,k,i)-uu_b.block(4,4,k,i))+S82*(uu_b.block(4,4+2,k,i)-uu_b.block(4,4-1,k,i)) +
            S83*(uu_b.block(4,4+3,k,i)-uu_b.block(4,4-2,k,i))+S84*(uu_b.block(4,4+4,k,i)-uu_b.block(4,4-3,k,i)));
        zz_b.block(4,4,k,i)=temp.block(4,4,k,i).cwiseProduct(zz_b.block(4,4,k,i))-dtx*(S81*(ww_b.block(4,4,k,i)-ww_b.block(4-1,4,k,i))+S82*(ww_b.block(4+1,4,k,i)-ww_b.block(4-2,4,k,i)) +
            S83*(ww_b.block(4+2,4,k,i)-ww_b.block(4-3,4,k,i))+S84*(ww_b.block(4+3,4,k,i)-ww_b.block(4-4,4,k,i)));
        xz_b.block(4,4,k,i)=temp.block(4,4,k,i).cwiseProduct(xz_b.block(4,4,k,i)) -
            dtx*(S81*(uu_b.block(4+1,4,k,i)-uu_b.block(4,4,k,i))+S82*(uu_b.block(4+2,4,k,i)-uu_b.block(4-1,4,k,i))+S83*(uu_b.block(4+3,4,k,i)-uu_b.block(4-2,4,k,i))+S84*(uu_b.block(4+4,4,k,i)-uu_b.block(4-3,4,k,i))) -
            dtx*(S81*(ww_b.block(4,4,k,i)-ww_b.block(4,4-1,k,i))+S82*(ww_b.block(4,4+1,k,i)-ww_b.block(4,4-2,k,i))+S83*(ww_b.block(4,4+2,k,i)-ww_b.block(4,4-3,k,i))+S84*(ww_b.block(4,4+3,k,i)-ww_b.block(4,4-4,k,i)));}

        zz_b.row(pad_top).setZero();
    }
};

// Save only the halo needed to reconstruct the undamped physical domain.
// In particular the free-surface projection and sponge are never inverted.
struct Boundary {
    std::vector<int> index;
    mutable std::vector<float> buffer;
    std::FILE* file=nullptr;
    mutable int cached=-2;
    explicit Boundary(const Elastic& e,bool save=true,const std::string& path="",bool read=false) {
        const int r=e.pad_top,z0=r+1,z1=z0+e.nz,x0=e.nbc,x1=x0+e.nx;
        for(int x=x0-r;x<x1+r;++x) for(int z=z0-r;z<z1+r;++z)
            if(z<z0 || z>=z1 || x<x0 || x>=x1) index.push_back(z+e.nzbc*x);
        if(!save) return;
        buffer.resize(5*index.size());
        file=path.empty()?std::tmpfile():std::fopen(path.c_str(),read?"rb":"w+b");
        if(!file) throw std::runtime_error("Cannot open the temporary reconstruction boundary file.");
        if(read) {
            seek(0,SEEK_END);
            if(tell()!=static_cast<long long>(buffer.size()*sizeof(float))*e.nt) {
                std::fclose(file);file=nullptr;
                throw std::runtime_error("Boundary file does not match nt or grid size.");
            }
        }
    }
    Boundary(const Boundary&)=delete;
    ~Boundary() {if(file) std::fclose(file);}
    void seek(long long offset,int origin) const {
#ifdef _WIN32
        int status=_fseeki64(file,offset,origin);
#else
        int status=fseeko(file,offset,origin);
#endif
        if(status) throw std::runtime_error("Cannot seek in the reconstruction boundary file.");
    }
    long long tell() const {
#ifdef _WIN32
        return _ftelli64(file);
#else
        return ftello(file);
#endif
    }
    void save(State& v,int it) {
        auto f=v.fields();
        for(int j=0;j<5;++j) for(size_t p=0;p<index.size();++p) buffer[j*index.size()+p]=f[j]->data()[index[p]];
        if(std::fwrite(buffer.data(),sizeof(float),buffer.size(),file)!=buffer.size())
            throw std::runtime_error("Cannot write the reconstruction boundary file.");
    }
    void restore(State& v,int it,bool stress) const {
        if(it>=0 && cached!=it) {
            seek(static_cast<long long>(buffer.size()*sizeof(float))*it,SEEK_SET);
            if(std::fread(buffer.data(),sizeof(float),buffer.size(),file)!=buffer.size())
                throw std::runtime_error("Cannot read the reconstruction boundary file.");
            cached=it;
        }
        auto f=v.fields();
        for(int j=stress?2:0;j<(stress?5:2);++j) for(size_t p=0;p<index.size();++p)
            f[j]->data()[index[p]]=it<0?0:buffer[j*index.size()+p];
    }
};
inline void reverse(const Elastic& e,State& v,const Strain& strain,const Boundary& boundary,
                    const Point& source,float amplitude,float source_before,int source_type,int it) {
    const int z=e.pad_top+1,x=e.nbc,nz=e.nz,nx=e.nx;
    // The interior damping is exactly one. Halo values come from the saved boundary.
    State increment(e.nzbc,e.nxbc);
    e.constitutive(strain,increment);
    v.xx.block(z,x,nz,nx)-=increment.xx.block(z,x,nz,nx)*e.dtx;
    v.zz.block(z,x,nz,nx)-=increment.zz.block(z,x,nz,nx)*e.dtx;
    v.xz.block(z,x,nz,nx)-=increment.xz.block(z,x,nz,nx)*e.dtx;
    if(source_type==5) source.prescribe(v.ww,source_before);
    else e.source(v,source,-amplitude,source_type);
    boundary.restore(v,it-1,true);
    increment.uu.setZero();increment.ww.setZero();
    increment.xx=v.xx;increment.zz=v.zz;increment.xz=v.xz;
    e.velocity(increment);
    v.uu.block(z,x,nz,nx)-=increment.uu.block(z,x,nz,nx);
    v.ww.block(z,x,nz,nx)-=increment.ww.block(z,x,nz,nx);
    boundary.restore(v,it-1,false);
}
} // namespace swd

#ifdef MATLAB_MEX_FILE
#include "mex.h"
namespace swd {
inline double number(const mxArray* a,size_t n) {
    if(mxIsSingle(a)) return static_cast<const float*>(mxGetData(a))[n];
    if(mxIsDouble(a)) return mxGetPr(a)[n];
    throw std::runtime_error("Expected a single or double array.");
}
inline MatrixXf matrix(const mxArray* a) {
    if(!a || mxIsComplex(a) || (!mxIsSingle(a)&&!mxIsDouble(a))) throw std::runtime_error("Expected a real numeric array.");
    MatrixXf m(mxGetM(a),mxGetN(a));
    for(size_t j=0;j<mxGetNumberOfElements(a);++j) m.data()[j]=float(number(a,j));
    return m;
}
inline mxArray* array(const MatrixXf& a) {
    mxArray* out=mxCreateNumericMatrix(a.rows(),a.cols(),mxSINGLE_CLASS,mxREAL);
    std::copy(a.data(),a.data()+a.size(),static_cast<float*>(mxGetData(out))); return out;
}
inline std::string filename(const mxArray* a) {
    if(!a||!mxIsChar(a)) throw std::runtime_error("Invalid temporary boundary filename.");
    char* p=mxArrayToString(a);if(!p) throw std::runtime_error("Cannot read the boundary filename.");
    std::string out(p);mxFree(p);return out;
}
inline Elastic model(const mxArray* const* in,bool forward) {
    Elastic e;
    e.nt=int(number(in[0],0));e.nzbc=int(number(in[0],1));e.nxbc=int(number(in[0],2));
    e.dtx=float(number(in[0],3));e.ng=int(number(in[0],4));
    e.fd_order_num=int(number(in[0],forward?11:9));
    e.nt_interval=int(number(in[0],forward?12:10));
    e.nz=int(number(in[0],forward?13:11));e.nx=int(number(in[0],forward?14:12));
    e.nbc=(e.nxbc-e.nx)/2;e.dt=forward?0:float(number(in[0],13));e.dimensions();
    e.temp=matrix(in[1]);e.ca=matrix(in[2]);e.cl=matrix(in[3]);e.cm=matrix(in[4]);
    e.cm1=matrix(in[5]);e.b=matrix(in[6]);e.b1=matrix(in[7]);
    if(e.nt<1||e.nt_interval<1||e.nbc<e.pad_top||e.nz<2||e.nx<2) throw std::runtime_error("Invalid grid or time dimensions.");
    for(const MatrixXf* p:{&e.temp,&e.ca,&e.cl,&e.cm,&e.cm1,&e.b,&e.b1})
        if(p->rows()!=e.nzbc||p->cols()!=e.nxbc) throw std::runtime_error("Material array dimensions do not match.");
    return e;
}
inline void geometry(Elastic& e,const mxArray* a) {
    if(!a || !mxIsStruct(a)) throw std::runtime_error("Invalid fixed grid.");
    e.mapped=mxGetScalar(mxGetField(a,0,"mapped"))!=0;
    if(!e.mapped) return;
    const mxArray* g=mxGetField(a,0,"G");
    if(!g||!mxIsSparse(g)||mxGetM(g)!=size_t(3*e.nzbc*e.nxbc)||mxGetN(g)!=size_t(2*e.nzbc*e.nxbc)) throw std::runtime_error("Grid dimensions do not match the model.");
    std::vector<Eigen::Triplet<float>> triplets;
    const mwIndex* row=mxGetIr(g);const mwIndex* col=mxGetJc(g);const double* val=mxGetPr(g);
    triplets.reserve(col[mxGetN(g)]);
    for(size_t x=0;x<mxGetN(g);++x) for(size_t p=col[x];p<col[x+1];++p) triplets.emplace_back(row[p],x,float(val[p]));
    e.G.resize(mxGetM(g),mxGetN(g));e.G.setFromTriplets(triplets.begin(),triplets.end());
    e.jn=matrix(mxGetField(a,0,"JN"));e.js=matrix(mxGetField(a,0,"JS"));
    e.invU=matrix(mxGetField(a,0,"invU"));e.invW=matrix(mxGetField(a,0,"invW"));
    float vmax=std::sqrt(e.ca.maxCoeff()),scale=float(mxGetScalar(mxGetField(a,0,"dt_scale")));
    if(e.dtx*vmax>0.5f*scale*1.001f) throw std::runtime_error("Time step is too large for the deformed grid. Multiply dt by grid.dt_scale before building the wavelet.");
}
inline std::vector<Point> points(const mxArray* a,const Elastic& e,bool vertical=true) {
    if(mxGetM(a)!=2 || int(mxGetN(a))!=e.ng) throw std::runtime_error("Expected 2-by-ng receiver coordinates.");
    std::vector<Point> out;
    for(int j=0;j<e.ng;++j) out.emplace_back(number(a,2*j)-1-(e.mapped&&vertical?.5:0),number(a,2*j+1)-1-(e.mapped&&!vertical?.5:0),e.nzbc,e.nxbc,e.mapped?e.pad_top+1:0);
    return out;
}
inline mxArray* stateArray(State& v) {
    mwSize dim[]={mwSize(v.uu.rows()),mwSize(v.uu.cols()),5};
    mxArray* out=mxCreateNumericArray(3,dim,mxSINGLE_CLASS,mxREAL);
    float* p=static_cast<float*>(mxGetData(out));
    for(auto f:v.fields()) {std::copy(f->data(),f->data()+f->size(),p);p+=f->size();}
    return out;
}
inline void readState(State& v,const mxArray* a) {
    if(!a||mxGetNumberOfElements(a)!=size_t(5*v.uu.size())) throw std::runtime_error("Invalid reconstruction state.");
    size_t j=0;for(auto f:v.fields()) for(int n=0;n<f->size();++n) f->data()[n]=float(number(a,j++));
}
}
#endif
