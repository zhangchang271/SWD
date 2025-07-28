/*
 *zhangchang23@mails.jlu.edu.cn
 *
 *
 *2021/11/4  Changchun
* staggerfd
*/
/* License:
 *
 * GNU GENERAL PUBLIC LICENSE
 * Version 3, 29 June 2007
 *
 * Copyright (C) 2024 Zhang Chang
 *
 * This file is part of SWD.
 *
 * SWD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#define EIGEN_USE_MKL_ALL
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include <iostream>
#include "type_traits"


// using namespace matlab::data;
using namespace std;
using namespace Eigen;
using matlab::mex::ArgumentList;
class MexFunction : public matlab::mex::Function {
// // // // // // // // // // // // // // // // // // // // // //     
    	//! Extracts the pointer to underlying data from the non-const iterator (`TypedIterator<T>`).
	/*! This function does not throw any exceptions. */
	template <typename T>
	inline T* toPointer(const matlab::data::TypedIterator<T>& it) MW_NOEXCEPT {
		static_assert(std::is_arithmetic<T>::value && !std::is_const<T>::value,
			"Template argument T must be a std::is_arithmetic and non-const type.");
		return it.operator->();
	}
	template <typename T>
	inline T* getPointer(matlab::data::TypedArray<T>& arr) MW_NOEXCEPT {
		static_assert(std::is_arithmetic<T>::value, "Template argument T must be a std::is_arithmetic type.");
		return toPointer(arr.begin());
	}
	template <typename T>
	inline const T* getPointer(const matlab::data::TypedArray<T>& arr) MW_NOEXCEPT {
		return getPointer(const_cast<matlab::data::TypedArray<T>&>(arr));
	}
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //     
public:

    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        matlab::data::ArrayFactory factory;
/*  input:  nt
 *          nzbc
 *          nxbc
 *          dtx
 *          temp
 *          ng
 *          ca
 *          cl
 *          cm
 *          b
 *          s
 *          sourcetype
 *Eigen init:
 *          uu
 *          ww
 *          xx
 *          xz
 *          zz
 *          fux
 *          fuz
 *          bwx
 *          bwz
 *return:
 *          seismo_w
 *          seismo_u
 *          wavefield_gradient
 */
//      Input constant from matlab workspace: nt, nzbc, nxbc, dtx, ng, source_type_num, fd_order_num
        int nt = inputs[0][0];
        int nzbc = inputs[0][1];
        int nxbc = inputs[0][2];
        float dtx = inputs[0][3];
        int ng = inputs[0][4];
        int sz = inputs[0][5];sz--;
        int sx = inputs[0][6];sx--;
        int gz = inputs[0][7];gz--;
        int gx = inputs[0][8];gx--;
        int dg = inputs[0][9];
        int source_type_num = inputs[0][10];
        int fd_order_num = inputs[0][11];
        int number_elements = nt*ng;
        int length_geophone = ng*dg;
        int nt_interval = inputs[0][12];
        int nz = inputs[0][13];
        int nx = inputs[0][14];
        int format_num = inputs[0][15];
        int nbc = (nxbc-nx)/2;
        int num_nt_record = nt/nt_interval;
        int wavefield_elements = num_nt_record*nx*nz;
//      Input variables from matlab workspace: temp ca cl cm b s
        matlab::data::TypedArray<float> temp_input = std::move(inputs[1]);
        auto temp_ptr = getPointer(temp_input);
        std::vector<size_t> size_input;
        size_input = temp_input.getDimensions();
        MatrixXf temp = Map<MatrixXf>(temp_ptr,nzbc,nxbc);
        
        matlab::data::TypedArray<float> ca_input = std::move(inputs[2]);
        auto ca_ptr = getPointer(ca_input);
        MatrixXf ca = Map<MatrixXf>(ca_ptr,nzbc,nxbc);
        
        matlab::data::TypedArray<float> cl_input = std::move(inputs[3]);
        auto cl_ptr = getPointer(cl_input);        
        MatrixXf cl = Map<MatrixXf>(cl_ptr,nzbc,nxbc);
        
        matlab::data::TypedArray<float> cm_input = std::move(inputs[4]);
        auto cm_ptr = getPointer(cm_input);
        MatrixXf cm = Map<MatrixXf>(cm_ptr,nzbc,nxbc);
        
        matlab::data::TypedArray<float> cm1_input = std::move(inputs[5]);
        auto cm1_ptr = getPointer(cm1_input);
        MatrixXf cm1 = Map<MatrixXf>(cm1_ptr,nzbc,nxbc);
        
        matlab::data::TypedArray<float> b_input = std::move(inputs[6]);
        auto b_ptr = getPointer(b_input);        
        MatrixXf b = Map<MatrixXf>(b_ptr,nzbc,nxbc);
        
        matlab::data::TypedArray<float> b1_input = std::move(inputs[7]);
        auto b1_ptr = getPointer(b1_input);        
        MatrixXf b1 = Map<MatrixXf>(b1_ptr,nzbc,nxbc);
       
        matlab::data::TypedArray<float> s_input = std::move(inputs[8]);
        auto s_ptr = getPointer(s_input);        
        VectorXf s = Map<VectorXf>(s_ptr,nt,1);
//      Eigen Initialising input variables: uu, ww, xx, xz, zz
        MatrixXf uu(nzbc,nxbc);  
        uu << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf ww(nzbc,nxbc);  
        ww << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf xx(nzbc,nxbc);  
        xx << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf xz(nzbc,nxbc);  
        xz << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf zz(nzbc,nxbc);  
        zz << MatrixXf::Zero(nzbc,nxbc);
//      Eigen Initialising input variables: fux, fuz, bwx, bwz
        MatrixXf fux(nzbc,nxbc);  
        fux << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf fuz(nzbc,nxbc);  
        fuz << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf bwx(nzbc,nxbc);  
        bwx << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf bwz(nzbc,nxbc);  
        bwz << MatrixXf::Zero(nzbc,nxbc);
//      Eigen Initialising output variables: seismo_w, seismo_u       
        MatrixXf seismo_w(nt,ng);  
        seismo_w << MatrixXf::Zero(nt,ng);
        MatrixXf seismo_u(nt,ng);  
        seismo_u << MatrixXf::Zero(nt,ng);
        
        MatrixXf wavefield_gradient_fux(nz,nx*num_nt_record);
        wavefield_gradient_fux << MatrixXf::Zero(nz,nx*num_nt_record);
        MatrixXf wavefield_gradient_fuz(nz,nx*num_nt_record);
        wavefield_gradient_fuz << MatrixXf::Zero(nz,nx*num_nt_record);
        MatrixXf wavefield_gradient_bwx(nz,nx*num_nt_record);
        wavefield_gradient_bwx << MatrixXf::Zero(nz,nx*num_nt_record);
        MatrixXf wavefield_gradient_bwz(nz,nx*num_nt_record);
        wavefield_gradient_bwz << MatrixXf::Zero(nz,nx*num_nt_record);
       
//      Eigen zero_vector for free surface zz
        VectorXf zero_vector(nxbc);
        zero_vector << VectorXf::Zero(nxbc);
        VectorXf geophone_vector(nxbc);
        geophone_vector << VectorXf::Zero(nxbc);
        int k;int i;int pad_top;
        if(fd_order_num==22){
            k = nzbc-2; i = nxbc-2; pad_top = 1;}
        else if(fd_order_num==24){
            k = nzbc-4; i = nxbc-4; pad_top = 2;}
        else if(fd_order_num==26){
            k = nzbc-6; i = nxbc-6; pad_top = 3;}
        else if(fd_order_num==28){
            k = nzbc-8; i = nxbc-8; pad_top = 4;}

        float S41 = 1.1250; float S42 = -0.0416666667;
        float S61 = 1.17187; float S62 = -6.51042E-2; float S63 = 4.68750E-3;
        float S81 = 1.19629; float S82 = -7.97526E-2; float S83 = 9.57031E-3; float S84 = -6.97545E-4;
// forward differential iterations
   for( int it = 0; it < nt; it++ ){
//        if(source_type_num == 1){
//            xx(sz,sx)=xx(sz,sx)+s(it);
//            zz(sz,sx)=zz(sz,sx)+s(it);}
//        else if(source_type_num == 2){
//            uu(sz,sx)=uu(sz,sx)+s(it);
//            uu(sz-1,sx)=uu(sz-1,sx)-s(it);
//            ww(sz,sx)=ww(sz,sx)+s(it);
//            ww(sz,sx+1)=ww(sz,sx+1)+s(it);}
//        else if(source_type_num == 3){
//            uu(sz,sx)=uu(sz,sx)+s(it);
//            ww(sz,sx)=ww(sz,sx)+s(it);}
//        else if(source_type_num == 4){
//            zz(sz,sx)=zz(sz,sx)+s(it);}
//        else if(source_type_num == 5){
//            ww(sz,sx)=s(it);}
       
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
// if (fd_order==22)
//         uu(k,i)=temp(k,i).*uu(k,i)+b(k,i).*(S21*(xx(k,i+1)-xx(k,i)+...
//             xz(k,i)-xz(k-1,i)));
//         ww(k,i)=temp(k,i).*ww(k,i)+b1(k,i).*(S21*(xz(k,i)-xz(k,i-1)+...
//             zz(k+1,i)-zz(k,i)));
//     elseif (fd_order==24)
//         uu(k,i)=temp(k,i).*uu(k,i)+b(k,i).*(S41*(xx(k,i)-xx(k,i-1))+S42*(xx(k,i+1)-xx(k,i-2))+...
//             S41*(xz(k,i)-xz(k-1,i))+S42*(xz(k+1,i)-xz(k-2,i)));
//         
//         ww(k,i)=temp(k,i).*ww(k,i)+b(k,i).*(S41*(xz(k,i+1)-xz(k,i))+S42*(xz(k,i+2)-xz(k,i-1))+...
//             S41*(zz(k+1,i)-zz(k,i))+S42*(zz(k+2,i)-zz(k-1,i)));
//     elseif (fd_order==26)
//         uu(k,i)=temp(k,i).*uu(k,i)+b(k,i).*(S61*(xx(k,i)-xx(k,i-1))+S62*(xx(k,i+1)-xx(k,i-2))+...
//             +S63*(xx(k,i+2)-xx(k,i-3))+S61*(xz(k,i)-xz(k-1,i))+S62*...
//             (xz(k+1,i)-xz(k-2,i))+S63*(xz(k+2,i)-xz(k-3,i)));
//         
//         ww(k,i)=temp(k,i).*ww(k,i)+b(k,i).*(S61*(xz(k,i+1)-xz(k,i))+S62*(xz(k,i+2)-xz(k,i-1))+...
//             +S63*(xz(k,i+3)-xz(k,i-2))+S61*(zz(k+1,i)-zz(k,i))+S62*...
//             (zz(k+2,i)-zz(k-1,i))+S63*(zz(k+3,i)-zz(k-2,i)));
//     elseif (fd_order==28)
//         uu(k,i)=temp(k,i).*uu(k,i)+b(k,i).*(S81*(xx(k,i)-xx(k,i-1))+S82*(xx(k,i+1)-xx(k,i-2))+...
//             +S83*(xx(k,i+2)-xx(k,i-3))+S84*(xx(k,i+3)-xx(k,i-4))+...
//             S81*(xz(k,i)-xz(k-1,i))+S82*(xz(k+1,i)-xz(k-2,i))+S83*(xz(k+2,i)-xz(k-3,i))+...
//             S84*(xz(k+3,i)-xz(k-4,i)));
//         
//         ww(k,i)=temp(k,i).*ww(k,i)+b(k,i).*(S81*(xz(k,i+1)-xz(k,i))+S82*(xz(k,i+2)-xz(k,i-1))+...
//             S83*(xz(k,i+3)-xz(k,i-2))+S84*(xz(k,i+4)-xz(k,i-3))+...
//             S81*(zz(k+1,i)-zz(k,i))+S82*(zz(k+2,i)-zz(k-1,i))+S83*(zz(k+3,i)-zz(k-2,i))+...
//             S84*(zz(k+4,i)-zz(k-3,i)));
//     end
       if(source_type_num == 1){
           xx(sz,sx)=xx(sz,sx)+s(it);
           zz(sz,sx)=zz(sz,sx)+s(it);}
       else if(source_type_num == 2){
           uu(sz,sx)=uu(sz,sx)+s(it);
           uu(sz-1,sx)=uu(sz-1,sx)-s(it);
           ww(sz,sx)=ww(sz,sx)+s(it);
           ww(sz,sx+1)=ww(sz,sx+1)+s(it);}
       else if(source_type_num == 3){
           uu(sz,sx)=uu(sz,sx)+s(it);
           ww(sz,sx)=ww(sz,sx)+s(it);}
       else if(source_type_num == 4){
           zz(sz,sx)=zz(sz,sx)+s(it);}
       else if(source_type_num == 5){
           ww(sz,sx)=s(it);}
       
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
//         fux(k,i)=S21*(uu(k,i)-uu(k,i-1));
//         fuz(k,i)=S21*(uu(k+1,i)-uu(k,i));
//         bwx(k,i)=S21*(ww(k,i+1)-ww(k,i));
//         bwz(k,i)=S21*(ww(k,i)-ww(k-1,i)); 
       xx=temp.cwiseProduct(xx) + (ca.cwiseProduct(fux) + cl.cwiseProduct(bwz))*dtx;
       zz=temp.cwiseProduct(zz) + (ca.cwiseProduct(bwz) + cl.cwiseProduct(fux))*dtx;
       xz=temp.cwiseProduct(xz) + (cm1.cwiseProduct(fuz + bwx))*dtx;
//     xx=temp.*xx+(ca.*fux+cl.*bwz)*dtx;
//     zz=temp.*zz+(ca.*bwz+cl.*fux)*dtx;
//     xz=temp.*xz+cm.*(fuz+bwx)*dtx;
       zz.row(pad_top) = zero_vector;//临时修改5.20
//        zz(pad_top,:)=0.0;
       geophone_vector=ww.row(gz);
       seismo_w.row(it) = geophone_vector(seq(gx,gx+length_geophone,dg));
       geophone_vector=uu.row(gz);
       seismo_u.row(it) = geophone_vector(seq(gx,gx+length_geophone,dg));
       if(it%nt_interval==0){
       wavefield_gradient_fux.block(0,nx*it/nt_interval,nz,nx)=fux.block(pad_top+1,nbc,nz,nx);
       wavefield_gradient_fuz.block(0,nx*it/nt_interval,nz,nx)=fuz.block(pad_top+1,nbc,nz,nx);
       wavefield_gradient_bwx.block(0,nx*it/nt_interval,nz,nx)=bwx.block(pad_top+1,nbc,nz,nx);
       wavefield_gradient_bwz.block(0,nx*it/nt_interval,nz,nx)=bwz.block(pad_top+1,nbc,nz,nx);
       }
//        for ig=1:ng
//           seismo_w(it,ig)=ww(gz(ig),gx(ig));
//           seismo_u(it,ig)=uu(gz(ig),gx(ig));
//        end
       
   }
          
       float* eig_seismo_u_ptr = seismo_u.data();
       float* eig_seismo_w_ptr2 = seismo_w.data();
//        float* eig_ww_ptr2 = ww.data();
       std::vector<size_t> size_output(1,nt);
//        size_output.insert(size_output.begin(),nt);
       size_output.insert(size_output.end(),ng);
       std::vector<size_t> size_output_wavefield(1,nz);
       if(format_num==2){
           size_output_wavefield.insert(size_output_wavefield.end(),nx*num_nt_record);}   //export output wavefield gradient as 2D matrix; easy to calculate adjoint source with eigen;
       else if(format_num==3){
       size_output_wavefield.insert(size_output_wavefield.end(),nx);
       size_output_wavefield.insert(size_output_wavefield.end(),num_nt_record);}//export output wavefield gradient as 3D matrix; It has same format as original matlab code; 
       float* wavefield_gradient_fux_ptr1 = wavefield_gradient_fux.data();
       float* wavefield_gradient_fuz_ptr2 = wavefield_gradient_fuz.data();
       float* wavefield_gradient_bwx_ptr3 = wavefield_gradient_bwx.data();
       float* wavefield_gradient_bwz_ptr4 = wavefield_gradient_bwz.data();
       matlab::data::StructArray S = factory.createStructArray({ 1,1 }, { "fux", "fuz","bwx","bwz" });
       S[0]["fux"] = factory.createArray<float>(size_output_wavefield,wavefield_gradient_fux_ptr1,wavefield_gradient_fux_ptr1+wavefield_elements);//It is column-major so that the shape of fux can be automaticly reshaped in 3D matrix.
       S[0]["fuz"] = factory.createArray<float>(size_output_wavefield,wavefield_gradient_fuz_ptr2,wavefield_gradient_fuz_ptr2+wavefield_elements);
       S[0]["bwx"] = factory.createArray<float>(size_output_wavefield,wavefield_gradient_bwx_ptr3,wavefield_gradient_bwx_ptr3+wavefield_elements);
       S[0]["bwz"] = factory.createArray<float>(size_output_wavefield,wavefield_gradient_bwz_ptr4,wavefield_gradient_bwz_ptr4+wavefield_elements);
       
       outputs[0] = factory.createArray<float>(size_output,eig_seismo_u_ptr,eig_seismo_u_ptr + number_elements);
       outputs[1] = factory.createArray<float>(size_output,eig_seismo_w_ptr2,eig_seismo_w_ptr2 + number_elements);
       outputs[2] = S;
    } 

    void arrayProduct(matlab::data::TypedArray<float>& matrix, float multiplier) {
//     cout <<"matrix:" <<  matrix << endl;
//     cout <<"&matrix:" << & matrix << endl;
    }


};