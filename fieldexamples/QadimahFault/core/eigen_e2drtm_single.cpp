/*
 *zhangchang23@mails.jlu.edu.cn
 *
 *
 *2021/11/6 Changchun
* e2drtm
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
class MexFunction : public matlab::mex::Function{
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
        int nt = inputs[0][0];
        int nzbc = inputs[0][1];
        int nxbc = inputs[0][2];
        float dtx = inputs[0][3];
        int ng = inputs[0][4];
        int nbc = inputs[0][5];
        int gz = inputs[0][6];gz--;
        int gx = inputs[0][7];gx--;
        int dg = inputs[0][8];
        int fd_order_num = inputs[0][9];
        
        int length_geophone = ng*dg;
        int nt_interval = inputs[0][10];
        int nz = inputs[0][11];
        int nx = inputs[0][12];
        float dt = inputs[0][13];
        int num_nt_record = nt/nt_interval;
        int wavefield_elements = num_nt_record*nx*nz;
        int number_elements = nz*nx;
        
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
       
        matlab::data::TypedArray<float> seismo_w_input = std::move(inputs[8]);
        auto seismo_w_ptr = getPointer(seismo_w_input);        
        MatrixXf seismo_w = Map<MatrixXf>(seismo_w_ptr,nt,ng);
        
        
        matlab::data::TypedArray<float> wavefield_gradient_fux_input = std::move(inputs[9]);
        auto wavefield_gradient_fux_ptr = getPointer(wavefield_gradient_fux_input);
        MatrixXf wavefield_gradient_fux = Map<MatrixXf>(wavefield_gradient_fux_ptr,nz,nx*num_nt_record);
        matlab::data::TypedArray<float> wavefield_gradient_fuz_input = std::move(inputs[10]);
        auto wavefield_gradient_fuz_ptr = getPointer(wavefield_gradient_fuz_input);
        MatrixXf wavefield_gradient_fuz = Map<MatrixXf>(wavefield_gradient_fuz_ptr,nz,nx*num_nt_record);
        matlab::data::TypedArray<float> wavefield_gradient_bwx_input = std::move(inputs[11]);
        auto wavefield_gradient_bwx_ptr = getPointer(wavefield_gradient_bwx_input);
        MatrixXf wavefield_gradient_bwx = Map<MatrixXf>(wavefield_gradient_bwx_ptr,nz,nx*num_nt_record);
        matlab::data::TypedArray<float> wavefield_gradient_bwz_input = std::move(inputs[12]);
        auto wavefield_gradient_bwz_ptr = getPointer(wavefield_gradient_bwz_input);
        MatrixXf wavefield_gradient_bwz = Map<MatrixXf>(wavefield_gradient_bwz_ptr,nz,nx*num_nt_record);
        
        
//      Eigen Initialising input variables: uu, ww, xx, xz, zz
        MatrixXf uu_b(nzbc,nxbc);  
        uu_b << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf ww_b(nzbc,nxbc);  
        ww_b << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf xx_b(nzbc,nxbc);  
        xx_b << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf xz_b(nzbc,nxbc);  
        xz_b << MatrixXf::Zero(nzbc,nxbc);
        MatrixXf zz_b(nzbc,nxbc);  
        zz_b << MatrixXf::Zero(nzbc,nxbc);
//      Eigen Initialising input variables: cl_img, cm_img
        MatrixXf cl_img(nz,nx);  
        cl_img << MatrixXf::Zero(nz,nx);
        MatrixXf cm_img(nz,nx);  
        cm_img << MatrixXf::Zero(nz,nx);
        MatrixXf illum_div(nz,nx);
        illum_div << MatrixXf::Zero(nz,nx);
       
//      Eigen zero_vector for free surface zz
        VectorXf zero_vector(nxbc);
        zero_vector << VectorXf::Zero(nxbc);
        VectorXf geophone_vector(nxbc);
        geophone_vector << VectorXf::Zero(nxbc);
        VectorXf residual_vector(nx);
        
        int k;int i;int pad_top; 
        if(fd_order_num==22){
            k = nzbc-2; i = nxbc-2; pad_top = 1;}
        else if(fd_order_num==24){
            k = nzbc-4; i = nxbc-4; pad_top = 2;}
        else if(fd_order_num==26){
            k = nzbc-6; i = nxbc-6; pad_top = 3;}
        else if(fd_order_num==28){
            k = nzbc-8; i = nxbc-8; pad_top = 4;}
        float S21 = 1;
        float S41 = 1.1250; float S42 = -0.0416666667;
        float S61 = 1.17187; float S62 = -6.51042E-2; float S63 = 4.68750E-3;
        float S81 = 1.19629; float S82 = -7.97526E-2; float S83 = 9.57031E-3; float S84 = -6.97545E-4;    
        
//         MatrixXf wavefield_gradient_wwb(nzbc,nxbc*nt);
//         wavefield_gradient_wwb << MatrixXf::Zero(nzbc,nxbc*nt);  
   for( int it = nt-2; it > -1; it-- ){
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
        
        
       geophone_vector = ww_b.row(gz);
       
       residual_vector = seismo_w.row(it);
       geophone_vector(seq(gx,gx+length_geophone-1,dg)) += dt*residual_vector;//should be dt*residual_vector/den; but den matrix is ones(nz,nx) as default. 
       ww_b.row(gz) = geophone_vector;
       
//      update stress
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
     
        zz_b.row(pad_top) = zero_vector;
        
//                    wavefield_gradient_wwb.block(0,nxbc*it,nzbc,nxbc)=ww_b;//
                   
       if(it%nt_interval == 0){
            cl_img+=(wavefield_gradient_fux.block(0,nx*it/nt_interval,nz,nx) + wavefield_gradient_bwz.block(0,nx*it/nt_interval,nz,nx)).cwiseProduct((xx_b.block(pad_top+1,nbc,nz,nx)+zz_b.block(pad_top+1,nbc,nz,nx)));
            cm_img+=2*(xx_b.block(pad_top+1,nbc,nz,nx).cwiseProduct(wavefield_gradient_fux.block(0,nx*it/nt_interval,nz,nx)) + 
                    zz_b.block(pad_top+1,nbc,nz,nx).cwiseProduct(wavefield_gradient_bwz.block(0,nx*it/nt_interval,nz,nx)) + 
                    0.5*xz_b.block(pad_top+1,nbc,nz,nx).cwiseProduct((wavefield_gradient_fuz.block(0,nx*it/nt_interval,nz,nx) + wavefield_gradient_bwx.block(0,nx*it/nt_interval,nz,nx))));  
            illum_div = illum_div+(wavefield_gradient_fux.block(0,nx*it/nt_interval,nz,nx) + wavefield_gradient_bwz.block(0,nx*it/nt_interval,nz,nx)).cwiseProduct((wavefield_gradient_fux.block(0,nx*it/nt_interval,nz,nx) + wavefield_gradient_bwz.block(0,nx*it/nt_interval,nz,nx)));
            

       }
}
//        std::vector<size_t> size_output_wavefield(1,nzbc);//
//        size_output_wavefield.insert(size_output_wavefield.end(),nxbc);//
//        size_output_wavefield.insert(size_output_wavefield.end(),1700);//
       
       std::vector<size_t> size_output(1,nz);
       size_output.insert(size_output.end(),nx);
       float* cl_img_ptr1 = cl_img.data();
       float* cm_img_ptr2 = cm_img.data();
       float* illum_div_ptr3 = illum_div.data();
//        float* wavefield_gradient_wwb_ptr = wavefield_gradient_wwb.data();//
       outputs[0] = factory.createArray<float>(size_output,cl_img_ptr1,cl_img_ptr1 + number_elements);
       outputs[1] = factory.createArray<float>(size_output,cm_img_ptr2,cm_img_ptr2 + number_elements);
       outputs[2] = factory.createArray<float>(size_output,illum_div_ptr3,illum_div_ptr3 + number_elements);
//        outputs[3] = factory.createArray<float>(size_output_wavefield,wavefield_gradient_wwb_ptr,wavefield_gradient_wwb_ptr + 23120000);//
    }};
