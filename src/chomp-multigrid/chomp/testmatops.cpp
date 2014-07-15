/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#include "Chomp.h"
#include "../mzcommon/gauss.h"

#include <assert.h>
#include <mzcommon/TimeUtil.h>

using namespace chomp;

void regularChol(const MatX& A, MatX& L) {

  int n = A.rows();
  L = MatX(n,n);
  L.setZero();
  
  for (int j=0; j<n; ++j) {
    for (int i=j; i<n; ++i) {
      double sum = 0;
      for (int k=0; k<j; ++k) { 
        sum += L(i,k) * L(j,k); 
      }
      if (i == j) {
        L(j,j) = sqrt(A(j,j) - sum);
      } else {
        L(i,j) = ( A(i,j) - sum ) / L(j,j);
      }
    }
  }
 

}

// Ax = b
//   L L^T x = b
// multiply both sides by L^-1
//   L^T x = L^-1 b
// then multiply both sides by L^-T
//  x = L^-T L^-1 b

void regularCholSolve(const MatX& L, MatX& x) {

  assert( L.rows() == L.cols() );

  int n = L.rows();
  
  assert(x.rows() == n);

#if 0

  for (int c=0; c<x.cols(); ++c) {

    for (int i=0; i<n; ++i) {
      for (int j=0; j<i; ++j) {
        x(i,c) -= L(i,j)*x(j,c); // here j < i so col < row
      }
      x(i,c) /= L(i,i);
    }

    for (int i=n-1; i>=0; --i) {
      for (int j=i+1; j<n; ++j) {
        x(i,c) -= L(j,i) * x(j,c); // here j > i so col < row
      }
      x(i,c) /= L(i,i);
    }

  }

#else

  for (int i=0; i<n; ++i) {
    for (int j=0; j<i; ++j) {
      x.row(i) -= L(i,j)*x.row(j); // here j < i so col < row
    }
    x.row(i) /= L(i,i);
  }

  for (int i=n-1; i>=0; --i) {
    for (int j=i+1; j<n; ++j) {
      x.row(i) -= L(j,i) * x.row(j); // here j > i so col < row
    }
    x.row(i) /= L(i,i);
  }

#endif

}

double relErr(const MatX& m1, const MatX& m2) {

  return ( (m2-m1).lpNorm<Eigen::Infinity>() /
           std::max(m1.lpNorm<Eigen::Infinity>(), 
                    m2.lpNorm<Eigen::Infinity>()) );

}

void testMaps() {

  MatX A(4,4);

  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      A(i,j) = 4*i + j;
    }
  }

  MatX At = A.transpose();


  typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> StrideX;

  Eigen::Map< MatX, 0, StrideX > At2(A.data(), 4, 4, StrideX(1,4));

  std::cout << "A   = \n" << A << "\n\n";
  std::cout << "At  = \n" << A.transpose() << "\n\n";
  std::cout << "At2 = \n" << At2 << "\n";

  std::cout << "A.data()   = " << A.data() << "\n";
  std::cout << "At.data()  = " << At.data() << "\n";
  std::cout << "At2.data() = " << At2.data() << "\n";

  std::cout << "A*At = \n" << A*At << "\n";
  std::cout << "A*At2 = \n" << A*At2 << "\n";
  
}

int primary( int argc, char ** argv ){
  int n = 7;

  if (argc > 1) { 
    int nn = atoi(argv[1]);
    if (nn) { n = nn; }
  }

  MatX A(n,n);
  A.setZero();


  for (int i=0; i<n; ++i) {
    if (i+2 < n) { A(i,i+2) = 1; }
    if (i+1 < n) { A(i,i+1) = -4; }
    A(i,i) = 6;
    if (i > 0) { A(i,i-1) = -4; }
    if (i > 1) { A(i,i-2) = 1; }
  }

  MatX coeffs(1,3);
  
  coeffs << 1, -4, 6;

  MatX Ainv = A.inverse();

  MatX L, Ls;
  MatX testvec = MatX::Random(n, 2 );
  MatX testvec2 = testvec;
  MatX m1(n,n); m1.setIdentity();
  MatX m2(n,n); m2.setIdentity();
  MatX m3(n,n); m3.setIdentity();
  MatX m4(n,n); m4.setIdentity();
  MatX work;

  Eigen::LLT<MatX> cholSolver(A);

  TimeStamp t0 = TimeStamp::now();

  //regularChol(A, L);
  cholSolver.compute(A);

  TimeStamp t1 = TimeStamp::now();

  //regularCholSolve(L, m1);
  m1 = cholSolver.solve(m1);

  TimeStamp t2 = TimeStamp::now();

  skylineChol(n, coeffs, Ls);

  TimeStamp t3 = TimeStamp::now();
  
  skylineCholSolve(Ls, m2);

  TimeStamp t4 = TimeStamp::now();

  diagMul(coeffs, m4, m3);

  TimeStamp t5 = TimeStamp::now();
    
  MatX x0(1,1), x1(1,1);
  x0(0) = -1;
  x1(0) = 1;

  //MatX b(n,1);
  //createBMatrix(n, coeffs, x0, x1, b);

  double d1 = (t1-t0).toDouble();
  double d2 = (t2-t1).toDouble();
  double d3 = (t3-t2).toDouble();
  double d4 = (t4-t3).toDouble();
  double d5 = (t5-t4).toDouble();

  
  if (n < 20) { 
    std::cout << "A =\n" << A << "\n";
    std::cout << "A =\n" << m3 << "\n";
    std::cout << "A.inv() =\n" << A.inverse() << "\n";
    std::cout << "m1 =\n" << m1 << "\n";
    std::cout << "errm1 = \n" << m1-A.inverse() << "\n";
    std::cout << "m2 =\n" << m2 << "\n";
    std::cout << "errm2 = \n" << m2-A.inverse() << "\n";
  }
  
  assert( relErr(A, m3) < 1e-5 );
  assert( relErr(m1, m2) < 1e-5 );

  std::cout << "with n=" << n << ", regular cholesky decomp took " << d1 << "s.\n";
  std::cout << "with n=" << n << ", regular cholesky solve took " << d2 << "s.\n";
  std::cout << "with n=" << n << ", band cholesky decomp took " << d3 << "s. (" << 100*d3/d1 << "% of regular)\n";
  std::cout << "with n=" << n << ", band cholesky solve took " << d4 << "s. (" << 100*d4/d2 << "% of regular)\n";
  std::cout << "with n=" << n << ", band multiply took " << d5 << ".s\n";

  return 0;


}

int secondary( int argc, char ** argv ){
  int n = 7;

  if (argc > 1) { 
    int nn = atoi(argv[1]);
    if (nn) { n = nn; }
  }

  MatX A(n,n);
  A.setZero();


  for (int i=0; i<n; ++i) {
    if (i+2 < n) { A(i,i+2) = 1; }
    if (i+1 < n) { A(i,i+1) = -4; }
    A(i,i) = 6;
    if (i > 0) { A(i,i-1) = -4; }
    if (i > 1) { A(i,i-2) = 1; }
  }

  MatX coeffs(1,3);
  
  coeffs << 1, -4, 6;

  MatX Ainv = A.inverse();

  MatX Ls;
  MatX testvec = MatX::Random(n, 2 );
  for ( int i = 0; i < testvec.size() ; i ++ ){
      testvec(i) = gauss_ziggurat_standard();
  }

  std::cout << "Testvec= \n" << testvec << "\n\n";

  MatX testvec2 = testvec;
  MatX m1(n,n); m1.setIdentity();
  MatX m2(n,n); m2.setIdentity();
  MatX m3(n,n); m3.setIdentity();
  MatX m4(n,n); m4.setIdentity();
  MatX m5(n,n); m5.setIdentity();
  MatX work;

  Eigen::LLT<MatX> cholSolver(A);
  Eigen::LLT<MatX> cholSolver2(Ainv);

  //regularCholSolve
  m1 = cholSolver.solve(m1);

  //our cholesky solver
  skylineChol(n, coeffs, Ls);
  skylineCholSolve(Ls, m2);

    
  diagMul(coeffs, m4, m3);


  MatX ell = cholSolver2.matrixL();
  MatX ell_mult = ell * testvec;
  

  MatX L = cholSolver.matrixL();
  MatX L_inv_trans = L.inverse().transpose();
  MatX L_mult = L_inv_trans * testvec;

  
  skylineCholMultiplyInverse( Ls, m5 );
  skylineCholMultiplyInverseTranspose( Ls, testvec);
  

  L = cholSolver.matrixL();
  if (n < 20) { 
    std::cout << "A =\n" << A << "\n";
    std::cout << "A =\n" << m3 << "\n";
    std::cout << "A.inv() =\n" << A.inverse() << "\n";
    std::cout << "m1 =\n" << m1 << "\n";
    std::cout << "errm1 = \n" << m1-A.inverse() << "\n";
    std::cout << "m2 =\n" << m2 << "\n";
    std::cout << "errm2 = \n" << m2-A.inverse() << "\n";
    std::cout << "L-1 * identity = \n" << m5 <<"\n\n";
    std::cout << "testvec = \n" << testvec <<"\n\n";
    std::cout << "ell_mult = \n" << ell_mult <<"\n\n";
    std::cout << "L_mult = \n" << L_mult <<"\n\n";

    std::cout << "L = \n" << L  << "\n\n";
    std::cout << "ELL = \n" << ell <<"\n\n";
    
    std::cout << "L^-1 = \n" << L.inverse()  << "\n\n";
    std::cout << "ELL^-1 = \n" << ell.inverse() <<"\n\n";

    std::cout << "L^-T = \n" << L_inv_trans  << "\n\n";
    std::cout << "ELL^-T = \n" << ell.inverse().transpose()  << "\n\n";
    
    std::cout << "ELL*ELL^T = \n" << ell*ell.transpose()  << "\n\n";
    std::cout << "(ELL*ELL^T)^-1 = \n" << (ell*ell.transpose()).inverse()  << "\n\n";
    std::cout << "L*L^T = \n" << L*L.transpose()  << "\n\n";
    std::cout << "Diff = \n" << L - ell.inverse().transpose() << "\n\n";

  }
  
  assert( relErr(A, m3) < 1e-5 );
  assert( relErr(m1, m2) < 1e-5 );

  assert( relErr(m1, m2) < 1e-5 );
  assert( relErr(testvec, ell_mult) < 1e-5 );

  return 0;

}


//this tests goal set operations
int tertiary( int argc, char ** argv ){
  int n = 7;

  if (argc > 1) { 
    int nn = atoi(argv[1]);
    if (nn) { n = nn; }
  }

  MatX A(n,n);
  A.setZero();
  
  MatX coeffs, gs_coeffs;

  bool acceleration = false;

  if ( acceleration ){
    coeffs.resize(1,3);
    gs_coeffs.resize(2,2);
  
    coeffs << 1, -4, 6;
    gs_coeffs << 6, -3,
                -3,  2;

    for (int i=0; i<n; ++i) {
      if (i+2 < n) { A(i,i+2) = 1; }
      if (i+1 < n) { A(i,i+1) = -4; }
      A(i,i) = 6;
      if (i > 0) { A(i,i-1) = -4; }
      if (i > 1) { A(i,i-2) = 1; }
    }
  }else {
    coeffs.resize(1,2);
    gs_coeffs.resize(1,1);
  
    coeffs << -1, 2;
    gs_coeffs << 1;

    for (int i=0; i<n; ++i) {
      if (i+1 < n) { A(i,i+1) = -1; }
      A(i,i) = 2;
      if (i > 0) { A(i,i-1) = -1; }
    }
  }

  
  for ( int i = 0; i < gs_coeffs.rows() ; i ++ ){ 
    int A_row = A.rows() - gs_coeffs.rows() + i;

    for ( int j = 0; j < gs_coeffs.cols() ; j ++ ){ 
      int A_col = A.cols() - gs_coeffs.cols() + j;

      A(A_row, A_col) = gs_coeffs(i,j);
    }
  }


  MatX Ainv = A.inverse();

  MatX L, Ls;
  MatX testvec = MatX::Random(n, 2 );
  MatX testvec2 = testvec;
  MatX m1(n,n); m1.setIdentity();
  MatX m2(n,n); m2.setIdentity();
  MatX m3(n,n); m3.setIdentity();
  MatX m4(n,n); m4.setIdentity();
  MatX work;

  Eigen::LLT<MatX> cholSolver(A);
  m1 = cholSolver.solve(m1);


  skylineChol(n, coeffs, gs_coeffs, Ls);
  skylineCholSolve(Ls, m2);

  diagMul(coeffs, gs_coeffs, m4, m3);

    
  
  if (n < 20) { 
    std::cout << "A =\n" << A << "\n";
    std::cout << "A =\n" << m3 << "\n";
    std::cout << "A.inv() =\n" << A.inverse() << "\n";
    std::cout << "m1 =\n" << m1 << "\n";
    std::cout << "errm1 = \n" << m1-A.inverse() << "\n";
    std::cout << "m2 =\n" << m2 << "\n";
    std::cout << "errm2 = \n" << m2-A.inverse() << "\n";
  }
  
  assert( relErr(A, m3) < 1e-5 );
  assert( relErr(m1, m2) < 1e-5 );

  return 0;


}
int main(int argc, char** argv) { 
    //primary( argc, argv );
    //secondary( argc, argv );
    tertiary( argc, argv );
}


