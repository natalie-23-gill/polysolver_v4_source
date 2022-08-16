// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// CodeMin : lightweight minimization library from the
//           CodeAxe phylogenetic analysis package
//
//   http://www.phrap.org
//
// Copyright 2009 Christopher T Saunders (ctsa@u.washington.edu)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

/// \file


#include <minfunc_interface.h>
#include <minimize_conj_direction.h>

#include <iostream>
#include <vector>


const unsigned N_DIM(5);


struct test_minfunc : public codemin::minfunc_interface<double> {

  test_minfunc(const std::vector<double>& a) : _a(a) {}

  unsigned dim() const { return N_DIM; }

  double val(const double* x) {
    const unsigned nv(dim());

    double val(0.);
    for(unsigned i(0);(i+1)<nv;++i){
      const double x1(x[i]-_a[i]);
      const double t1((x[i+1]-_a[i+1])-(x1*x1));
      const double t2(1.-(x1));
      val += 100.*t1*t1+t2*t2;
    }
    return val;
  }

private:
  const std::vector<double> _a;
};



int
main() {

  static const double line_tol(1e-7);
  static const double start_ratio(0.05);
  static const double min_start_dist(1e-6);
  static const double end_tol(1e-7);
  static const unsigned max_iter(200);

  double x[N_DIM] = {2.0, 4.0, -11.0, -0.9, 0.2};
  
  static const unsigned N_DIM2(N_DIM*N_DIM);
  double conj_dir[N_DIM2];
  std::fill(conj_dir,conj_dir+N_DIM2,0.);
  for(unsigned i(0);i<N_DIM;++i) {
    const double start_dir( std::max(std::fabs(x[i]*start_ratio),min_start_dist) );
    conj_dir[i*(N_DIM+1)] = start_dir;
  }

  std::vector<double> offsets(N_DIM);
  for(unsigned i(0);i<N_DIM;++i){ offsets[i] = i; }

  double start_tol(end_tol);
  unsigned iter;
  double func_min;
  double final_iter_delta;
  test_minfunc mf(offsets);

  codemin::minimize_conj_direction(x,conj_dir,mf,start_tol,end_tol,line_tol,
                                   func_min,iter,final_iter_delta,max_iter);

  std::cout << "\nminimize_conj_direction:\n\n";
  std::cout << "f(x): " << func_min << "\n\n";
  for(unsigned i(0);i<N_DIM;++i){
    std::cout << "x(" << i << "): " << x[i] << "\n";
  }
  std::cout << "\n";
}
