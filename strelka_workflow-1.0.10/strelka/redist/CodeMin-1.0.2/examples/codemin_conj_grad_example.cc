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
#include <minimize_conj_gradient.h>

#include <iostream>
#include <vector>


const unsigned N_DIM(2);


struct test_minfunc : public codemin::minfunc_gradient_interface<double> {

  test_minfunc(const std::vector<double>& a) : _a(a) {}

  unsigned dim() const { return N_DIM; }

  double val(const double* x) {

    const double x0(x[0]-_a[0]);
    const double t1((x[1]-_a[1])-x0);
    const double t2(1.-x0);

    return 100.*t1*t1+t2*t2;
  }

  double dval(const double* x,
              double* dx) {

    const double x0(x[0]-_a[0]);
    const double t1((x[1]-_a[1])-x0);
    const double t2(1.-x0);

    dx[0] = -400.*t1*x0-2*t2;
    dx[1] = 200*t1;
    
    return 100.*t1*t1+t2*t2;
  }

private:
  const std::vector<double> _a;
};



int
main() {

  static const double tol(1e-7);
  static const unsigned max_iter(200);

  double x[N_DIM] = {-0.9, 0.2};
  
  std::vector<double> offsets(N_DIM);
  for(unsigned i(0);i<N_DIM;++i){ offsets[i] = i; }

  unsigned iter;
  double func_min;
  double final_iter_delta;
  test_minfunc mf(offsets);

  codemin::minimize_conj_gradient(x,mf,tol,func_min,iter,final_iter_delta,max_iter);

  std::cout << "\nminimize_conj_gradient:\n\n";
  std::cout << "f(x): " << func_min << "\n\n";
  for(unsigned i(0);i<N_DIM;++i){
    std::cout << "x(" << i << "): " << x[i] << "\n";
  }
  std::cout << "\n";
}
