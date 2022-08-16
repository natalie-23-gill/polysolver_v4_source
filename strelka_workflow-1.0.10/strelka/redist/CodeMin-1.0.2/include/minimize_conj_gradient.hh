// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// CodeMin : lightweight minimization library from the
//           CodeAxe phylogenetic analysis package
//
//   http://www.phrap.org
//
// Copyright 2007, 2009 Christopher T Saunders (ctsa@u.washington.edu)
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

#include "minimize_array_util.h"
#include "minimize_line.h"
#include "minimize_util.h"
#include "workspace.h"


namespace codemin {


template <typename FloatType,
          typename MinFunc>
void
minimize_conj_gradient(FloatType* vec_current,
                       MinFunc& mf,
                       const FloatType tol,
                       FloatType& f_min,
                       unsigned& iter,
                       FloatType& final_iter_delta_f,
                       const unsigned max_iter){

  const unsigned n(mf.dim());

  workspace<FloatType> ws(4*n);

  FloatType* vec_deriv(ws.ptr());
  FloatType* vec_grad(ws.ptr()+n);
  FloatType* vec_conj_grad(ws.ptr()+2*n);
  FloatType* linespace(ws.ptr()+3*n);

  FloatType gdot(0);

  std::fill(vec_conj_grad,vec_conj_grad+n,FloatType(0.));

  // more intuitive internal name:
  FloatType& iter_delta_f(final_iter_delta_f);

  iter=0;
  while(iter<max_iter){
    iter++;

    f_min = mf.dval(vec_current,vec_deriv);

    FloatType gnum(0);
    if(iter>1){
      for(unsigned i(0);i<n;++i){
        gnum += (vec_deriv[i]+vec_grad[i])*vec_deriv[i];
      }
      gnum /= gdot;
    }

    for(unsigned i(0);i<n;++i){
      vec_grad[i] = -vec_deriv[i];
      vec_conj_grad[i] = vec_grad[i]+gnum*vec_conj_grad[i];
    }

    const FloatType f_min_prev(f_min);

    try {
      minimize_line(vec_current,vec_conj_grad,f_min,mf,f_min,linespace);
    } catch (std::exception& e) {
      min_log_os() << "EXCEPTION: " << e.what() << "\n"
                   << "...caught in " << __FILE__ << ":" << __LINE__ << " gnum,gdot " << gnum << " " << gdot << "\n";
      for(unsigned i(0);i<n;++i){
        min_log_os() << "i, x, -grad(x), conj_grad: " << i << "\t"
                     << vec_current[i] << "\t" << vec_grad[i] << "\t" << vec_conj_grad[i] << "\n";
      }
      throw;
    }

    iter_delta_f=(f_min_prev-f_min);

    if(iter_delta_f <= tol) { // converging on absolute tolerance
      if(iter_delta_f<0.) throw minimize_exception("minimize_conj_grad(): value of f increased in iteration");
      break;
    }

    gdot = array_dot(vec_grad,vec_grad,n);
    if(gdot == 0.) break;
  }
}

}
