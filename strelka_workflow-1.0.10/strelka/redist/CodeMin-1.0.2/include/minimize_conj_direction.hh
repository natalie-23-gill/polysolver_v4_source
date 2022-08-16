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

#include <limits>

namespace codemin {

#ifdef CODEMIN_DEBUG
#define MIN_CJ_DEBUG
#endif



template <typename FloatType,
          typename MinFunc>
void
minimize_conj_direction(FloatType* vec_current,
                        FloatType* conj_direction_set,
                        MinFunc& mf,
                        FloatType& start_tol,
                        const FloatType end_tol,
                        const FloatType line_tol,
                        FloatType& f_min,
                        unsigned& iter,
                        FloatType& final_iter_delta_f,
                        const unsigned max_iter){

  static const FloatType ZERO_VECTOR_NORM_THRESHOLD(std::numeric_limits<FloatType>::epsilon()*1e2);

  const unsigned n(mf.dim());

  workspace<FloatType> ws(3*n);
  FloatType* vec_start(ws.ptr());
  FloatType* vec_newdir(ws.ptr()+n);
  FloatType* linespace(ws.ptr()+2*n);

  FloatType tol(start_tol);

  f_min = mf.val(vec_current);

  // rename this so that it's less confusing within the function
  FloatType& iter_delta_f(final_iter_delta_f);

  iter=0;
  while(iter<max_iter){
    iter++;

    const FloatType f_min_start(f_min);
    std::copy(vec_current,vec_current+n,vec_start);

    unsigned max_drop_index(0);
    FloatType max_drop(0.);

    for(unsigned i(0);i<n;++i){
      const FloatType f_min_prev(f_min);
      try{
        minimize_line(vec_current,conj_direction_set+i*n,f_min,mf,f_min,linespace,line_tol,tol);
      } catch (std::exception& e){
        min_log_os() << "EXCEPTION: " << e.what() << "\n"
                     << "...caught in " << __FILE__ << ":" << __LINE__ << "\n";
        throw;
      }

      const FloatType f_drop(f_min_prev-f_min);
      if(f_drop > max_drop){
        max_drop = f_drop;
        max_drop_index = i;
      }
#ifdef MIN_CJ_DEBUG
      min_log_os() << "CONJD: i,f_min_prev,f_min: " << i << " " << f_min_prev << " " << f_min << "\n";
#endif
    }

    for(unsigned i(0);i<n;++i){
      vec_newdir[i]  = vec_current[i] - vec_start[i];
    }

    if(array_norm2(vec_newdir,n) >= ZERO_VECTOR_NORM_THRESHOLD){

      FloatType* vec_project(linespace);
      for(unsigned i(0);i<n;++i){
        vec_project[i] = vec_newdir[i] + vec_current[i];
      }

      const FloatType f_min_project(mf.val(vec_project));

      if(f_min_project < f_min_start){
        const FloatType tmp1(f_min_start-f_min-max_drop);
        const FloatType tmp2(f_min_start-f_min_project);
        const FloatType test(2.*(f_min_start-2.*f_min+f_min_project)*tmp1*tmp1-max_drop*tmp2*tmp2);

        if(test < 0.){
          for(unsigned i(0);i<n;++i){
            conj_direction_set[i+max_drop_index*n] = conj_direction_set[i+(n-1)*n];
            conj_direction_set[i+(n-1)*n] = vec_newdir[i];
          }
        }

        if(f_min_project < f_min){ // don't waste f(vec_project) if it's lower:
          std::copy(vec_project,vec_project+n,vec_current);
          f_min=f_min_project;
        }
      }

      try {
        minimize_line(vec_current,vec_newdir,f_min,mf,f_min,linespace,line_tol,tol); // try vec_newdir no matter what test says
      } catch (std::exception& e){
        min_log_os() << "EXCEPTION: " << e.what() << "\n"
                     << "...caught in " << __FILE__ << ":" << __LINE__ << "\n";
        throw;
      }

#ifdef MIN_CJ_DEBUG
    } else {
      min_log_os() << "CONJD: fmin,fmin_start: " << f_min << " " << f_min_start << "\n"
                   << "CONJD: tol,start_tol,end_tol: " << tol << " " << start_tol << " " << end_tol << "\n"
                   << "CONJD: anorm2: " << array_norm2(vec_newdir,n) << "\n";
      throw minimize_exception("minimize_conj_dir(): vec_newdir zero norm\n");
#endif
    }

#ifdef MIN_CJ_DEBUG
    min_log_os() << "CONJD: iter, tol,end_tol fabs(f_min-f_min_start): "
                 << iter << " " << tol << " " << end_tol << " " << std::fabs(f_min-f_min_start) << "\n";

#endif

    iter_delta_f=(f_min_start-f_min);

    if(iter_delta_f <= tol) { // converging on absolute tolerance
      if(iter_delta_f<0.) throw minimize_exception("minimize_conj_dir(): value of f increased in iteration");
      if(tol*0.999<=end_tol) break;
      tol *= .1;
      if((tol*0.9)<=end_tol) tol=end_tol; // tol <- end_tol as soon as tol gets reasonably close
    }
  }
}

}
