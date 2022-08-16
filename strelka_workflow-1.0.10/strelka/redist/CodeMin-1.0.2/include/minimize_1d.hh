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

#include "minimize_util.h"
#include "test_float.h"

#include <cmath>

#include <algorithm>
#include <limits>

#ifdef CODEMIN_USE_BOOST
#include "boost/math/special_functions/sign.hpp"
using boost::math::copysign;
#endif

#ifdef CODEMIN_DEBUG
#define MIN_1D_DEBUG
#endif

namespace codemin {

template <typename FloatType,
          typename MinFunc>
void
get_computable_projection(FloatType& d_project,
                          const FloatType d_base,
                          MinFunc& mf,
                          const char* project_label,
                          const FloatType offset=0,
                          const FloatType scale=1){

  static const unsigned MAX_REDUCTION(10);
  static const FloatType REDUCTION_FACTOR(0.1);

  const FloatType orig_d_project(d_project);

  for(unsigned i(0);i<MAX_REDUCTION;++i){
    try{
      const bool ivc(mf.is_val_computable(offset+d_project*scale));
      if(ivc) return;
    } catch(std::exception& e){
      min_log_os() << "EXCEPTION: " << e.what() << "\n"
                   << "...caught in " << __FILE__ << ":" << __LINE__
                   << " i,o,dp,s: " << i << " " << offset << " " << d_project << " " << scale << "\n";
      throw;
    }

    d_project=d_base+(d_project-d_base)*REDUCTION_FACTOR;

    min_log_os() << "WARNING: adjusting " << project_label << " from/to: "
                 << orig_d_project << " " << d_project << "\n";
  }

  // ARG!!!
  //pass_away("can't find computable projection in minimize_1d");
}



// given points d1 and d2, and f1=f(d1), find a d1,d2,d3
// bracketed minimum
//
template <typename FloatType,
          typename MinFunc>
void
bracket_minimum(FloatType& d1,
                FloatType& d2,
                FloatType& d3,
                FloatType& f1,
                FloatType& f2,
                FloatType& f3,
                MinFunc& mf){

  static const FloatType GOLD_RATIO((3.-std::sqrt(5.))/2.);
  static const FloatType ADVANCE_RATIO((1.-GOLD_RATIO)/GOLD_RATIO);
  static const FloatType MAX_ADVANCE_RATIO(100.);

  get_computable_projection(d2,d1,mf,"d2");

  try {
    f2 = mf.val(d2);
  } catch(std::exception& e) {
    min_log_os() << "EXCEPTION: " << e.what() << "\n"
                 << "...caught in bracket_minimum, call to val(d2), d2: " << d2 << "\n";
    throw;
  }

#ifdef MIN_1D_DEBUG
  min_log_os() << "bracket_minimum, init: f2=val(d2), d2: " << f2 << " " << d2 << "\n";
#endif

  // setup so 1->2 is downhill; also, no use in continuing in orig
  // 1->2 direction if f(1) == f(2)
  if(f1 <= f2){
    std::swap(f1,f2);
    std::swap(d1,d2);
  }

  const FloatType offset(d1);
  const FloatType scale(d2-d1);

  FloatType wd1(0.),wd2(1.),wd3;

  // goal is to flip between golden ratio steps and parabolic steps
  // until we bracket a minimum
  //
  while(true){
    // golden ratio advance:
    wd3 = wd2+ADVANCE_RATIO*(wd2-wd1);

    get_computable_projection(wd3,wd2,mf,"wd3",offset,scale);

    try {
      f3 = mf.val(offset+wd3*scale);
    } catch(std::exception& e) {
      min_log_os() << "EXCEPTION: " << e.what() << "\n"
                   << "...caught in bracket_minimum, call to val(offset+wd3*scale), offset,wd3,scale,wd2,wd1: "
                   << offset << " " << wd3 << " " << scale  << " " << wd2 << " " << wd1 << "\n";
      throw;
    }
#ifdef MIN_1D_DEBUG
    min_log_os() << "bracket_minimum, f3=val(offset+wd3*scale), offset,wd3,scale,wd2,wd1: "
                 << f3 << " " << offset << " " << wd3 << " " << scale  << " " << wd2 << " " << wd1 << "\n";
#endif

    if(f2 <= f3) break;  // found bracketed minimum

    // from d1,d2,d3 formulate parabolic advance to d4:
    const FloatType r((wd2-wd1)*(f2-f3));
    const FloatType q((wd3-wd2)*(f2-f1));
    const FloatType wd4_max(wd2+MAX_ADVANCE_RATIO*(wd3-wd2));
    const FloatType wd4_tmp(wd2-((wd3-wd2)*q-(wd2-wd1)*r)/(2.*(q-r)));
    FloatType wd4(std::min(wd4_tmp,wd4_max));

    if(wd4 > wd2){ // parabolic advance is reasonable:
      if(wd4 > wd3){
        get_computable_projection(wd4,wd3,mf,"wd4",offset,scale);
      }

      FloatType f4;
      try {
        f4 = mf.val(offset+wd4*scale);
      } catch(std::exception& e) {
        min_log_os() << "EXCEPTION: " << e.what() << "\n"
                     << "...caught in bracket_minimum, call to val(offset+wd4*scale), offset,wd4,scale: "
                     << offset << " " << wd4 << " " << scale  << " "
                     << r << " " << q << " " << wd4_max << " " << wd4_tmp << " "
                     << wd3 << "\n";
        throw;
      }
#ifdef MIN_1D_DEBUG
      min_log_os() << "bracket_minimum, f4=val(offset+wd4*scale), offset,wd4,scale: "
                   << f4 << " " << offset << " " << wd4 << " " << scale  << " "
                   << r << " " << q << " " << wd4_max << " " << wd4_tmp << " "
                   << wd3 << "\n";
#endif

      if(wd4 < wd3){
        if(f4 < f3) { // found bracketed minimum!
          f1=f2; wd1=wd2;
          f2=f4; wd2=wd4;
          break;
        } else { // d1=d4,d2=d3, advance...
          f1=f4; wd1=wd4;
          f2=f3; wd2=wd4;
        }
      } else { // d1=d3,d2=d4, advance...
        if(f4 < f3){
          f1=f3; wd1=wd3;
          f2=f4; wd2=wd4;
        } else { // found  bracketed minimum!
          f1=f2; wd1=wd2;
          f2=f3; wd2=wd3;
          f3=f4; wd3=wd4;
          break;
        }
      }
    } else { // parabolic advance is unreasonable:
      f1=f2; wd1=wd2;
      f2=f3; wd2=wd3;
    }
  }

  d3 = offset+wd3*scale;
  d2 = offset+wd2*scale;
  d1 = offset+wd1*scale;
}




// start from bracketed interval d1,d2,d3, where, for fx=f(dx),
// f2 < f3 && f2 < f1
//
template <typename FloatType,
          typename MinFunc>
void
brent_line_minimizer(FloatType d1,
                     FloatType d2,
                     FloatType d3,
                     FloatType f1,
                     FloatType f2,
                     FloatType f3,
                     MinFunc& mf,
                     const FloatType dtol,
                     const FloatType ftol,
                     FloatType& dmin,
                     FloatType& fmin,
                     const unsigned max_iter = 100){

  static const FloatType GOLD_RATIO((3.-std::sqrt(5.))/2.);
  static const FloatType SEPS(1e-10);

  FloatType d_delta(0);
  FloatType d_delta_prev(0);
  FloatType d_low0(d2),d_low1(d1),d_low2(d3);
  FloatType f_low0(f2),f_low1(f1),f_low2(f3);

  if(f_low1>f_low2){
    std::swap(f_low1,f_low2);
    std::swap(d_low1,d_low2);
  }

  FloatType da(d1),db(d3);
  if(da>db) std::swap(da,db);

  { // setup d_delta_prev to start off with a parabolic fit:
    const FloatType d_middle(0.5*(da+db));
    if( d_low0 >= d_middle ) { d_delta_prev=da-d_low0; } // assign to larger interval between a and b
    else                     { d_delta_prev=db-d_low0; }
  }

  for(unsigned i(0);i<max_iter;++i){
    const FloatType d_middle(0.5*(da+db));
    const FloatType tol1(dtol*std::fabs(d_low0)+SEPS);
    if(2.*std::fabs(d_low0-d_middle) <= (4.*tol1-db+da)){
      dmin=d_low0;
      fmin=f_low0;
      return;
    }

    // add secondary tolerance on the change in f(x):
    if(std::fabs(f_low0-f_low2) <= ftol){
      dmin=d_low0;
      fmin=f_low0;
      return;
    }

    bool is_bad_fit(false);

    if(std::fabs(d_delta_prev) > tol1){  // solve for d_delta by parabolic fit
      const FloatType r((d_low0-d_low1)*(f_low0-f_low2));
      const FloatType q((d_low0-d_low2)*(f_low0-f_low1));
      if( q == r ) {
        is_bad_fit=true;
      } else {
        const FloatType p(-((d_low0-d_low2)*q-(d_low0-d_low1)*r)/(2.*(q-r)));

        const FloatType d_delta_prev2(d_delta_prev);
        d_delta_prev=d_delta;

        is_bad_fit = is_float_invalid(p) || std::fabs(p) >= std::fabs(0.5*d_delta_prev2) ||
          p >= (d_low0-da) || p >= (db-d_low0);

        if(! is_bad_fit){
          d_delta = p;
          const FloatType d_tmp = d_low0+d_delta;
          if((d_tmp-da) < 2.*tol1 || (db-d_tmp) < 2.*tol1) d_delta = copysign(tol1,d_middle-d_low0);
        }
      }
    } else { is_bad_fit=true; }

    if(is_bad_fit){ // solve for d_delta by golden ratio:
      if( d_low0 >= d_middle ) { d_delta_prev=da-d_low0; } // assign to larger interval between a and b
      else                     { d_delta_prev=db-d_low0; }
      d_delta = GOLD_RATIO*d_delta_prev;
    }

    d_delta=copysign(std::max(std::fabs(d_delta),tol1),d_delta);
    const FloatType d_next(d_low0+d_delta);

    FloatType f_next;
    try {
      f_next = mf.val(d_next);
    } catch(std::exception& e) {
      min_log_os() << "EXCEPTION: " << e.what() << "\n"
                   << "...caught in brent_line_minimizer, call to val(d_next), d_next,d_low0,d_delta: "
                   << d_next << " " << d_low0 << " " << d_delta << " "
                   << d_low1 << " " << d_low2 << " "
                   << f_low0 << " " << f_low1 << " " << f_low2 << " "
                   << da << " " << db << " " << d_delta_prev << " " << is_bad_fit << "\n";
      throw;
    }
#ifdef MIN_1D_DEBUG
    min_log_os() << "brent_line_minimizer, f_next=val(d_next), d_next,d_low0,d_delta: "
                 << f_next << " " << d_next << " " << d_low0 << " " << d_delta << " "
                 << d_low1 << " " << d_low2 << " "
                 << f_low0 << " " << f_low1 << " " << f_low2 << " "
                 << da << " " << db << " " << d_delta_prev << " " << is_bad_fit << "\n";
#endif

    // update bracket:
    if(f_next <= f_low0){
      if(d_next >= d_low0) da=d_low0;
      else                 db=d_low0;
    } else {
      if(d_next < d_low0)  da=d_next;
      else                 db=d_next;
    }

    // update low values:
    if       (f_next <= f_low0){
      d_low2 = d_low1; f_low2 = f_low1;
      d_low1 = d_low0; f_low1 = f_low0;
      d_low0 = d_next; f_low0 = f_next;
    } else if(f_next <= f_low1){
      d_low2 = d_low1; f_low2 = f_low1;
      d_low1 = d_next; f_low1 = f_next;
    } else if(f_next <= f_low2){
      d_low2 = d_next; f_low2 = f_next;
    }
  }
  throw minimize_exception("brent line minimizer failed to converge");
}




template <typename FloatType,
          typename MinFunc>
void
minimize_1d(FloatType x1,
            FloatType x2,
            FloatType f1,
            MinFunc& mf,
            FloatType& x_min,
            FloatType& f_min,
            FloatType x_tol,
            FloatType f_tol){

  static const FloatType default_x_tol(std::sqrt(std::numeric_limits<FloatType>::epsilon()));

  x_tol = std::max(x_tol,default_x_tol);

  FloatType x3;
  FloatType f2,f3;

  bracket_minimum(x1,x2,x3,f1,f2,f3,mf);
  brent_line_minimizer(x1,x2,x3,f1,f2,f3,mf,x_tol,f_tol,x_min,f_min);
}

}
