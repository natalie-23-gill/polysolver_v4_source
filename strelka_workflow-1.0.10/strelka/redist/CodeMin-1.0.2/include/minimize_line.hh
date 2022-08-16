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

#include "minfunc_interface.h"
#include "minimize_1d.h"

#include <algorithm>


namespace codemin {


template <typename FloatType>
void
array_scale_plus_vector(const FloatType* v1,
                        const FloatType s,
                        FloatType* v2,
                        const unsigned N){
  for(unsigned i(0);i<N;++i) v2[i] += s*v1[i];
}



/// \brief wraps a multidimensional function to expose a 1d interface
/// for line minimization
///
template <typename FloatType,
          typename MinFunc>
struct minfunc_1d_wrapper : public minfunc_1d_interface<FloatType> {

  minfunc_1d_wrapper(MinFunc& init_mf,
                     FloatType* scratch_space, // 1*mf.dim()
                     const FloatType* init_v_base,
                     const FloatType* init_v_dir,
                     const FloatType init_offset = 0)
    : n(init_mf.dim()), v_base(init_v_base), v_dir(init_v_dir),
      v_tmp(scratch_space), mf(init_mf), offset(init_offset) {}

  bool
  is_val_computable(const FloatType x) const {
    setup_vt(x);
    return mf.is_val_computable(v_tmp);
  }

  FloatType
  val(const FloatType x) const {
    setup_vt(x);

    FloatType a;
    try {
      a = mf.val(v_tmp);
    } catch (std::exception& e){
      min_log_os() << "EXCEPTION:: minfunc_1d_wrapper, call to val(v_tmp)\n";
      for(unsigned i(0);i<n;++i){
        min_log_os() << " i,v_tmp[i]: " << i << " " << v_tmp[i] << "\n";
      }
      throw;
    }
    return a;
  }

private:
  void
  setup_vt(const FloatType x) const {
    std::copy(v_base,v_base+n,v_tmp);
    array_scale_plus_vector(v_dir,(x-offset),v_tmp,n);
 }

public:

  unsigned n;
  const FloatType* v_base;
  const FloatType* v_dir;
  FloatType* v_tmp;
  MinFunc& mf;
  FloatType offset;  ///<< here to prevent dealing with all d hanging around 0 when initial guess is accurate
};




template <typename FloatType,
          typename MinFunc>
void
minimize_line(FloatType* vec_base,
              FloatType* vec_dir,
              FloatType f_base,
              MinFunc& mf,
              FloatType& f_min,
              FloatType* scratch_space,
              FloatType x_tol,
              FloatType f_tol){

  static const FloatType x_offset(1);
  static const FloatType x1(x_offset), x2(x_offset+1);

  minfunc_1d_wrapper<FloatType,MinFunc> mf1d(mf,scratch_space,vec_base,vec_dir,x_offset);

  FloatType x_min;

  minimize_1d(x1,x2,f_base,mf1d,x_min,f_min,x_tol,f_tol);

  x_min -= x_offset;

  array_scale_plus_vector(vec_dir,x_min,vec_base,mf.dim());
}

}
