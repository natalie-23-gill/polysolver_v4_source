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
///

#ifndef __MINIMIZE_ARRAY_UTIL_H
#define __MINIMIZE_ARRAY_UTIL_H


#include <cmath>

namespace codemin {


template <typename FloatType>
inline
FloatType
array_dot(const FloatType v1[],
          const FloatType v2[],
          const unsigned N){
  FloatType d = 0;
  for(unsigned i(0);i<N;++i){
    d += v1[i]*v2[i];
  }
  return d;
}


#ifdef USE_CBLAS
template <>
inline
float
array_dot(const float v1[],
          const float v2[],
          const unsigned N){
  return cblas_sdot(N,v1,1,v2,1);
}

template <>
inline
double
array_dot(const double v1[],
          const double v2[],
          const unsigned N){
  return cblas_ddot(N,v1,1,v2,1);

}
#endif


template <typename FloatType>
inline
FloatType
array_norm2(const FloatType* v,
            const unsigned N){
  FloatType sum(0);
  for(unsigned i(0);i<N;++i) sum += v[i]*v[i];
  return std::sqrt(sum);
}

#ifdef USE_CBLAS
template <>
inline
float
array_norm2(const float* v,
            const unsigned N){
  return cblas_snrm2(N,v,1);
}

template <>
inline
double
array_norm2(const double* v,
            const unsigned N){
  return cblas_dnrm2(N,v,1);
}
#endif

}

#endif

