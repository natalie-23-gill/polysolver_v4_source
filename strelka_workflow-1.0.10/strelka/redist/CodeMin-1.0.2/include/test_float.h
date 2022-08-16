// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// CodeMin : lightweight minimization library from the
//           CodeAxe phylogenetic analysis package
//
//   http://www.phrap.org
//
// Copyright 2007, 2009, 2010 Christopher T Saunders (ctsa@u.washington.edu)
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
/// \brief portable floating point tests
///
/// This code uses simple floating point tests by default which may
/// not be accurate on all systems. For greater reliability, define
/// CODEMIN_USE_BOOST if boost is available, or CODEMINE_USE_MATH_H to
/// attempt to use the clib floating point tests on systems where
/// these tests are available.
/// 

#ifndef __TEST_FLOAT_H
#define __TEST_FLOAT_H

#ifdef CODEMIN_USE_BOOST
#include "boost/math/special_functions.hpp"
#else 
#ifdef CODEMIN_USE_MATH_H
#include <cmath>
#endif
#endif

namespace codemin {

#ifdef CODEMIN_USE_BOOST

template <typename FloatType>
bool
is_float_nan(FloatType x){ return (boost::math::isnan)(x); }

template <typename FloatType>
bool
is_float_inf(FloatType x){ return (boost::math::isinf)(x); }

#else
#ifdef CODEMIN_USE_MATH_H 

// use macros in math.h:

template <typename FloatType>
bool
is_float_nan(FloatType x){ return std::isnan(x); }

template <typename FloatType>
bool
is_float_inf(FloatType x){ return std::isinf(x); }

#else 

// portable tests, may be less reliable than boost or clib macros

// these tests not working??? check against -ffast-math (and friends) in gcc
//
// add this if sassert is available: STATIC_ASSERT( (0./0.) != (0./0.) )

template <typename FloatType>
bool
is_float_nan(FloatType x){ return x != x; }

template <typename FloatType>
bool
is_float_inf(FloatType x){ return (is_float_nan(x*0) and (! is_float_nan(x))); }

#endif
#endif


template <typename FloatType>
bool
is_float_invalid(FloatType x){
  return (is_float_nan(x) or is_float_inf(x));
}

}

#endif

