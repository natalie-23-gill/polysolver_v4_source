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

#ifndef __MINIMIZE_1D_H
#define __MINIMIZE_1D_H

//#include "minfunc_interface.h"

namespace codemin {

/// \brief Find local minimum in a 1d function
///
/// f1 must be set to mf.val(x1) at call
/// x2 is the first argument evaluated in the search for the minimum --
///    this is useful for suggesting the range of the 1d search
///
/// default minimization tolerance (x_tol) is sqrt(epsilon(FloatType))
/// as suggested by RP Brent in "Minimization without derivatives"
///
/// an optional secondary tolerance on the change to the function
/// value at each iteration may also be set
///
/// \param x1 starting function arg value
///
/// \param x2 next function arg to evaluate. abs(x1-x2) suggests the
/// scale of the domain to search over -- this can speed up the
/// minimization if known
///
/// \param f1 must be set to f(x1) [mf.val(x1)]
///
/// \param mf the minimization function. Must conform to the
/// interface of: minfunc_1d_interface<FloatType> >
///
/// \param x_min returns the argument of minimum function value
///
/// \param f_min returns the minimum function value
///
/// \param x_tol "normal" minimization tolerance: the minimum change
/// in x per iteration
///
/// \param f_tol secondary minimization tolerance: the minimum change
/// in f(x) per iteration, if both tolerances are set the routine
/// stops when either criteria is not met in a single iteration
///
template <typename FloatType,
          typename MinFunc> // = minfunc_1d_interface<FloatType> >
void
minimize_1d(FloatType x1,
            FloatType x2,
            FloatType f1,
            MinFunc& mf,
            FloatType& x_min,
            FloatType& f_min,
            FloatType x_tol=0,
            FloatType f_tol=0);

}

#include "minimize_1d.hh"


#endif
