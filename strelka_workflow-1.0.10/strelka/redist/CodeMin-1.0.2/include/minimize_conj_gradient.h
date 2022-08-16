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

#ifndef __MINIMIZE_CONJ_GRADIENT_H
#define __MINIMIZE_CONJ_GRADIENT_H

#include "minfunc_interface.h"

namespace codemin {


/// \brief polak-ribiere conjugate gradient minimizer
///
/// \param vec_current starting function argument; on exit set to
/// minimum function argument
///
/// \param mf minimization function; must fulfill
/// minfunc_interface<FloatType>
///
/// \param tol final tolerance (on min function value improvement in
/// each iteration), minimization terminates when min function
/// improvement does not exceed this value in one iteration
///
/// \param f_min final minimum function value
///
/// \param iter number of iterations required to reach minimum
///
/// \param final_iter_delta_f decrease in f_min in the final iteration
///
/// \param max_iter maximum number of conj dir iterations to attempt
///
template <typename FloatType,
          typename MinFunc> // = minfunc_gradient_interface<FloatType> >
void
minimize_conj_gradient(FloatType* vec_current,
                       MinFunc& mf,
                       const FloatType tol,
                       FloatType& f_min,
                       unsigned& iter,
                       FloatType& final_iter_delta_f,
                       const unsigned max_iter);

}

#include "minimize_conj_gradient.hh"

#endif
