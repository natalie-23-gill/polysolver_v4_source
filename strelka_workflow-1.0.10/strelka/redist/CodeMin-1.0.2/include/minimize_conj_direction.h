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

#ifndef __MINIMIZE_CONJ_DIRECTION_H
#define __MINIMIZE_CONJ_DIRECTION_H

//#include "minfunc_interface.h"

namespace codemin {

/// \brief minimize function using modified conjugate direction set
/// method from (Powell, 1964)
///
/// \param vec_current starting function argument; on exit set to
/// minimum function argument
///
/// \param conj_direction_set starting conjugate direction set, scale
/// of conjugate vectors provides a hint to the minimizer of the scale
/// of parameter space to explore
///
/// \param mf minimization function; must fulfill
/// minfunc_interface<FloatType>
///
/// \param start_tol starting tolerance (on min function value
/// improvement in each iteration) this will be decreased by one order
/// of magnitude at each iteration where tolerance is not exceeded
/// until it reaches (or falls below) end_tol, at which point
/// minimization will terminate when end_tol is not exceeded. This is
/// a heuristic speed-up for long minimizations.
///
/// \param end_tol final tolerance (on min function value improvement
/// in each iteration)
///
/// \param line_tol minimum change in argument value per iteration during
/// line minimization
///
/// \param f_min final minimum function value
///
/// \param iter number of iterations required to reach minimum
///
/// \param final_iter_delta_f function improvement in final iteration
///
/// \param max_iter maximum number of conj dir iterations to attempt
///
template <typename FloatType,
          typename MinFunc> // = minfunc_interface<FloatType>
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
                        const unsigned max_iter);

}

#include "minimize_conj_direction.hh"

#endif
