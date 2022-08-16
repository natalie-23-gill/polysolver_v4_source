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

#ifndef __MINIMIZE_LINE_H
#define __MINIMIZE_LINE_H


//#include "minfunc_interface.h"

namespace codemin {


/// \brief minimize function along a line in a multidimensional parameter space
///
/// \param vec_base min line defined by vec_base+C*vec_dir; returns minimum
/// function argument
///
/// \param vec_dir min line defined by vec_base+C*vec_dir; length of
/// vec_dir gives minimizer a hint of the scale to search over
///
/// \param f_base must be set to f(vec_base) [ mf.val(vec_base ]
///
/// \param mf minimization function; must fulfill
/// minfunc_interface<FloatType>
///
/// \param f_min minimum function value
///
/// \param scratch_space tmp storage for routine: size=mf.dim()*sizeof(FloatType)
///
/// \param x_tol minimum change in argument value per iteration
///
/// \param f_tol minimum change to function value per iteration
///
template <typename FloatType,
          typename MinFunc> // = minfunc_interface<FloatType> >
void
minimize_line(FloatType* vec_base,
              FloatType* vec_dir,
              FloatType f_base,
              MinFunc& mf,
              FloatType& f_min,
              FloatType* scratch_space,
              FloatType x_tol=0,
              FloatType f_tol=0);

}

#include "minimize_line.hh"

#endif
