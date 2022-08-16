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

#ifndef __MINFUNC_INTERFACE_H
#define __MINFUNC_INTERFACE_H

namespace codemin {

/// \brief min function interface for line minimizers
///
template <typename FloatType>
struct minfunc_1d_interface {
  virtual ~minfunc_1d_interface() {}

  virtual FloatType val(const FloatType x) const  = 0;

  virtual bool is_val_computable(const FloatType) const { return true; }
};


/// \brief min function interface for multidimensional minimizers
///
template <typename FloatType>
struct minfunc_interface {
  virtual ~minfunc_interface() {}

  virtual unsigned dim() const = 0;

  virtual FloatType val(const FloatType* v) = 0;

  // this is not guaranteed to be called before val(), it only gives
  // minimizers the option to handle out-of-bounds parameter space
  // more gracefully than would be possible with a penalized score
  //
  virtual bool is_val_computable(const FloatType*) { return true; }
};


template <typename FloatType>
struct minfunc_gradient_interface : public minfunc_interface<FloatType> {

  virtual
  FloatType dval(const FloatType* v,
                 FloatType* dv) = 0;
};

}


#endif
