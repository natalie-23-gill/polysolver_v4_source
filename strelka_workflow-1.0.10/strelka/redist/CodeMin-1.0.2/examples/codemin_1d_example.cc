// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// CodeMin : lightweight minimization library from the
//           CodeAxe phylogenetic analysis package
//
//   http://www.phrap.org
//
// Copyright 2009 Christopher T Saunders (ctsa@u.washington.edu)
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


#include <minfunc_interface.h>
#include <minimize_1d.h>

#include <iostream>


struct test_minfunc : public codemin::minfunc_1d_interface<double> {

  test_minfunc(const double a) : _a(a) {}

  double val(const double x) const {
      return ((x-_a)*(x-_a));
  }

private:
  const double _a;
};



int
main() {

  static const double x_min_target(3.2);
  static const double x1(0.5);
  static const double x2(0.4);

  test_minfunc mf(x_min_target);

  double x_min;
  double func_min;

  codemin::minimize_1d(x1,x2,mf.val(x1),mf,x_min,func_min);

  std::cout << "\nminimize_1d:\n\n";
  std::cout << "f(x): " << func_min << "\n\n";
  std::cout << "x: " << x_min << "\n\n";
}
