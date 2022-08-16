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

#ifndef __WORKSPACE_H
#define __WORKSPACE_H


namespace codemin {


struct uncopyable {
protected:
    uncopyable() {}
    ~uncopyable() {}

private:
    uncopyable(const uncopyable&);
    uncopyable& operator=(const uncopyable&);
};



template <typename T>
struct workspace : private uncopyable {

  explicit
  workspace(const unsigned n = 0) : _ptr(0), _size(0) {
    if(n) resize(n);
  }

  ~workspace() { if(_ptr) delete [] _ptr; }

  void
  resize(const unsigned s){
    if(_size<s){
      _size=s;
      if(_ptr) delete [] _ptr;
      _ptr = new T[_size];
    }
  }

  void
  resize_with_data(const unsigned s){
    if(_size<s){
      T* new_ptr(new T[s]);
      std::copy(_ptr,_ptr+_size,new_ptr);
      _size=s;
      if(_ptr) delete [] _ptr;
      _ptr = new_ptr;
    }
  }

  T* ptr() { return _ptr; }

private:
  T* _ptr;
  unsigned _size;
};

}

#endif
