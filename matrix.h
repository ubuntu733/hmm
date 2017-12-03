/*！
 * Copyright (c) 2017 by ChengEn Shen
 * \file matrix.h
 * \brief 
 * \author ChenEn Shen
 */
#ifndef HMM_MATRIX_H
#define HMM_MATRIX_H

#include <cstdint>
#include <memory>
#include "real.h"

class matrix {
 private:
  real *data_;
  int64_t m_;
  int64_t n_;

 public:
  matrix();
  matrix(int64_t m, int64_t n);
  ~matrix();

  real& at(int64_t i, int64_t j);

  void set(int64_t i, int64_t j, real value);

  void zero();

  int64_t row();
  int64_t col();
};

typedef std::shared_ptr<matrix> matrixPtr;


#endif //HMM_MATRIX_H
