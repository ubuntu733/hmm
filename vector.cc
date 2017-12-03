/*！
 * Copyright (c) 2017 by ChengEn Shen
 * \file vector.cc
 * \brief 
 * \author ChenEn Shen
 */
#include "vector.h"

vector::vector(int64_t dim) {
  dim_ = dim;
  data_ = new real[dim_];
}

vector::~vector() {
  delete[] data_;
}

int64_t vector::size() const {
  return dim_;
}

void vector::zero() {
  for (int64_t i=0; i<dim_; ++i) {
    data_[i] = 0;
  }
}

real& vector::operator[](int64_t index) {
  return data_[index];
}