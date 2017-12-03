/*！
 * Copyright (c) 2017 by ChengEn Shen
 * \file matrix.cc
 * \brief 
 * \author ChenEn Shen
 */
#include "matrix.h"

matrix::matrix() {
  m_ = 0;
  n_ = 0;
  data_ = nullptr;
}

matrix::matrix(int64_t m, int64_t n) {
  m_ = m;
  n_ = n;
  data_ = new real[m*n];
}

matrix::~matrix() {
  delete[] data_;
}

void matrix::set(int64_t i, int64_t j, real value) {
  data_[i*n_+j] = value;
}

real& matrix::at(int64_t i, int64_t j) {
  return data_[i*n_+j];
}

void matrix::zero() {
  for (int64_t i=0; i<m_*n_; ++i) {
    data_[i*n_+m_] = 0;
  }
}

int64_t matrix::row() {
  return m_;
}

int64_t matrix::col() {
  return n_;
}