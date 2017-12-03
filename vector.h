/*！
 * Copyright (c) 2017 by ChengEn Shen
 * \file vector.h
 * \brief 
 * \author ChenEn Shen
 */
#ifndef HMM_VECTOR_H
#define HMM_VECTOR_H

#include <cstdint>
#include <memory>
#include "real.h"

class vector {
 public:
  explicit vector(int64_t dim);
  ~vector();

  real& operator[](int64_t index);

  int64_t size() const;

  void zero();

 private:
  int64_t dim_;
  real *data_;
};

typedef std::shared_ptr<vector> vectorPtr;

#endif //HMM_VECTOR_H
