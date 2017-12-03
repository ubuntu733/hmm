/*！
 * Copyright (c) 2017 by ChengEn Shen
 * \file model.h
 * \brief hmm模型
 * \author ChenEn Shen
 */
#ifndef HMM_MODEL_H
#define HMM_MODEL_H

#include <cstdio>
#include <vector>
#include <iostream>

#include "matrix.h"
#include "vector.h"

class model {
 public:
  model(int64_t, int64_t);
  model(int64_t n,
        int64_t m,
        matrixPtr transitionMatrix,
        matrixPtr emitMatrix,
        vectorPtr initStates);

  void load(FILE *fp);
  void save(FILE *fp);

  void initModel();

  real viterbi(const std::vector<int64_t>& observation, //观察序列
                std::vector<int64_t>& states);  //隐态序列
  real Forward(const std::vector<int64_t>& observation);
  void BaumWelch(const std::vector<int64_t> &observation);

 private:
  void Forward(const std::vector<int64_t>& observation,
                matrixPtr& alpha, real& probs);
  void Backward(const std::vector<int64_t>& observation,
                matrixPtr& belta);
  void ComputeGamma(matrixPtr& gamma, const matrixPtr& alpha, const matrixPtr& belta);
  void ComputeXi(const std::vector<int64_t>& observation,
                 std::vector<matrixPtr> &xi,
                 const matrixPtr& alpha,
                 const matrixPtr& belta);
  std::vector<matrixPtr> allocXi(int64_t T, int64_t N);
 private:
  int64_t N_;  //隐态数量
  int64_t M_;  //观察态数量
  matrixPtr A_;  //转移概率
  matrixPtr B_;  //发射概率
  vectorPtr pi_;  //初始概率

};


#endif //HMM_MODEL_H
