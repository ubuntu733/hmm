/*！
 * Copyright (c) 2017 by ChengEn Shen
 * \file model.cc
 * \brief 
 * \author ChenEn Shen
 */
#include "model.h"

model::model(int64_t m, int64_t n):M_(m), N_(n) {
  A_.reset(new matrix(N_, N_));
  B_.reset(new matrix(N_, M_));
  pi_.reset(new vector(N_));
}

model::model(int64_t n,
             int64_t m,
             matrixPtr transitionMatrix,
             matrixPtr emitMatrix,
             vectorPtr initStates):N_(n), M_(m), A_(transitionMatrix), B_(emitMatrix), pi_(initStates) {

}

void model::initModel() {
  A_->zero();
  B_->zero();
  pi_->zero();
}

void model::BaumWelch(const std::vector<int64_t> &observation) {
  if (observation.empty()) {
    std::cerr << "观察序列为空" << "\n";
    return;
  }
  int64_t i,j,k;
  int64_t t,l;
  int64_t T = static_cast<int64_t >(observation.size());
  matrixPtr alpha(new matrix(T, N_));
  matrixPtr belta(new matrix(T, N_));
  matrixPtr gamma(new matrix(T, N_));
  std::vector<matrixPtr> xi = allocXi(T, N_);
  real probs, preprobs,delta;
  real numeratorA, denominatorA;
  real numeratorB, denominatorB;
  Forward(observation,alpha, probs);
  Backward(observation, belta);
  ComputeGamma(gamma, alpha, belta);
  ComputeXi(observation, xi, alpha, belta);
  do {
    for (i=0; i<N_; ++i) {
      pi_->operator[](i) = 0.001 + 0.999*gamma->at(0,i);
    }

    for (i=0; i<N_; ++i) {
      denominatorA = 0.0;
      for (t=0; t<T; ++t) {
        denominatorA += gamma->at(t,i);
      }

      for (j=0; j<N_; ++j) {
        numeratorA = 0.0;
        for (t=0; t<T; ++t) {
          numeratorA += gamma->at(t,j);
        }
        A_->set(i,j, 0.001 + 0.999*numeratorA/denominatorA);
      }
      denominatorB = denominatorA + gamma->at(T-1, i);
      for (k = 0; k < M_; k++) {
        numeratorB = 0.0;
        for (t = 0; t < T; t++) {
          if (observation[t] == k)
            numeratorB += gamma->at(t,i);
        }

        B_->set(i, k, 0.001+0.999*numeratorB/denominatorB);
      }
      delta = probs-preprobs;
      preprobs = preprobs;
    }
    Forward(observation,alpha, probs);
    Backward(observation,belta);
    ComputeGamma(gamma, alpha, belta);
    ComputeXi(observation, xi, alpha, belta);
  } while (probs-preprobs > delta);
}

void model::ComputeGamma(matrixPtr &gamma, const matrixPtr &alpha, const matrixPtr &belta) {
  int64_t i,j;
  int64_t t;
  real sum;

  int64_t T = static_cast<int64_t>(alpha->row());
  for (t=0; t<T; ++t) {
    sum = 0.0;
    for (j=0; j<N_; ++j) {
      gamma->set(t,j,alpha->at(t,j)*belta->at(t,j));
      sum += gamma->at(t,j);
    }

    for (i=0; i<N_; ++i) {
      gamma->set(t, i, gamma->at(t,i)/sum);
    }
  }
}

void model::ComputeXi(const std::vector<int64_t>& observation,
                      std::vector<matrixPtr> &xi,
                      const matrixPtr &alpha,
                      const matrixPtr &belta) {
  int64_t i,j,t;
  real sum;
  int64_t T = static_cast<int64_t >(observation.size());
  for (t=0; t<T; ++t) {
    sum = 0.0;
    for (i=0; i<N_; ++i) {
      for (j=0; j<N_; ++j) {
        xi[t]->set(i, j, alpha->at(t,i)*belta->at(t+1,j)*A_->at(i,j)*B_->at(j,observation[t+1]));
        sum += xi[t]->at(i,j);
      }
    }
    for (i=0; i<N_; ++i) {
      for (j = 0; j <N_ ; ++j) {
        real tmp = xi[t]->at(i,j);
        xi[t]->set(i,j, tmp/sum);
      }
    }
  }
}

void model::Forward(const std::vector<int64_t> &observation, matrixPtr& alpha, real &probs) {
  if (observation.empty()) {
    std::cerr << "观察序列为空" << "\n";
    return;
  }
  alpha->zero();
  int64_t i,j;
  int64_t T = static_cast<int64_t >(observation.size());
  for (i=0; i<N_; ++i) {
    alpha->set(0, i, pi_->operator[](i)*B_->at(i, observation[0]));
  }

  for (int64_t t=0; t<T-1; ++t) {
    for (i=0; i<N_; ++i) {
      real sumprobs = 0.0;
      for (j=0; j<N_; ++j) {
        sumprobs += alpha->at(t,j)*A_->at(j,i);
      }
      alpha->set(t+1,i, sumprobs*B_->at(i, observation[t+1]));
    }
  }

  probs = 0.0;
  for (i=0; i<N_; i++) {
    probs += alpha->at(T-1,i);
  }
  return;
}


void model::Backward(const std::vector<int64_t> &observation, matrixPtr &belta) {
  if (observation.empty()) {
    std::cerr << "观察序列为空" << "\n";
    return;
  }
  belta->zero();
  int64_t i,j;
  int64_t T = static_cast<int64_t >(observation.size());
  for (int64_t i=0; i<N_; ++i) {
    belta->set(T-1, i, 1.0);
  }

  for (int64_t t=T-2; t>=0; --t) {
    for (int64_t i=0; i<N_; ++i) {
      real sum = 0.0;
      for (int64_t j=0; j<N_; ++j) {
        sum += A_->at(i,j)*B_->at(j,observation[t+1])*belta->at(t+1, j);
      }
      belta->set(t, i, sum);
    }
  }
}

real model::Forward(const std::vector<int64_t>& observation) {
  if (observation.empty()) {
    std::cerr << "观察序列为空" << "\n";
    return 0.0;
  }
  int64_t i,j;
  int64_t T = static_cast<int64_t >(observation.size());
  matrixPtr alpha(new matrix(T, N_));
  for (i=0; i<N_; ++i) {
    alpha->set(0, i, pi_->operator[](i)*B_->at(i, observation[0]));
  }

  real sumprobs;
  for (int64_t t=0; t<T-1; ++t) {
    for (i=0; i<N_; ++i) {
      sumprobs = 0.0;
      for (j=0; j<N_; ++j) {
        sumprobs += alpha->at(t,j)*A_->at(j,i);
      }
      alpha->set(t+1,i, sumprobs*B_->at(i, observation[t+1]));
    }
  }
  sumprobs = 0.0;
  for (i=0; i<N_; i++) {
    sumprobs += alpha->at(T-1,i);
  }
  return  sumprobs;
}

real model::viterbi(const std::vector<int64_t> &observation, std::vector<int64_t> &states) {
  if (observation.empty()) {
    std::cerr<< "观察序列为空" << "\n";
    return 0.0;
  }
  int64_t i,j;
  int64_t	maxvalind;
  real	maxval, val;

  int64_t T = static_cast<int64_t>(observation.size());
  matrixPtr delta(new matrix(T, N_));
  matrixPtr psi(new matrix(T, N_));

  for (i=0; i<N_; ++i) {
    delta->set(0, i, pi_->operator[](i)*B_->at(i,observation[0]));
    psi->set(0, i, 0);
  }

  for (int64_t t=1; t<T; ++t) {
    for (j=0; j<N_; ++j) {
      maxval = 0.0;
      maxvalind = 0;
      for (i=0; i<N_; ++i) {
        val = delta->at(t-1,i)*A_->at(i,j);
        if (val>maxval) {
          maxval = val;
          maxvalind = i;
        }
      }
      delta->set(t, j, maxval*B_->at(j, observation[t]));
      psi->set(t, j, maxvalind);
    }
  }

  real probs = 0.0;
  states.clear();
  states.resize(T);
 // states = new std::vector(T);
  states[T-1] = 1;
  for (i=0; i<N_; ++i) {
    if (delta->at(T-1, i) > probs) {
      probs = delta->at(T-1, i);
      states[T-1] = i;
    }
  }

  for (int64_t t=T-2; t>=0; --t) {
    states[t] = psi->at(t+1, states[t+1]);
  }

  std::cerr << "viterbi解码成功"<<"\n";
  return probs;
}

std::vector<matrixPtr> model::allocXi(int64_t T, int64_t N) {
  std::vector<matrixPtr> xi(T);
  for (int64_t i=0; i<T; ++i) {
    matrixPtr matrixPtr1(new matrix(N,N));
    xi.push_back(matrixPtr1);
  }
  return xi;
}


