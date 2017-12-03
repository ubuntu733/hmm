#include <iostream>
#include <fstream>

#include "model.h"
#include "vector.h"
#include "matrix.h"
int main() {
  vectorPtr initStates(new vector(2));
  matrixPtr transitionProbability(new matrix(2,2));
  matrixPtr emitProbability(new matrix(2,3));
 // model hmmmodel(3,2);
  initStates->operator[](0) = 0.6;
  initStates->operator[](1) = 0.4;
  transitionProbability->set(0, 0 , 0.7);
  transitionProbability->set(0, 1 , 0.3);
  transitionProbability->set(1, 0 , 0.6);
  transitionProbability->set(1, 1 , 0.4);
  emitProbability->set(0, 0, 0.5);
  emitProbability->set(0, 1, 0.4);
  emitProbability->set(0, 2, 0.1);
  emitProbability->set(1, 0, 0.1);
  emitProbability->set(1, 1, 0.3);
  emitProbability->set(1, 2, 0.6);
  model hmmmodel(3, 2, transitionProbability, emitProbability, initStates);
  std::vector<int64_t> observation{0, 1, 2};
  std::cout << hmmmodel.Forward(observation) << "\n";
  std::vector<int64_t> states;
  std::cout << hmmmodel.viterbi(observation, states) << "\n";
  for (auto iter : states) {
    std::cout << iter << "\n";
  }

}