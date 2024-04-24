#if !defined EH_H
#define EH_H
#include "AO.h"
#include <armadillo>
#include <cassert>

void Generate_Hmat(arma::mat &OV_mat, std::vector<Atomic_orbital> &AOs, arma::mat &H_mat);

double Solve_EH(arma::mat &OV_mat, arma::mat &H_mat, arma::mat &C_mat,
                arma::vec &energy_vec, int num_ele);

#endif // EH_H
