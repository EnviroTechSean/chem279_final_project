#ifndef AO_H
#define AO_H

#include <armadillo>
#include <iostream>

double Overlap_onedim(double xa, double xb, double alphaa, double alphab,
                      int la, int lb);

class Atomic_orbital {
private:
  arma::vec center;
  arma::uvec lmn;
  arma::vec primitive_exponents;
  arma::vec primitive_dcoefs;
  int len;
  std::string atomic_orbital_label;

public:
  Atomic_orbital(arma::vec &R0_input, arma::vec &alpha_input, arma::vec &d_input,
     arma::uvec &lmn_input, std::string lable_input);
  ~Atomic_orbital() {}
  void printinfo();
  arma::uvec get_lmn() { return lmn; }
  arma::vec get_alpha() { return primitive_exponents; }
  arma::vec get_d_coe() { return primitive_dcoefs; }
  arma::vec get_R0() { return center; }
  int get_len() { return len; }
  std::string get_lable() { return atomic_orbital_label; }
};
int GenerateAOs(std::vector<Atomic_orbital> &AOs, const std::string &fname, const arma::mat &H_basis,
                const arma::mat &C_basis);
void PrintAOs(std::vector<Atomic_orbital> &AOs);
double Overlap_3d(arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab,
                  arma::uvec &lmna, arma::uvec &lmnb);

double Eval_Ov_AOs(Atomic_orbital &sh1, Atomic_orbital &sh2);

void Eval_OV_mat(std::vector<Atomic_orbital> &MoleculeAOs, arma::mat &OV_mat);

#endif // AO_H
