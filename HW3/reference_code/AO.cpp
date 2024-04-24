#include "AO.h"
#include "util.h"
#include <cassert>
#include <fstream>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

// return number of electrons
int GenerateAOs(vector<Atomic_orbital> &AOs, const string &fname, const arma::mat &H_basis,
                const arma::mat &C_basis) {
  int basislen = H_basis.n_rows;
  assert(C_basis.n_rows == basislen);
  int num_charge;
  int num_Atoms;

  ifstream in(fname, ios::in);
  // cout << fname;

  string line;
  getline(in, line);
  istringstream iss(line);
  if (!(iss >> num_Atoms >> num_charge))
    throw invalid_argument("There is some problem with AO format.");
  int count_atoms = 0;

  while (getline(in, line)) {
    istringstream iss(line);
    arma::vec R0(3);
    int AtomicN = 0;
    if (!(iss >> AtomicN >> R0[0] >> R0[1] >> R0[2]))
      throw invalid_argument("There is some problem with AO format.");

    arma::uvec lmn = {0, 0, 0};
    arma::vec alpha(basislen);
    arma::vec d_coe(basislen);
    if (AtomicN == 1) {
      alpha = H_basis.col(0);
      d_coe = H_basis.col(1);
      string lable("H1s");
      Atomic_orbital readedAO(R0, alpha, d_coe, lmn, lable);
      AOs.push_back(readedAO);
    } else {
      if (AtomicN == 6) {
        alpha = C_basis.col(0);
        d_coe = C_basis.col(1);
        string lable("C2s");
        Atomic_orbital readedAO(R0, alpha, d_coe, lmn, lable);
        AOs.push_back(readedAO);
        for (size_t j = 0; j < 3; j++) {
          d_coe = C_basis.col(2);
          lmn.zeros();
          lmn(j) = 1;
          string lable("C2p");
          Atomic_orbital readedAOp(R0, alpha, d_coe, lmn, lable);
          AOs.push_back(readedAOp);
        }
      } else
        throw invalid_argument("There are AOs other than H and C.");
    }
    // cout << readAO << std::endl;
    count_atoms++;
  }
  if (count_atoms != num_Atoms) {
    throw invalid_argument("Number of AOs are not consistent ");
  }
  in.close();
  return AOs.size() - num_charge;
}

// void AO::Reset(double x0_input, double y0_input, double z0_input, double
// alpha_input, int l_input){
//   R0(0)=x0_input; R0(1)=y0_input; R0(2)=z0_input; alpha=alpha_input;
//   l=l_input;
// }
void PrintAOs(std::vector<Atomic_orbital> &AOs) {
  for (auto ao : AOs)
    ao.printinfo();
}

void Atomic_orbital::printinfo() {
  printf("This AO info: %s, R( %1.2f, %1.2f, %1.2f), with angular momentum: "
         "%lld %lld %lld\n",
         atomic_orbital_label.c_str(), center(0), center(1), center(2), lmn(0), lmn(1), lmn(2));
  primitive_dcoefs.print("d_coe");
  primitive_exponents.print("alpha");
}

double Overlap_onedim(double xa, double xb, double alphaa, double alphab,
                      int la, int lb) {
  // double x = Combination(3, 1);
  // double y = DoubleFactorial(4);
  double prefactor =
      exp(-alphaa * alphab * (xa - xb) * (xa - xb) / (alphaa + alphab)) *
      sqrt(M_PI / (alphaa + alphab));
  double xP = (alphaa * xa + alphab * xb) / (alphaa + alphab);

  double result = 0.0;
  for (int i_index = 0; i_index <= la; i_index++)
    for (int j_index = 0; j_index <= lb; j_index++) {
      if ((i_index + j_index) % 2 == 1)
        continue;
      double C_part = Combination(la, i_index) * Combination(lb, j_index);
      double DF_part = DoubleFactorial(i_index + j_index - 1);
      double numerator =
          pow(xP - xa, la - i_index) * pow(xP - xb, lb - j_index);
      // Caution: convert i_index + j_index to float!
      double dominator =
          pow(2 * (alphaa + alphab), double(i_index + j_index) / 2.0);
      double temp = C_part * DF_part * numerator / dominator;
      result += temp;
      // printf("%f %f %f %f  %f\n", C_part, DF_part, numerator, dominator,
      // temp);
    }

  result *= prefactor;
  return result;
}

double Overlap_3d(arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab,
                  arma::uvec &lmna, arma::uvec &lmnb) {
  double Overlap =
      Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, lmna(0), lmnb(0)) *
      Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, lmna(1), lmnb(1)) *
      Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, lmna(2), lmnb(2));
  return Overlap;
}

Atomic_orbital::Atomic_orbital(arma::vec &R0_input, arma::vec &alpha_input, arma::vec &d_input,
       arma::uvec &lmn_input, string lable_input)
    : center(R0_input), primitive_exponents(alpha_input), primitive_dcoefs(d_input), lmn(lmn_input),
      atomic_orbital_label(lable_input) {
  assert(center.n_elem == 3);
  assert(lmn.n_elem == 3);
  len = primitive_exponents.n_elem;
  assert(primitive_dcoefs.n_elem == len);
  for (size_t k = 0; k < len; k++) {
    double Overlap_Self = Overlap_3d(center, center, primitive_exponents(k), primitive_exponents(k), lmn, lmn);
    primitive_dcoefs(k) /= sqrt(Overlap_Self);
  }
}

double Eval_Ov_AOs(Atomic_orbital &sh1, Atomic_orbital &sh2) {

  int len = sh1.get_len();
  assert(sh2.get_len() == len);
  arma::vec alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();
  arma::vec Ra = sh1.get_R0(), Rb = sh2.get_R0();
  arma::uvec la = sh1.get_lmn(), lb = sh2.get_lmn();
  arma::vec da = sh1.get_d_coe(), db = sh2.get_d_coe();

  double sum = 0.;
  for (size_t k = 0; k < len; k++) {
    double alpha_k = alphaa(k);
    for (size_t j = 0; j < len; j++) {
      double alpha_j = alphab(j);
      double Overlap = Overlap_3d(Ra, Rb, alpha_k, alpha_j, la, lb);
      // printf("%ld %ld = %1.10f\n", k, j, Overlap);
      sum += da(k) * db(j) * Overlap;
    }
  }
  return sum;
}

void Eval_OV_mat(vector<Atomic_orbital> &MoleculeAOs, arma::mat &OV_mat) {
  int dim = MoleculeAOs.size();
  assert(OV_mat.n_rows == dim && OV_mat.n_cols == dim);
  for (size_t k = 0; k < dim; k++) {
    for (size_t j = 0; j <= k; j++) {
      double OV_1AO = Eval_Ov_AOs(MoleculeAOs[k], MoleculeAOs[j]);
      OV_mat(k, j) = OV_1AO;
      OV_mat(j, k) = OV_1AO;
    }
  }
}
