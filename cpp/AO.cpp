#include "AO.h"
#include "util.h"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

using namespace std;
double hartree_to_ev = 27.211396641308;

// AO functions
void Atomic_orbital::printinfo() {
  printf("This AO info: %s, R( %1.4f, %1.4f, %1.4f), with angular momentum: "
         "%lld %lld %lld\n",
         label.c_str(), center(0), center(1), center(2), lmn(0), lmn(1),
         lmn(2));
  d_coefs.print("d_coefs");
  exponents.print("exponents");
}

Atomic_orbital::Atomic_orbital(arma::vec &center_input, arma::vec &exponents_input,
                               arma::vec &d_coeffs_input, arma::uvec &lmn_input,
                               string label_input)
    : center(center_input), exponents(exponents_input), d_coefs(d_coeffs_input),
      lmn(lmn_input), label(label_input) {
  assert(center.n_elem == 3);
  assert(lmn.n_elem == 3);
  assert(d_coefs.n_elem == exponents.n_elem);
  num_primatives = exponents.n_elem;
  for (size_t k = 0; k < num_primatives; k++) {
    double Overlap_Self =
        Overlap_3d(center, center, exponents(k), exponents(k), lmn, lmn);
    d_coefs(k) /= sqrt(Overlap_Self);
  }
}

double Atomic_orbital::evaluate(const arma::vec &point) const
{
  double total_value = 0.0; // The return value, after updates in loop

  arma::vec distance = point - center;
  double distance_from_center_squared = arma::dot(distance, distance);
  double x_diff = distance[0]; 
  double y_diff = distance[1]; 
  double z_diff = distance[2]; 

  double primitive_contribution;

  for (size_t i = 0; i < num_primatives; i++)
  {
    // Start by multiplying the x y and z distances to the primitive contribution
    primitive_contribution = pow(x_diff, lmn[0]) * pow(y_diff, lmn[1]) * pow(z_diff, lmn[2]);
    // Now add normalized d coefficient for the relevant gaussian
    primitive_contribution *= d_coefs[i];
    // Now do the exponential term
    primitive_contribution *= exp(-exponents[i] * distance_from_center_squared);

    total_value += primitive_contribution;
  }
  return total_value;
}


std::map<std::string, int> VAN_map{
    {"H", 1}, {"C", 4}, {"N", 5}, {"O", 6}, {"F", 7},
};
Atom GenerateAtom(std::string atomname, arma::vec center) {
  std::string basisname = "../basis/" + atomname + "_STO3G.txt";
  arma::mat basis;
  basis.load(basisname);

  arma::uvec lmn = {0, 0, 0};
  arma::vec alpha = basis.col(0);
  arma::vec d_coe = basis.col(1);
  string label("s");

  // Convert each double in the center vector to string with precision and add to label
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2); // Set precision to 2 decimal places
  for (size_t i = 0; i < center.size(); ++i) {
    oss << "_" << center(i);
  }
  // Append coordinates to label
  label += oss.str();

  Atomic_orbital AO_s(center, alpha, d_coe, lmn, label);
  if (atomname == string("H")) {
    Atom aatom(atomname, 1);
    aatom.addAO(AO_s);
    return aatom;
  }
  if (VAN_map.find(atomname) == VAN_map.end()) {
    throw invalid_argument("Do not support this kind of atom.");
  }
  int atomicnumber = VAN_map[atomname];
  Atom aatom(atomname, atomicnumber);
  aatom.addAO(AO_s);
  for (size_t j = 0; j < 3; j++) {
    d_coe = basis.col(2);
    lmn.zeros();
    lmn(j) = 1;
    string label("p");
    Atomic_orbital readedAOp(center, alpha, d_coe, lmn, label);
    aatom.addAO(readedAOp);
  }
  return aatom;
}

double Eval_Ov_AOs(Atomic_orbital &ao1, Atomic_orbital &ao2) {
  int len = ao1.get_num_primatives();
  assert(ao2.get_num_primatives() == len);
  arma::vec alphaa = ao1.get_alpha(), alphab = ao2.get_alpha();
  arma::vec Ra = ao1.get_center(), Rb = ao2.get_center();
  arma::uvec la = ao1.get_lmn(), lb = ao2.get_lmn();
  arma::vec da = ao1.get_d_coe(), db = ao2.get_d_coe();

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

void Eval_Ov1st_AOs(arma::vec &OV1st, Atomic_orbital &ao1,
                    Atomic_orbital &ao2) {
  int len = ao1.get_num_primatives();
  assert(ao2.get_num_primatives() == len);
  arma::vec alphaa = ao1.get_alpha(), alphab = ao2.get_alpha();
  arma::vec Ra = ao1.get_center(), Rb = ao2.get_center();
  arma::uvec la = ao1.get_lmn(), lb = ao2.get_lmn();
  arma::vec da = ao1.get_d_coe(), db = ao2.get_d_coe();

  OV1st.zeros(3);
  arma::vec Overlap(3);
  for (size_t k = 0; k < len; k++) {
    double alpha_k = alphaa(k);
    for (size_t j = 0; j < len; j++) {
      double alpha_j = alphab(j);
      Overlap1st_3d(Overlap, Ra, Rb, alpha_k, alpha_j, la, lb);
      // Overlap.print("Overlap");
      OV1st += da(k) * db(j) * Overlap;
    }
  }
  // OV1st.print("OV1st in ao");
}

double Eval_2eI_sAO(Atomic_orbital &ao1, Atomic_orbital &ao2) {
  // ao1.printinfo();
  // ao2.printinfo();
  arma::uvec la = ao1.get_lmn(), lb = ao2.get_lmn();
  if (!(arma::accu(la) == 0 && arma::accu(lb) == 0))
    throw invalid_argument("Eval_2eI_sAO is only used for s Orbitals.");
  int len = ao1.get_num_primatives();
  assert(ao2.get_num_primatives() == len);
  arma::vec Ra = ao1.get_center(), Rb = ao2.get_center();
  arma::vec alphaa = ao1.get_alpha(), alphab = ao2.get_alpha();
  arma::vec da = ao1.get_d_coe(), db = ao2.get_d_coe();

  double sum = 0.;
  for (size_t k1 = 0; k1 < len; k1++)
    for (size_t k2 = 0; k2 < len; k2++) {
      double sigmaA = 1.0 / (alphaa(k1) + alphaa(k2));
      for (size_t j1 = 0; j1 < len; j1++)
        for (size_t j2 = 0; j2 < len; j2++) {
          double sigmaB = 1.0 / (alphab(j1) + alphab(j2));
          double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB);
          // printf("%ld %ld %ld %ld = %1.10f\n", k1, k2, j1, j2, I2e);
          // return I2e;
          sum += da(k1) * da(k2) * db(j1) * db(j2) * I2e;
        }
    }
  return sum;
}

void Eval_2eI1st_sAO(arma::vec &Gamma1st, Atomic_orbital &ao1,
                     Atomic_orbital &ao2) {
  // ao1.printinfo();
  // ao2.printinfo();
  arma::uvec la = ao1.get_lmn(), lb = ao2.get_lmn();
  if (!(arma::accu(la) == 0 && arma::accu(lb) == 0))
    throw invalid_argument("Eval_2eI_sAO is only used for s Orbitals.");
  int len = ao1.get_num_primatives();
  assert(ao2.get_num_primatives() == len);
  arma::vec Ra = ao1.get_center(), Rb = ao2.get_center();
  arma::vec alphaa = ao1.get_alpha(), alphab = ao2.get_alpha();
  arma::vec da = ao1.get_d_coe(), db = ao2.get_d_coe();

  Gamma1st.zeros(3);
  arma::vec gam1st(3);
  for (size_t k1 = 0; k1 < len; k1++)
    for (size_t k2 = 0; k2 < len; k2++) {
      double sigmaA = 1.0 / (alphaa(k1) + alphaa(k2));
      for (size_t j1 = 0; j1 < len; j1++)
        for (size_t j2 = 0; j2 < len; j2++) {
          double sigmaB = 1.0 / (alphab(j1) + alphab(j2));
          I2e_pG_1st(gam1st, Ra, Rb, sigmaA, sigmaB);
          Gamma1st += da(k1) * da(k2) * db(j1) * db(j2) * gam1st;
        }
    }
}

std::map<std::string, std::string> AN_map{
    {"1", "H"}, {"6", "C"}, {"7", "N"}, {"8", "O"}, {"9", "F"}};
Molecule_basis::Molecule_basis(const string &fname) {
  int num_charge, num_Atoms;

  ifstream in(fname, ios::in);
  // cout << fname;

  string line;
  getline(in, line);
  istringstream iss(line);
  if (!(iss >> num_Atoms >> num_charge))
    throw invalid_argument("There is some problem with molecule format.");
  int count_atoms = 0;

  while (getline(in, line)) {
    istringstream iss(line);
    arma::vec center(3);
    // int AtomicN = 0;
    string atomnumber;
    if (!(iss >> atomnumber >> center[0] >> center[1] >> center[2]))
      throw invalid_argument("There is some problem with AO format.");

    if (AN_map.find(atomnumber) == AN_map.end()) {
      throw invalid_argument("Do not support this kind of atom.");
    }

    string atomname = AN_map[atomnumber];
    arma::uvec lmn = {0, 0, 0};
    Atom readAtom = GenerateAtom(atomname, center);
    mAtoms.push_back(readAtom);
    // cout << readAO << std::endl;
    count_atoms++;
  }
  if (count_atoms != num_Atoms) {
    throw invalid_argument("Number of AOs are not consistent ");
  }
  in.close();
  num_ele = 0;
  for (auto atom : mAtoms)
    num_ele += atom.VAN;
  num_ele -= num_charge;
}

void Atom::PrintAOs() {
  printf("This is a %s atom\n", name.c_str());
  for (auto ao : mAOs)
    ao.printinfo();
  printf("\n");
}
void Atom::set_center(arma::vec &centeri) {
  for (auto &ao : mAOs)
    ao.set_center(centeri);
}

void Molecule_basis::PrintAtoms() {
  for (auto atom : mAtoms)
    atom.PrintAOs();
}
int Molecule_basis::getnumAOs() {
  int numAOs = 0;
  for (auto atom : mAtoms)
    numAOs += atom.mAOs.size();
  return numAOs;
}

void Molecule_basis::eval_OVmat(arma::mat &OV_mat) {
  int dim = getnumAOs();
  assert(OV_mat.n_rows == dim && OV_mat.n_cols == dim);
  vector<Atomic_orbital> allAOs;
  for (auto atom : mAtoms)
    for (auto ao : atom.mAOs)
      allAOs.push_back(ao);
  for (size_t k = 0; k < dim; k++) {
    for (size_t j = 0; j <= k; j++) {
      double OV_1AO = Eval_Ov_AOs(allAOs[k], allAOs[j]);
      OV_mat(k, j) = OV_1AO;
      OV_mat(j, k) = OV_1AO;
    }
  }
}

// Fast index is xyz, then atom A, finally B
void Molecule_basis::eval_OV1stmat(arma::mat &OV_1stmat) {
  int dim = getnumAOs();
  // cout<< OV_1stmat.n_rows <<" " <<OV_1stmat.n_rows;
  assert(OV_1stmat.n_rows == 3 && OV_1stmat.n_cols == dim * dim);
  vector<Atomic_orbital> allAOs;
  for (auto atom : mAtoms)
    for (auto ao : atom.mAOs)
      allAOs.push_back(ao);
  for (size_t k = 0; k < dim; k++)
    for (size_t j = 0; j < dim; j++) {
      arma::vec OV_1AO(OV_1stmat.colptr(k * dim + j), 3, false, true);
      if (arma::approx_equal(allAOs[k].get_center(), allAOs[j].get_center(), "absdiff",
                             1e-5))
        OV_1AO.zeros();
      else
        Eval_Ov1st_AOs(OV_1AO, allAOs[k], allAOs[j]);
      // OV_1AO.print("OV_1AO");
    }
  // OV_1stmat.print("OV_1stmat in func");
}

void Molecule_basis::eval_gammamat(arma::mat &gamma_mat) {
  int dim = mAtoms.size();
  assert(gamma_mat.n_rows == dim && gamma_mat.n_rows == dim);
  for (size_t k = 0; k < dim; k++) {
    for (size_t j = 0; j <= k; j++) {
      double OV_1AO =
          Eval_2eI_sAO(mAtoms[k].mAOs[0], mAtoms[j].mAOs[0]);
      gamma_mat(k, j) = OV_1AO;
      gamma_mat(j, k) = OV_1AO;
    }
  }
  gamma_mat *= hartree_to_ev;
}

void Molecule_basis::eval_gamma1stmat(arma::mat &gamma1st_mat) {
  int dim = mAtoms.size();
  assert(gamma1st_mat.n_rows == 3 && gamma1st_mat.n_cols == dim * dim);
  for (size_t k = 0; k < dim; k++) {
    for (size_t j = 0; j < dim; j++) {
      arma::vec Ga_1AO(gamma1st_mat.colptr(k * dim + j), 3, false, true);
      if (k == j)
        Ga_1AO.zeros();
      else
        Eval_2eI1st_sAO(Ga_1AO, mAtoms[k].mAOs[0],
                        mAtoms[j].mAOs[0]);
    }
  }
  gamma1st_mat *= hartree_to_ev;
}

double Molecule_basis::eval_nuc_E() {
  double Ec = 0;
  int dim = mAtoms.size();
  for (size_t k = 0; k < dim; k++) {
    arma::vec Ra = mAtoms[k].mAOs[0].get_center();
    for (size_t j = 0; j < k; j++) {
      arma::vec Rb = mAtoms[j].mAOs[0].get_center();
      double Rd = arma::norm(Ra - Rb, 2);
      Ec += mAtoms[k].VAN * mAtoms[j].VAN / Rd;
    }
  }
  return Ec * hartree_to_ev;
}
void Molecule_basis::eval_nuc_1st(arma::mat &nuc_1st) {
  nuc_1st.zeros();
  int dim = mAtoms.size();
  assert(nuc_1st.n_rows == 3 && nuc_1st.n_cols == dim);
  for (size_t k = 0; k < dim; k++) {
    arma::vec gradient_A(nuc_1st.colptr(k), 3, false, true);
    arma::vec Ra = mAtoms[k].mAOs[0].get_center();
    for (size_t j = 0; j < dim; j++) {
      if (j == k)
        continue;
      arma::vec Rb = mAtoms[j].mAOs[0].get_center();
      double Rd = arma::norm(Ra - Rb, 2);
      gradient_A -= mAtoms[k].VAN * mAtoms[j].VAN / pow(Rd, 3) * (Ra - Rb);
    }
  }
  nuc_1st *= hartree_to_ev;
}
