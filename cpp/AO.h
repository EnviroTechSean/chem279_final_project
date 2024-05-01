#ifndef AO_H
#define AO_H

#include <armadillo>
#include <iostream>
#include <vector>

class Atomic_orbital {
private:
  arma::vec center;
  arma::vec exponents;
  arma::vec d_coefs;
  arma::uvec lmn;
  std::string label;
  int num_primatives;
  // Atom *belong;
public:
  Atomic_orbital(arma::vec &center_input, arma::vec &alpha_input, arma::vec &d_input,
     arma::uvec &lmn_input, std::string label_input);
  ~Atomic_orbital() {}
  void printinfo();
  arma::uvec get_lmn() { return lmn; }
  arma::vec get_alpha() { return exponents; }
  arma::vec get_d_coe() { return d_coefs; }
  arma::vec get_center() { return center; }
  int get_num_primatives() { return num_primatives; }
  std::string get_label() { return label; }
  void set_center(arma::vec &centeri) { center = centeri; }
  // void set_belongAtom(Atom *belongi){ belong = belongi;}

  // Adding this for my final project: will be used to calculate many many grid points
  double evaluate(const arma::vec &point) const;

};

class Atom {
public:
  std::vector<Atomic_orbital> mAOs;
  std::string name;
  int VAN; // Valence atomic number
  Atom() : name("0") {}
  Atom(std::vector<Atomic_orbital> AOs) : mAOs(AOs), name("0") {}
  Atom(std::vector<Atomic_orbital> AOs, std::string atomname, int VAN_i)
      : mAOs(AOs), name(atomname), VAN(VAN_i) {}
  Atom(std::string atomname, int VAN_i) : name(atomname), VAN(VAN_i) {}
  void addAO(Atomic_orbital aAO) { mAOs.push_back(aAO); }
  void PrintAOs();
  void set_center(arma::vec &centeri);
};

class Molecule_basis {
public:
  std::vector<Atom> mAtoms;
  int num_ele;
  Molecule_basis() : num_ele(0) {}
  Molecule_basis(std::vector<Atom> Atoms) : mAtoms(Atoms), num_ele(0) {}
  Molecule_basis(std::vector<Atom> Atoms, int cha)
      : mAtoms(Atoms), num_ele(cha) {}
  Molecule_basis(const std::string &fname);
  void addAtom(Atom aAtom) { mAtoms.push_back(aAtom); }
  void setnum_ele(int cha) { num_ele = cha; }
  void PrintAtoms();
  int getnum_ele() { return num_ele; }
  int getnum_atoms() { return mAtoms.size(); }
  int getnumAOs();
  void eval_OVmat(arma::mat &OV_mat);
  void eval_OV1stmat(arma::mat &OV1st_mat);
  // void eval_Hmat(arma::mat &OV_mat, arma::mat &H_mat);
  void eval_gammamat(arma::mat &gamma_mat);
  void eval_gamma1stmat(arma::mat &gamma1st_mat);
  double eval_nuc_E();
  void eval_nuc_1st(arma::mat &nuc_1st);
};

double Eval_Ov_AOs(Atomic_orbital &sh1, Atomic_orbital &sh2);

double Eval_2eI_sAO(Atomic_orbital &sh1, Atomic_orbital &sh2);

#endif // AO_H
