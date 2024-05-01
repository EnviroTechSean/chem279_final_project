#if !defined CNDO_H
#define CNDO_H
#include "AO.h"
#include <armadillo>
#include <cassert>

class SCF {
public:
  virtual int init() = 0;
  virtual int run() = 0;
  virtual arma::mat getPa() = 0;
  virtual arma::mat getPb() = 0;
  virtual double getEnergy() = 0;
};

class CNDO : public SCF {
private:
  Molecule_basis &mol;
  int dim, max_iter;
  int p, q; // num of alpha and beta electrons
  double tol;
  arma::mat S;
  arma::mat gamma;
  arma::mat H_core;
  arma::mat Pa, Pb;
  arma::mat Ga, Gb;
  arma::mat Ca, Cb;
  arma::vec Ea, Eb;
  double Ee, Ec, Etotal;

public:
  CNDO(Molecule_basis &mol_i, int max_it, double tolerence);
  virtual int init();
  int updateG();
  virtual int run();
  virtual arma::mat getPa() { return Pa; }
  virtual arma::mat getPb() { return Pb; }
  virtual arma::mat getCa() { return Ca; }
  virtual arma::mat getCb() { return Cb; }
  virtual double getEnergy();
  arma::mat getGradient();
};

#endif // CNDO_H
