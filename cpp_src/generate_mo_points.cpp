#include "AO.h"
#include "CNDO.h"
#include <armadillo>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;


int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("usage generate_mo_points filename, for example ./generate_mo_points ../sample_input/H2.txt\n");
    return EXIT_FAILURE;
  }
  string fname(argv[1]);
  try {
    Molecule_basis mol(fname);
    mol.PrintAtoms();

  } catch (invalid_argument &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
