#include "AO.h"
#include "CNDO.h"
#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

// Function to evaluate a single atomic orbital on a 3D grid
arma::cube evaluate_grid_for_orbital(const Atomic_orbital& orbital, const arma::vec& lower_bounds, const arma::vec& upper_bounds, const arma::ivec& grid_points) {
    // Calculate the step size for each dimension
    arma::vec step_size = (upper_bounds - lower_bounds) / (grid_points - 1);

    // Create a cube to store the results
    arma::cube result(grid_points(0), grid_points(1), grid_points(2), arma::fill::zeros);

    // Iterate over each point in the grid
    for (int i = 0; i < grid_points(0); ++i) {
        for (int j = 0; j < grid_points(1); ++j) {
            for (int k = 0; k < grid_points(2); ++k) {
                // Compute the coordinates of the point
                arma::vec point = {
                    lower_bounds(0) + i * step_size(0),
                    lower_bounds(1) + j * step_size(1),
                    lower_bounds(2) + k * step_size(2)
                };
                // Evaluate the orbital at this point
                result(i, j, k) = orbital.evaluate(point);
            }
        }
    }
    return result;
}



int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("usage generate_mo_points filename, for example ./generate_mo_points ../sample_input/H2.txt\n");
    return EXIT_FAILURE;
  }
  string fname(argv[1]);
  try {
    Molecule_basis mol(fname);
    // mol.PrintAtoms();

    // Define the grid properties
    arma::vec lower_bounds = {-5.0, -5.0, -5.0}; // Example bounds
    arma::vec upper_bounds = {5.0, 5.0, 5.0};
    arma::ivec grid_points = {50, 50, 50}; // 50x50x50 grid

    // Evaluate each atomic orbital on the grid
    for (Atom& atom : mol.mAtoms) {
        for (Atomic_orbital& orbital : atom.mAOs) {
            arma::cube orbital_values = evaluate_grid_for_orbital(orbital, lower_bounds, upper_bounds, grid_points);
            std::string output_filename = orbital.get_label() + ".txt";
            orbital_values.save(output_filename, arma::raw_ascii);
        }
    }

  } catch (invalid_argument &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
