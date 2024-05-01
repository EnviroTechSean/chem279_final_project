#include "AO.h"
#include "CNDO.h"
#include <armadillo>
#include <cassert>
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

std::vector<arma::cube> evaluate_molecular_orbitals(Molecule_basis& mol, const arma::vec& lower_bounds, const arma::vec& upper_bounds, const arma::ivec& grid_points, const bool& use_alpha = true)
{
  try
    {
        // Molecule_basis mol(fname);
        CNDO ourSCF(mol, 50, 1e-5); // Need the MO coefficients according to CNDO/2 and SCF
        int ok = ourSCF.init();
        if(ok != 0) {
            throw std::runtime_error("SCF init failed");
        }
        ok = ourSCF.run();
        if(ok != 0) {
            throw std::runtime_error("SCF run failed");
        }
        double Energy = ourSCF.getEnergy();

        std::vector<arma::cube> atomic_orbital_results;
        for (int i = 0; i < mol.mAtoms.size(); i++)
        {
          Atom& A_Atom = mol.mAtoms[i];
          for (int j = 0; j < A_Atom.mAOs.size(); j++)
          {
            Atomic_orbital& A_mu = A_Atom.mAOs[j];
            arma::cube AO_point_evaluations = evaluate_grid_for_orbital(A_mu, lower_bounds, upper_bounds, grid_points);
            atomic_orbital_results.push_back(AO_point_evaluations);
          }
        }

        std::vector<arma::cube> molecular_orbital_results;
        arma::mat C;
        if (use_alpha)
        {
          C = ourSCF.getCa(); // Coefficient matrix for molecular orbitals
        }
        else
        {
          C = ourSCF.getCb();
        }
        assert(C.n_cols == atomic_orbital_results.size());

        for (int p = 0; p < C.n_cols; p++)
        {
            arma::cube molecular_orbital_combination(atomic_orbital_results[0].n_rows, atomic_orbital_results[0].n_cols, atomic_orbital_results[0].n_slices, arma::fill::zeros);
            for (size_t q = 0; q < atomic_orbital_results.size(); q++)
            {
                molecular_orbital_combination += C(q, p) * atomic_orbital_results[q];
            }
            molecular_orbital_results.push_back(molecular_orbital_combination);
        }

        return molecular_orbital_results;
    }
    catch (const invalid_argument &e)
    {
        cerr << "Invalid argument: " << e.what() << endl;
        throw; // rethrow the caught exception
    }
    catch (const std::exception& e) // Catch other standard exceptions
    {
        cerr << "Error: " << e.what() << endl;
        throw; // rethrow the caught exception
    }
}

arma::cube sum_electron_density_cube(const std::vector<arma::cube>& alpha_MOs, const std::vector<arma::cube>& beta_MOs) {
    arma::cube result(alpha_MOs[0].n_rows, alpha_MOs[0].n_cols, alpha_MOs[0].n_slices, arma::fill::zeros);

    // Sum the squares of all the cubes in alpha_MOs
    for (const auto& cube : alpha_MOs) {
        result += arma::square(cube);
    }

    // Sum the squares of all the cubes in beta_MOs
    for (const auto& cube : beta_MOs) {
        result += arma::square(cube);
    }

    // Multiply the accumulated sum by 2 at each point in the cube
    result *= 2;

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
    arma::vec lower_bounds = {-10.0, -10.0, -10.0}; // Example bounds
    arma::vec upper_bounds = {10.0, 10.0, 10.0};
    arma::ivec grid_points = {50, 50, 50}; // 50x50x50 grid
    std::string atom_string = "";

    // Evaluate each atomic orbital on the grid
    for (Atom& atom : mol.mAtoms) {
      atom_string += atom.name;
        for (Atomic_orbital& orbital : atom.mAOs) {
            arma::cube orbital_values = evaluate_grid_for_orbital(orbital, lower_bounds, upper_bounds, grid_points);
            std::string output_filename = orbital.get_label() + ".txt";
            orbital_values.save(output_filename, arma::raw_ascii);
        }
    }

    // Now make the molecular orbitals (linear combinations of atomic orbitals per MO coeffs)
    std::vector<arma::cube> alpha_MOs = evaluate_molecular_orbitals(mol, lower_bounds, upper_bounds, grid_points);
    for (int i = 0; i < alpha_MOs.size(); i++)
    {
      arma::cube current_MO = alpha_MOs[i];
      std::string output_filename = atom_string + "_alpha_MO_" + std::to_string(i) + ".txt";
      current_MO.save(output_filename, arma::raw_ascii);
    }

    // Now make the molecular orbitals (linear combinations of atomic orbitals per MO coeffs)
    std::vector<arma::cube> beta_MOs = evaluate_molecular_orbitals(mol, lower_bounds, upper_bounds, grid_points, false);
    for (int i = 0; i < beta_MOs.size(); i++)
    {
      arma::cube current_MO = beta_MOs[i];
      std::string output_filename = atom_string + "_beta_MO_" + std::to_string(i) + ".txt";
      current_MO.save(output_filename, arma::raw_ascii);
    }

    // Dang professor was right in office hours, once I checked my notes this step was indeed easier than I remembered
    arma::cube electron_densities = sum_electron_density_cube(alpha_MOs, beta_MOs);
    std::string output_filename = atom_string + "_electron_densities" + ".txt";
    electron_densities.save(output_filename, arma::raw_ascii);

  } catch (invalid_argument &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
