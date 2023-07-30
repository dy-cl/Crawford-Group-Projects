//molecule.cpp

#include "molecule.hpp" // Use molecule class
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Eigenvalues"
#include "eigen-3.4.0/Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
using namespace std;

// Function to print geometry (cartesian coordinates)
void Molecule::print_geometry(){
    for (int i = 0; i < num_atoms; i++) {

        // Specify output width and precision
        printf("%16.12f %16.12f %16.12f %16.12f\n", atom[i], geometry[i][0], geometry[i][1], geometry[i][2]);
    }     
    cout << "\n"; 
}

// Function to print Hessian matrix as a 9x9 matrix
void Molecule::print_hessian() {
    cout << "Hessian Matrix:" << endl;
        for (int i = 0; i < 3 * num_atoms; i++) {
            for (int j = 0; j < 3 * num_atoms; j++) {
                printf("%16.12f ", hessian[i][j]);
            }
            cout << "\n";
        }
        cout << "\n";
}

// Function to mass weight Hessian
void Molecule::weight_hessian() {

    // Atomic masses for oxygen and hydrogen
    const double oxygen_mass = 15.994914619;
    const double hydrogen_mass = 1.0078250322;

    double i_mass;
    double j_mass;
  
    for (int i = 0; i < 3 * num_atoms; i++) {
            for (int j = 0; j < 3 * num_atoms; j++) {

                // Get the corresponding atomic masses for atom i
                if (atom[i / 3] == 8) {
                    i_mass = oxygen_mass;
                } else if (atom[i / 3] == 1) {
                    i_mass = hydrogen_mass;
                }

                // Get the corresponding atomic masses for atom i
                if (atom[j / 3] == 8) {
                    j_mass = oxygen_mass;
                } else if (atom[j / 3] == 1) {
                    j_mass = hydrogen_mass;
                }

                // Calculate weighted hessian value
                mass_weighted_hessian[i][j] = hessian[i][j] / sqrt(i_mass * j_mass);
            }
    }


    cout << "Weighted Hessian Matrix:" << endl;
        for (int i = 0; i < 3 * num_atoms; i++) {
            for (int j = 0; j < 3 * num_atoms; j++) {
                printf("%16.12f ", mass_weighted_hessian[i][j]);
            }
            cout << "\n";
        }
        cout << "\n";
}

// Function to diagonalize hessian
void Molecule::diagonalize_hessian() {

    // Atomic masses for oxygen and hydrogen
    const double oxygen_mass = 15.994914619;
    const double hydrogen_mass = 1.0078250322;

    double i_mass;
    double j_mass;

    // Create a mass-weighted Hessian matrix copy in Eigen format
    Matrix mass_weighted_hessian(3 * num_atoms, 3 * num_atoms);
    for (int i = 0; i < 3 * num_atoms; i++) {
        for (int j = 0; j < 3 * num_atoms; j++) {
            
                // Get the corresponding atomic masses for atom i
                if (atom[i / 3] == 8) {
                    i_mass = oxygen_mass;
                } else if (atom[i / 3] == 1) {
                    i_mass = hydrogen_mass;
                }

                // Get the corresponding atomic masses for atom i
                if (atom[j / 3] == 8) {
                    j_mass = oxygen_mass;
                } else if (atom[j / 3] == 1) {
                    j_mass = hydrogen_mass;
                }

                // Calculate weighted hessian value
                mass_weighted_hessian(i, j) = hessian[i][j] / sqrt(i_mass * j_mass);    
        }
    }

    // Use Eigen library to compute eigenvalues and eigenvectors
    Eigen::SelfAdjointEigenSolver<Matrix> eigensolver(mass_weighted_hessian);

    // Extract eigenvalues and eigenvectors
    Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();

    // Print eigenvalues
    cout << "Eigenvalues:" << "\n";
    for (int i = eigenvalues.size() - 1; i >= 0; i--) {
        cout << i << ": ";
        printf("%16.12f", eigenvalues[i]);
        cout << "\n";
    }

    // Convert eigen values to vibrational frequencies
    double conversion_factor = 5140.484532;
    frequencies = conversion_factor * eigenvalues.cwiseSqrt();

    cout << "\n";

    // Print vibrational frequncies
    cout << "Vibrational Frequencies: " << "\n";
    for (int i = frequencies.size() - 1; i >= 0; i--) {
        if (std::isnan(frequencies(i))) {
            cout << i << ": ";
            printf("%16.12f\n", 0.0); // Print zero for NaN values
        } else {
            cout << i << ": ";
            printf("%16.12f\n", frequencies(i));
        }
    }

}

// Constructor with both geometry and hessian filenames
Molecule::Molecule(const char* geom_filename, const char* hess_filename) {
    // Open geometry file and check for success
    ifstream geom_input(geom_filename);
    if (!geom_input) {
        cerr << "Error opening geometry file: " << geom_filename << endl;
        // Perform error handling or return an error code if necessary
        return;
    }

    // Read number of atoms
    geom_input >> num_atoms;

    // Allocate memory for arrays
    atom = new double[num_atoms];
    geometry = new double*[num_atoms];
    for (int i = 0; i < num_atoms; i++) {
        geometry[i] = new double[3];
    }

    // Read geometry from file
    for (int i = 0; i < num_atoms; i++) {
        geom_input >> atom[i] >> geometry[i][0] >> geometry[i][1] >> geometry[i][2];
    }

    // Close the geometry file
    geom_input.close();

    // Open hessian file and check for success
    ifstream hess_input(hess_filename);
    if (!hess_input) {
        // Perform error handling or return an error code if necessary
        cerr << "Error opening hessian file: " << hess_filename << endl;

        // Clean up the memory allocated for geometry data
        delete[] atom;
        for (int i = 0; i < num_atoms; i++) {
            delete[] geometry[i];
        }
        delete[] geometry;

        return;
    }

    // Skip the first line in the hessian file (header line)
    string header_line;
    getline(hess_input, header_line);

    // Allocate memory for the hessian array
    hessian = new double*[3 * num_atoms];
    for (int i = 0; i < 3 * num_atoms; i++) {
        hessian[i] = new double[3 * num_atoms];
    }

    // Allocate memory for the mass_weighted_hessian array
    mass_weighted_hessian = new double*[3 * num_atoms];
    for (int i = 0; i < 3 * num_atoms; i++) {
        mass_weighted_hessian[i] = new double[3 * num_atoms];
    }

    // Initialize the mass_weighted_hessian to zero
    for (int i = 0; i < 3 * num_atoms; i++) {
        for (int j = 0; j < 3 * num_atoms; j++) {
            mass_weighted_hessian[i][j] = 0.0;
        }
    }

    //Allocate memory for vibrational frequencies
    frequencies.resize(3 * num_atoms);
    
    // Read hessian from file
    for (int i = 0; i < 3 * num_atoms; i++) {
        for (int j = 0; j < 3 * num_atoms; j++) {
            hess_input >> hessian[i][j];
        }
    }

    // Close the hessian file
    hess_input.close();

    
}

// Destructor
Molecule::~Molecule() {

    // Deallocate memory for atom data
    if (atom != nullptr) {
        delete[] atom;
    }

    // Deallocate memory for geometry data
    if (geometry != nullptr) {
        for (int i = 0; i < num_atoms; i++) {
            delete[] geometry[i];
        }
        delete[] geometry;
    }

    // Deallocate memory for hessian data
    if (hessian != nullptr) {
        for (int i = 0; i < 3 * num_atoms; i++) {
            delete[] hessian[i];
        }
        delete[] hessian;
    }

    // Deallocate memory for mass_weighted_hessian data
    if (mass_weighted_hessian != nullptr) {
        for (int i = 0; i < 3 * num_atoms; i++) {
            delete[] mass_weighted_hessian[i];
        }
        delete[] mass_weighted_hessian;
    }
}
