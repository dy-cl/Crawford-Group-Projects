//molecule.cpp

#include "molecule.hpp" // Use molecule class
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>

using namespace std;

// Function to print geometry (cartesian coordinates)
void Molecule::print_geometry(){
    for (int i = 0; i < num_atoms; i++) {

        // Specify output width and precision
        printf("%8.5f %8.5f %8.5f %8.5f\n", atom[i], geometry[i][0], geometry[i][1], geometry[i][2]);
    }     
    cout << "\n"; 
}

// Function to print Hessian matrix as a 9x9 matrix
void Molecule::print_hessian() {
    cout << "Hessian Matrix:" << endl;
        for (int i = 0; i < 3 * num_atoms; i++) {
            for (int j = 0; j < 3 * num_atoms; j++) {
                printf("%12.6f ", hessian[i][j]);
            }
            cout << "\n";
        }
}

// Constructor
Molecule::Molecule(const char *filename) { 

    // Open file and check for success
    ifstream input(filename); 
    if (!input) {
        cerr << "Error opening file: " << filename << endl;
        // Perform error handling or return an error code if necessary
        return;
    }

    string fileExtension = filename;

    if (fileExtension.find("geom.txt") != string::npos) {
        // Read number of atoms
        input >> num_atoms;

        // Allocate memory for arrays
        atom = new double[num_atoms];
        cout << "Atom memory allocated" << "\n";

        geometry = new double*[num_atoms];
        cout << "Geometry memory allocated" << "\n";
        
        for (int i = 0; i < num_atoms; i++) {
            geometry[i] = new double[3];
        }

        // Read geometry from file
        for (int i = 0; i < num_atoms; i++) {
            input >> atom[i] >> geometry[i][0] >> geometry[i][1] >> geometry[i][2];
        }

    } else if (fileExtension.find("hessian.txt") != string::npos) {

        input >> num_atoms;

        // Read Hessian matrix from file
        hessian = new double*[3 * num_atoms];
        cout << "Hessian memory allocated" << "\n";
        cout << "\n";

        for (int i = 0; i < 3 * num_atoms; i++) {
            hessian[i] = new double[3 * num_atoms];
        }

        for (int i = 0; i < 3 * num_atoms; i++) {
            for (int j = 0; j < 3 * num_atoms; j++) {
                input >> hessian[i][j];
            }
        }  

    } else {
        cerr << "Invalid file extension: " << fileExtension << endl;
        // Perform error handling or return an error code if necessary
        return;
    }
    
    input.close();

}

Molecule::~Molecule() {
    cout << "Beginning deallocation" << "\n";

    if (atom != nullptr) {
        cout << "Beginning deallocation of atom" << "\n";
        delete[] atom;
        cout << "Atom deleted" << "\n";
    } else {
        cout << "Atom pointer is null" << "\n";
    }

    if (geometry != nullptr) {
        cout << "Beginning deallocation of geometry" << "\n";
        for (int i = 0; i < num_atoms; i++) {
            delete[] geometry[i];
            cout << "Geometry " << i << "deleted" << "\n";
        }
        delete[] geometry;
        cout << "Geometry deleted" << "\n";
    }

    if (hessian != nullptr) {
        cout << "Beginning deallocation of hessian" << "\n";
        for (int i = 0; i < 3 * num_atoms; i++) {
            delete[] hessian[i];
            cout << "Hessian " << i << "deleted" << "\n"; 
        }
        delete[] hessian;
        cout << "Hessian deleted" << "\n";
    }
}
