#include "molecule.h" // Use molecule class
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>

// Function to print geometry (cartesian coordinates)
void Molecule::print_geometry(){
    for (int i = 0; i < num_atoms; i++) {

        // Specify output width and precision
        printf("%d %8.5f %8.5f %8.5f\n", atom[i], geometry[i][0], geometry[i][1], geometry[i][2]);
    }      
}

// 
Molecule::Molecule(const char *filename) { 

    // Open file and check for success
    ifstream input(filename);
    assert(input.good());

    // Read first line (number of atoms)
    input >> num_atoms;

    // Allocate memory for array storage
    atom = new int[num_atoms];
    geometry = new double* [num_atoms];
    for (int i = 0; i < num_atoms; i++) {
        geometry[i] = new double[3];
    }

    // Read geomtry from file
    for (unsigned int i = 0; i < num_atoms; i++) { // Assume num_atoms is non-negative
        input >> atom[i] >> geometry[i][0] >> geometry[i][1] >> geometry[i][2];
    }

    input.close();
}

Molecule::~Molecule(){

    // Deallocate memory
    delete[] atom;
    

    for (int i = 0; i < num_atoms; i++) {
        delete[] geometry[i];
    }

    delete[] geometry;
 }