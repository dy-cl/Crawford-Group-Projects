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
        printf("%d %8.5f %8.5f %8.5f\n", atom[i], geometry[i][0], geometry[i][1], geometry[i][2]);
    }      
}

// Function to calculate bond lengths (Assume bonds are present between all molecules)
double Molecule::bond(int i, int j){

    // Calculate distance according to formula
    double distance = sqrt(pow(geometry[i][0] - geometry[j][0], 2) +
                           pow(geometry[i][1] - geometry[j][1], 2) +
                           pow(geometry[i][2] - geometry[j][2], 2));

    return distance;
}

// Calculate unit vector 
double Molecule::unit_vector(int col, int p, int q) {
    // 0 = x, 1 = y, 2 = z
    return -(geometry[p][col] - geometry[q][col]) / bond(p,q);
}

// Calculate cross product
tuple<double, double, double> Molecule::cross_product(int a, int b, int c) {
    double cross_x = (unit_vector(1, b, a) * unit_vector(2, b, c)) - (unit_vector(2, b, a) * unit_vector(1, b, c));
    double cross_y = (unit_vector(2, b, a) * unit_vector(0, b, c)) - (unit_vector(0, b, a) * unit_vector(2, b, c));
    double cross_z = (unit_vector(0, b, a) * unit_vector(1, b, c)) - (unit_vector(1, b, a) * unit_vector(0, b, c));

    // cout << cross_x << " " << cross_y << " " << cross_z << "\n";

    // Return tuple
    return make_tuple(cross_x, cross_y, cross_z);
}

// Calculate bond angle in degrees
double Molecule::bond_angle(int i, int j, int k) {

    double dot_product = (unit_vector(0, i, j) * unit_vector(0, k, j)) + // eXij * eXjk
                         (unit_vector(1, i, j) * unit_vector(1, k, j)) + // eYij * eYjk
                         (unit_vector(2, i, j) * unit_vector(2, k, j));  // eZij * eZjk
    
    return (180 / acos(-1.0)) * acos(dot_product); 
}

// Calculate out of plane angle in degrees
double Molecule::ooPlane_angle(int i, int j, int k, int l) {
    
    double eXjkl, eYjkl, eZjkl;

    // Get X, Y, Z components of cross product
    tie(eXjkl, eYjkl, eZjkl) = cross_product(j, k, l);

    // cout << eXjkl << " " << eYjkl << " " << eZjkl << "\n"; // Values are correctly assigned

    // Calculate dot product with eki
    double eX = eXjkl * unit_vector(0, k, i);
    double eY = eYjkl * unit_vector(1, k, i);
    double eZ = eZjkl * unit_vector(2, k, i);
    
    // Calculate sin of angle
    double sin_theta = (eX + eY + eZ) / sqrt(eXjkl * eXjkl + eYjkl * eYjkl + eZjkl * eZjkl);

    return (180 / acos(-1.0)) * asin(sin_theta);
}

// Calculate torsional angle in degrees
double Molecule::torsion_angle(int i, int j, int k, int l) {

    double eXijk, eYijk, eZijk;
    double eXjkl, eYjkl, eZjkl;

    // Get X, Y, Z components of both cross products
    tie(eXijk, eYijk, eZijk) = cross_product(i, j, k);
    tie(eXjkl, eYjkl, eZjkl) = cross_product(j, k, l);

    // Calculate dot product of cross products
    double dot_product1 = (eXijk * eXjkl) + (eYijk * eYjkl) + (eZijk * eZjkl);

    // Denominator
    double dot_product2 = sqrt((eXijk * eXijk) + (eYijk * eYijk) + (eZijk * eZijk)) * sqrt((eXjkl * eXjkl) + (eYjkl * eYjkl) + (eZjkl * eZjkl));

    // Calculate cos_tau in radians
    double cos_tau = dot_product1 / dot_product2;

    // cout << "Cos tau: " << cos_tau << "\n";

    if (cos_tau < -1.0) {
        cos_tau = -1.0;
    } else if (cos_tau > 1.0) {
        cos_tau = 1.0;
    } else {
        cos_tau = cos_tau;
    }

    return (180 / acos(-1.0)) * acos(cos_tau);
}

tuple<double, double, double, double> Molecule::cofMass(int atom, int i) {
    double masses[] = {
        0.000000, // Placeholder 0
        1.007825, // Hydrogen 1
        0.000000, // Placeholder 2
        0.000000, // Placeholder 3
        0.000000, // Placeholder 4
        0.000000, // Placeholder 5
        12.0096,  // Carbon 6
        0.000000, // Placeholder 7
        15.99903, // Oxygen 8
    };

    double mi = masses[atom];

    double Xi = geometry[i][0];
    double Yi = geometry[i][1];
    double Zi = geometry[i][2];

    return make_tuple(mi, Xi, Yi, Zi);
}

// Constructor
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

// Destructor
Molecule::~Molecule(){

    // Deallocate memory
    delete[] atom;
    
    for (int i = 0; i < num_atoms; i++) {
        delete[] geometry[i];
    }

    delete[] geometry;
 }