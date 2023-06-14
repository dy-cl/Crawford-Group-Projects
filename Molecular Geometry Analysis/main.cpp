#include "molecule.hpp" // Use molecule class
#include <iostream>
#include <cmath>
 
using namespace std;
 
int main(int argc, char *argv[]){

    // Create instance of molecule class with data file
    Molecule mol("acetaldehyde.dat");
    
    cout << "Num Atoms: " << mol.num_atoms << "\n";

    // Call print geometry function on the file
    mol.print_geometry();

    cout << "Interatomic distances: \n";

    // Calculate all possible bond distances
    for (int i = 0; i < mol.num_atoms; i++) {
        for (int j = 0; j < i; j++) {
            double distance = mol.bond(i, j);

            // Print the distance
            printf("%d %d %8.5f\n", i, j, distance);
        }
    }

    cout << "Bond angles: \n";

    // Calculate bond angles
    for (int i = 0; i < mol.num_atoms; i++) {
        for (int j = 0; j < i; j++) {
            for (int k = 0; k < j; k++) {
                if (mol.bond(i, j) < 4 && mol.bond(j, k) < 4) {
                    printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.bond_angle(i,j,k));
                }
            }
        }
    }
    return 0;
} 