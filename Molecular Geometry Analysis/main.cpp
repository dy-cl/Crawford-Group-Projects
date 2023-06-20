#include "molecule.hpp" // Use molecule class
#include <iostream>
#include <cmath>
#include <tuple>
 
using namespace std;
 
int main(int argc, char *argv[]){

    // Create instance of molecule class with data file
    Molecule mol("acetaldehyde.dat");
    
    cout << "Num Atoms: " << mol.num_atoms << "\n";

    // Print geometry 
    mol.print_geometry();

    cout << "Interatomic distances: \n";

    // Print all possible bond distances
    for (int i = 0; i < mol.num_atoms; i++) {
        for (int j = 0; j < i; j++) {
            double distance = mol.bond(i, j);

            // Print the distance
            printf("%d %d %8.5f\n", i, j, distance);
        }
    }

    cout << "Bond angles: \n";

    // Print bond angles
    for (int i = 0; i < mol.num_atoms; i++) {
        for (int j = 0; j < i; j++) {
            for (int k = 0; k < j; k++) {
                if (mol.bond(i, j) < 4 && mol.bond(j, k) < 4) {
                    printf("%2d-%2d-%2d %10.6f\n", i, j, k, mol.bond_angle(i,j,k));
                }
            }
        }
    }
    cout << "Out of plane angles: \n";

    // Print Out of Plane angles
    for(int i = 0; i < mol.num_atoms; i++) {
        for(int k = 0; k < mol.num_atoms; k++) {
            for(int j = 0; j < mol.num_atoms; j++) {
                for(int l = 0; l < j; l++) {
                    if(i!=j && i!=k && i!=l && j!=k && k!=l && mol.bond(i,k) < 4.0 && mol.bond(k,j) < 4.0 && mol.bond(k,l) < 4.0)
                        printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.ooPlane_angle(i,j,k,l));
                }
            }
        }
    }
    cout << "Torsional angles: \n";

    // Print torsional angles
    for(int i = 0; i < mol.num_atoms; i++) {
        for(int j = 0; j < i; j++) {
            for(int k = 0; k < j; k++) {
                for(int l = 0; l < k; l++) {
                    if(i!=j && i!=k && i!=l && j!=k && k!=l && mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0 && mol.bond(k,l) < 4.0)
                        printf("%2d-%2d-%2d-%2d %10.6f\n", i, j, k, l, mol.torsion_angle(i,j,k,l));
                }
            }
        }
    }
    
    // Print centre of mass
    double total_mass = 0.0;

    double Xcm = 0.0;
    double Ycm = 0.0;
    double Zcm = 0.0;

    for (int i = 0; i < mol.num_atoms; i++){
        tuple<double, double, double, double> result = mol.cofMass(mol.atom[i], i);

        double mi, Xi, Yi, Zi;
        tie(mi, Xi, Yi, Zi) = result;

        total_mass += mi;

        Xcm += (mi * Xi);
        Ycm += (mi * Yi);
        Zcm += (mi * Zi);
    }

    Xcm = Xcm / total_mass;
    Ycm = Ycm / total_mass;
    Zcm = Zcm / total_mass;

    cout << "Centre of mass coordinates: \n";

    printf("%10.6f %10.6f %10.6f\n", Xcm, Ycm, Zcm);

    return 0;
} 