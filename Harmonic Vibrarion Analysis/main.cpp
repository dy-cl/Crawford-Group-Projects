// main.cpp

#include "molecule.hpp" // Use molecule class
#include <iostream>
#include <cmath>
#include <tuple>
 
using namespace std;
 
int main(int argc, char *argv[]){

    // Create instance of molecule class with data file
    Molecule geom_mol("h2o_geom.txt");
    Molecule hessian_mol("h2o_hessian.txt");

    cout << "Molecular Geometry:" << "\n";
    cout << "Num Atoms: " << geom_mol.num_atoms << "\n";

    geom_mol.print_geometry();
    hessian_mol.print_hessian();

    return 0;
}