// main.cpp

#include "molecule.hpp" // Use molecule class
#include <iostream>
#include <cmath>
#include <tuple>
 
using namespace std;
 
int main(int argc, char *argv[]){

    // Create instance of molecule class with data file
    Molecule mol("h2o_geom.txt", "h2o_hessian.txt");

    cout << "Molecular Geometry:" << "\n";
    cout << "Num Atoms: " << mol.num_atoms << "\n";

    mol.print_geometry();
    mol.print_hessian();
    mol.weight_hessian();
    mol.diagonalize_hessian();

    return 0;
}