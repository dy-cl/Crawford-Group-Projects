// main.cpp

#include "molecule.hpp" // Use molecule class
#include <iostream>
#include <cmath>
#include <tuple>
 
using namespace std;
 
int main(int argc, char *argv[]){

    // Create instance of molecule class with data file
    Molecule mol("h2o_geom.txt");
    
    cout << "Num Atoms: " << mol.num_atoms << "\n";

    // Print geometry 
    mol.print_geometry();

    return 0;
}