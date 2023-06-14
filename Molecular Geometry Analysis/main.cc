#include "molecule.h" // Use molecule class
 
using namespace std;
 
int main(int argc, char *argv[]){

    // Create instance of molecule class with data file
    Molecule acetaldehyde("acetaldehyde.dat");
    
    // Call print geometry function on the file
    acetaldehyde.print_geometry();
    
    return 0;
}