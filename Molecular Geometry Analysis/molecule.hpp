#include <string>

using namespace std;

// Molecule class
class Molecule
{
    public:
        // Initialise variables
        int num_atoms;
        int *atom;
        double **geometry;

        // Define functions
        void print_geometry();
        double bond(int atom1, int atom2);
        double unit_vector(int col, int atom1, int atom2);
        double bond_angle(int atom1, int atom2, int atom3);
 
        // Construct and destruct
        Molecule(const char *filename);
        ~Molecule();
};


