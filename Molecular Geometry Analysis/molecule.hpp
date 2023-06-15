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
        tuple<double, double, double> cross_product(int atom1, int atom2, int atom3);
        double bond_angle(int atom1, int atom2, int atom3);
        double ooPlane_angle(int atom1, int atom2, int atom3, int atom4);
        double torsion_angle(int atom1, int atom2, int atom3, int atom4);
 
        // Construct and destruct
        Molecule(const char *filename);
        ~Molecule();
};


