//molecule.hpp

#include <string>

using namespace std;

// Molecule class
class Molecule
{
    public:
        // Initialise variables
        int num_atoms;
        double *atom;
        double **geometry;

        // Define functions
        void print_geometry();

        // Construct and destruct
        Molecule(const char *filename);
        ~Molecule();
};