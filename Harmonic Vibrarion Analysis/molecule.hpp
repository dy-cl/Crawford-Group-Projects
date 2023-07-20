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
        double ** hessian;

        // Define functions
        void print_geometry();
        void print_hessian();
        void weight_hessian();

        // Construct and destruct
        Molecule(const char* geom_filename, const char* hess_filename);
        ~Molecule();
};