//molecule.hpp

#include <string>
#include "eigen-3.4.0/Eigen/Dense"

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
        double ** mass_weighted_hessian;
        Eigen::VectorXd frequencies;

        // Define functions
        void print_geometry();
        void print_hessian();
        void weight_hessian();
        void diagonalize_hessian();

        // Construct and destruct
        Molecule(const char* geom_filename, const char* hess_filename);
        ~Molecule();
};