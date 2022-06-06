#ifndef OP_NEW
#define OP_NEW
double **cartesian_to_keplerian(double **);
double **J2_accelerations(double **);
double **state_derivative(double**);
double **RK4(double**, double , double, double);
double **orbit_propagate(double**, double, double, double);
#endif