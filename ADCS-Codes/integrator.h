#ifndef INTEGRATOR_H
#define INTEGRATOR_H
double** funcEval_dwdt(double** , double** , double** , double** , double** , double**);
double** funcEval_dqdt(double** , double** );
double **integrate_omega(double** , double** , double** , double** , double**, double**, double , double );
double **integrate_quaternion(double**, double**, double, double);
#endif  //INTEGRATOR_H