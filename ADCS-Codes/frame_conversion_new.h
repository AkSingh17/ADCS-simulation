#ifndef FRAME_CONVERSION_H
#define FRAME_CONVERSION_H
double **eci_to_ecef(double**,double);
double **ecef_to_lla(double **);
double **ecef_to_eci(double **,double,int ,int );
double **ned_to_ecef_to_eci(double **,double **,double);
double **eci_orbit(double **);
double **eci_body(double **, double**);
#endif // FRAME_CONVERSION_H