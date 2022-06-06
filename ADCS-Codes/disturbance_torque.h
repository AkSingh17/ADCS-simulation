#ifndef DISTURBANCE_TORQUE_H
#define DISTURBANCE_TORQUE_H
double  **drag_torque(double **, double **, double **, double **);
double **gravity_torque(double **, double **,double **);
double **solar_torque(double **, double **);
double **calc_disturbance_torque(double **, double **, double **, double **);
#endif // DISTURBANCE_TORQUE_H