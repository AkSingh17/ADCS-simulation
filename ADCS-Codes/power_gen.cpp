#include <bits/stdc++.h>
#include "programs_new.h"

using namespace std;

double solargen(double **sun,double eclipse,double intensity,double eol_efficiency){

    double **s = Unit_vector(sun,3,1);
    // double **total_gen= getzeromatrix(3,1);

    double generated_power_x = 0;
    double generated_power_y = 0;
    double generated_power_z = 0;

    //solar panel areas
    double A_x = 0.00562435*2;
    double A_y = 0.00562435*2;
    double A_z = 0.00562435;

    //check to get power only when the exposed side of the panel is getting the sun vector

    if((s[0][0]*A_x < 0))
    generated_power_x = -intensity*eol_efficiency*eclipse*(s[0][0]*A_x);

    if((s[1][0]*A_x < 0))
    generated_power_y = -intensity*eol_efficiency*eclipse*(s[1][0]*A_y);

    if((s[2][0]*A_x < 0))
    generated_power_z = -intensity*eol_efficiency*eclipse*(s[2][0]*A_z);

    double total_gen_power = generated_power_x + generated_power_y + generated_power_z;
    double total_gen_energy = total_gen_power*1;
    double total_gen_soc = total_gen_energy;

    // total_gen[0][0] = total_gen_power;
    // total_gen[1][0] = total_gen_energy;

    free_variable(s,3);

    return total_gen_power;
}