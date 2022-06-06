#include <bits/stdc++.h>
#include "programs_new.cpp"
#include "OP_new.cpp"
#include "frame_conversion_new.cpp"

int main(){

    double **state_onboard,**keplerian_elements,**orientation,**w_bi,**moi_satellite,**q_ri_init, **q_br_init,**q_bi_1,**prev_required_var, **trans;

    state_onboard = getzeromatrix(6,1);
    orientation = getzeromatrix(3,1);
    w_bi = getzeromatrix(3,1);
    moi_satellite = getzeromatrix(3,3);
    prev_required_var = getzeromatrix(12,1);

    // state_onboard[0][0] =  -5236.84633;
    // state_onboard[1][0] =  4124.17773;.
    // state_onboard[2][0] =  -1262.94137;

    // state_onboard[3][0] =  -3.86204515;
    // state_onboard[4][0] =  -3.12048032;
    // state_onboard[5][0] =  5.83839029;

    state_onboard[0][0] = 1238.488583514314;
    state_onboard[1][0] = 7080.353211119648;
    state_onboard[2][0] = 14.586295240982;

    state_onboard[3][0] = 1.111578762356251;
    state_onboard[4][0] = -0.204702311583077;
    state_onboard[5][0] = 7.364195645361482;

    keplerian_elements = cartesian_to_keplerian(state_onboard);

    orientation[0][0] = -1.5842;
    orientation[1][0] = 1.4185;
    orientation[2][0] = 0.1864;

    trans = eci_orbit(keplerian_elements);

    q_ri_init  = quatinv(trans_to_quaternion(trans));
    // 0.573816 -0.495255 0.505705 -0.411971

    cout<<q_ri_init[0][0]<<" "<<q_ri_init[1][0]<<" "<<q_ri_init[2][0]<<" "<<q_ri_init[3][0]<<endl;  
    
    return 0;
}