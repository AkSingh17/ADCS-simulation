import numpy as np
from sgp4.api import jday
from sgp4.api import Satrec
from sgp4.api import SGP4_ERRORS 
from typing import Dict

import time


def sgp4(init_time: float, read_interval: int, init_tle: list) -> Dict:
    """
        Simulate the orbital data (position and velocity) and act as a gps reading
        takes in an initial epoch, TLE value, and interval between successive gps readings
        :param init_time: initial epoch in seconds from 01.01.1970
        :param read_interval: time interval between successive gps readings (in int secs)
        :param init_tle: tle at the last/initial satellite position
        :returns : error detection flag, position and velocity in TEME frame
        
        """
    # A tuple giving the time in date hour min sec format which is required for the jday function
    time_next = time.localtime(init_time + read_interval)
    # jday gives the julian date from the above tuples time format
    j_date, precise_fraction = jday(*(time_next[0:6]))
    # the satellites TLEs are given next.
    s_tle, t_tle = init_tle
    # TLEs are converted to cartesian vectors
    satellite = Satrec.twoline2rv(s_tle, t_tle)

    # these are propagated till the required time
    error_message, position, velocity = satellite.sgp4(j_date, precise_fraction)
    #print(position)
    #print(velocity)
    return {"err": error_message,
            "pos": position,
            "vel": velocity}

out_file = open('ADCS-Codes/sgp4_position.txt','w')

for i in range(1,10000):
 
    dicty = sgp4(200000,i,['1 35932U 09051B   22045.84189200  .00000486  00000+0  11983-3 0  9994','2 35932  98.5829 258.2562 0008613  67.1042 293.1063 14.56688797658243'])
    position = str(dicty['pos'])
    vel = str(dicty['vel'])

    out_file.write(position)
    out_file.write('\n')