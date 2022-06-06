from astropy.time import Time

from astropy.coordinates import get_sun

t=Time('2014-03-17 23:03:30')

s=get_sun(t)

print(s)