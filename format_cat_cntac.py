from astLib.astCoords import decimal2dms, decimal2hms
from astropy.table import Table
from icecream import ic
import numpy as np


def main():
    act = Table.read('act-dr6-zrng.cat', format='ascii.basic')
    act = (act)
    act['block'] = ['A'] * act['name'].size
    act['name'] = [name.replace(' ', '_') for name in act['name']]
    act['hms'] = [decimal2hms(ra, ':') for ra in act['RADeg']]
    act['dms'] = [decimal2dms(dec, ':') for dec in act['decDeg']]
    act['mag'] = ['NA'] * act['name'].size
    act['name','hms','dms','mag'].write('act-dr6-zrng.cntac', format='ascii.csv', overwrite=True)
    return


if __name__ == "__main__":
    main()