#!/usr/bin/env python3

from astropy.wcs import WCS
from astropy.io import fits
import numpy as np
import argparse as ap

parser = ap.ArgumentParser(description='Fix WCS of ERIS sky images', formatter_class=ap.ArgumentDefaultsHelpFormatter)
parser.add_argument("in_f", type=str, help="Input file name.")
parser.add_argument("--path", type=str, default="./", help="Path to file.")
parser.add_argument("--inplace", action='store_true', help="Update in place.")
parser.add_argument("--out", type=str, default=None, help="Full path to output file. If --inplace, this ignored. If empty, [in_f] in current directory. Note: will overwrite the existing file.")

args = parser.parse_args()
f_name = args.in_f
path = args.path
out = args.out
inplace = args.inplace

SETUPS = {"SPIFFIER" : {
            "250mas" : {"PA" : 65.34, "crpix1" : 41.68, "crpix2" : 19.06, "sign" : 1, "pxscl" : 125, "raoff" : -1.3569e-4, "decoff" : 6.22951e-05}
            },
         "NIX" : {
            "13mas-JHK" : {"PA" : -63.17, "crpix1" : 1064.9, "crpix2" : 1099.1, "sign" : 1, "pxscl" : 13, "raoff" : -2.73765e-4, "decoff" : 1.812239e-4},
            "27mas-JHK" : {"PA" : -65.02, "crpix1" : 958.0, "crpix2": 1065.9, "sign" : -1, "pxscl" : 26.9, "raoff" : -2.73765e-4, "decoff" : 1.812239e-4}
            }
         }

def get_config(path, f_name):

    header = fits.getheader("%s/%s" % (path, f_name))
    INS = header["ESO INS1 SCSM NAME"]
    if INS == "NIX":
        config = SETUPS[INS][header["ESO INS2 NXCW NAME"]]
    elif INS == "SPIFFIER":
        config = SETUPS[INS][header["ESO INS3 SPXW NAME"]]
    else:
        raise ValueError("%s not a valid instrument" % INS)

    return config

def update_wcs(path, f_name, out=None):

    config = get_config(path, f_name)

    head = fits.getheader("%s/%s" % (path, f_name))
    data = fits.getdata("%s/%s" % (path, f_name))

    posang = (config["PA"]+head["ESO ADA POSANG"])*np.pi/180.

    head.set("CTYPE1", "RA--TAN")
    head.set("CTYPE2", "DEC-TAN")
    head.set("CRPIX1", config["crpix1"])
    head.set("CRPIX2", config["crpix2"])
    head.set("CRVAL1", head["RA"]+config["raoff"])
    head.set("CRVAL2", head["DEC"]+config["decoff"])

    rotM = lambda ang : np.array([[np.cos(ang), -np.sin(ang)],
        [np.sin(ang), np.cos(ang)]])

    cdu = np.array([[config["sign"], 0], [0, 1]])
    cd = np.matmul(rotM(posang), cdu)*config["pxscl"]*1e-3/60/60.

    head.set("CD1_1", cd[0][0])
    head.set("CD1_2", cd[0][1])
    head.set("CD2_1", cd[1][0])
    head.set("CD2_2", cd[1][1])

    if inplace:
        fits.PrimaryHDU(data, header=head).writeto("%s/%s" % (path, f_name), overwrite=True)
    else:
        if out is None:
            fits.PrimaryHDU(data, header=head).writeto("./%s" % f_name, overwrite=True)
        else:
            fits.PrimaryHDU(data, header=head).writeto(out, overwrite=True)

if __name__=="__main__":

    update_wcs(path, f_name, out=out)

