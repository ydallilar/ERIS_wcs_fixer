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

"""
SETUPS = {"SPIFFIER" : {
            "250mas" : {"PA" : 68.54, "crpix1" : 64-39.86, "crpix2" : 19.08, "sign" : -1, "pxscl" : 125, "raoff" : -2.0392e-4, "decoff" : 1.3612e-4}
            #"250mas" : {"PA" : 68.54, "crpix1" : 39.86, "crpix2" : 19.08, "sign" : 1, "pxscl" : 125, "raoff" : -2.0392e-4, "decoff" : 1.3612e-4}
            },
         "NIX" : {
            "13mas-JHK" : {"PA" : -63.17, "crpix1" : 1064.9, "crpix2" : 1099.1, "sign" : 1, "pxscl" : 13, "raoff" : -2.73765e-4, "decoff" : 1.812239e-4},
            "27mas-JHK" : {"PA" : -65.02, "crpix1" : 958.0, "crpix2": 1065.9, "sign" : -1, "pxscl" : 27, "raoff" : -2.4458e-4, "decoff" : 1.3427e-4}
            }
         }
"""

SETUPS = {
         "NIX" : {
            "13mas-JHK" : {"PA" : -0.4, "crpix1" : 1064.9, "crpix2" : 1099.1, "sign" : 1, "pxscl" : 13.09, "raoff" : 0, "decoff" : 0},
            "27mas-JHK" : {"PA" : 1.55, "crpix1" : 958.0, "crpix2": 1065.9, "sign" : -1, "pxscl" : 26.92, "raoff" : 0, "decoff" : 0},
            "13mas-LM" : {"PA" : -0.4, "crpix1" : 1064.9, "crpix2" : 1099.1, "sign" : 1, "pxscl" : 13.03, "raoff" : 0, "decoff" : 0}
            },
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

    return INS, config

def update_wcs(path, f_name, out=None):

    INS, config = get_config(path, f_name)

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

    if INS == "SPIFFIER": # can only do for SPIFFIER
        inv = np.linalg.inv(cd)
        cumoffsky = np.array([head.get("ESO OCS CUMOFFS RA"), head.get("ESO OCS CUMOFFS DEC")])/60./60.
        cumoffpix = -np.matmul(inv, cumoffsky)

        head.set("CUM X", cumoffpix[0])
        head.set("CUM Y", cumoffpix[1])

    if inplace:
        fits.PrimaryHDU(data, header=head).writeto("%s/%s" % (path, f_name), overwrite=True)
    else:
        if out is None:
            fits.PrimaryHDU(data, header=head).writeto("./%s" % f_name, overwrite=True)
        else:
            fits.PrimaryHDU(data, header=head).writeto(out, overwrite=True)

if __name__=="__main__":

    print(f_name)
    update_wcs(path, f_name, out=out)

