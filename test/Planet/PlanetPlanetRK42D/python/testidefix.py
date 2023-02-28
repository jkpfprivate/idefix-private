#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  21 12:04:41 2021

@author: gwafflard
"""

import os
import sys
import re
sys.path.append(os.getenv("IDEFIX_DIR"))
import numpy as np

def datafile(filename, *, directory=""):
    fullpath = os.path.join(directory, filename)
    with open(fullpath) as f1:
        data = f1.readlines()
    y = [[v for v in re.split(r"[\t ]+", r)] for r in data]
    columns = np.array(y, dtype="float64").T
    return(columns)

plot = False

planet0_ref = datafile("planet0.ref.dat")
planet0 = datafile("../planet0.dat")
planet1_ref = datafile("planet1.ref.dat")
planet1 = datafile("../planet1.dat")

dist0_ref = np.sqrt(planet0_ref[1]**2+planet0_ref[2]**2+planet0_ref[3]**2)
dist0 = np.sqrt(planet0[1]**2+planet0[2]**2+planet0[3]**2)
dist1_ref = np.sqrt(planet1_ref[1]**2+planet1_ref[2]**2+planet1_ref[3]**2)
dist1 = np.sqrt(planet1[1]**2+planet1[2]**2+planet1[3]**2)

# Compute the error on the planet distance
error_d0=np.max(np.abs((dist0-dist0_ref)/dist0_ref))
error_d1=np.max(np.abs((dist1-dist1_ref)/dist1_ref))

torque0_ref = datafile("tqwk0.ref.dat")
torque0 = datafile("../tqwk0.dat")
torque1_ref = datafile("tqwk1.ref.dat")
torque1 = datafile("../tqwk1.dat")

tq0_tot_ref = torque0_ref[1]+torque0_ref[2]
tq0_tot = torque0[1]+torque0[2]
tq1_tot_ref = torque1_ref[1]+torque1_ref[2]
tq1_tot = torque1[1]+torque1[2]

# Compute the error on the planet torque
error_t0=np.max(np.abs((tq0_tot-tq0_tot_ref)/tq0_tot_ref))
error_t1=np.max(np.abs((tq1_tot-tq1_tot_ref)/tq1_tot_ref))


if plot:
    import matplotlib.pyplot as plt
    plt.close("all")
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(planet0_ref[-1], dist0_ref, ls="-", c="k", label="dist 0 ref")
    ax.plot(planet0[-1], dist0, ls="--", c="r", label="dist 0 new")
    ax.plot(planet1_ref[-1], dist1_ref, ls="-", c="b", label="dist 1 ref")
    ax.plot(planet1[-1], dist1, ls="--", c="g", label="dist 1 new")
    ax.set(xlabel="time")
    fig.tight_layout()
    ax.legend(frameon=False)
    fig2, ax2 = plt.subplots(figsize=(6,6))
    ax2.plot(torque0_ref[-1], tq0_tot_ref, ls="-", c="k", label="torque 0 ref")
    ax2.plot(torque0[-1], tq0_tot, ls="--", c="r", label="torque 0 new")
    ax2.plot(torque1_ref[-1], tq1_tot_ref, ls="-", c="b", label="torque 1 ref")
    ax2.plot(torque1[-1], tq1_tot, ls="--", c="g", label="torque 1 new")
    ax2.set(xlabel="time")
    fig2.tight_layout()
    ax2.legend(frameon=False)
    plt.show()

diff_torque0 = np.max(np.abs(tq0_tot-tq0_tot_ref))
diff_torque1 = np.max(np.abs(tq1_tot-tq1_tot_ref))
print("diff_torque0=%e"%diff_torque0)
print("diff_torque1=%e"%diff_torque1)
print("Error_dist0=%e"%error_d0)
# print("Error_torque0=%e"%error_t0)
print("Error_dist1=%e"%error_d1)
# print("Error_torque1=%e"%error_t1)
# if error_d0<5.0e-2 and error_t0<5.0e-2 and error_d1<5.0e-2 and error_t1<5.0e-2 and error_rho<5.0e-2:
if error_d0<5.0e-2 and diff_torque0<1.0e-15 and error_d1<5.0e-2 and diff_torque1<1.0e-15:
    print("SUCCESS!")
    sys.exit(0)
else:
    print("FAILURE!")
    sys.exit(1)
