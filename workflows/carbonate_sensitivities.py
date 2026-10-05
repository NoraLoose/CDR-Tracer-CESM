#!/usr/bin/env python
# coding: utf-8
"""
Shared surface carbonate-sensitivity calculation.

Used by both
  * compute_sensitivities.py        - control run, global, one 2-D field per month
  * compute_truth_sensitivities.py  - DOR / OAE truth experiments, per polygon

Definition (PyCO2SYS, surface carbonate state):

    beta = (DIC - (HCO3 + 2*CO3) / Q) / CO2       -> stored as dDICdCO2
    eta  = 1 / Q,   Q = isocapnic_quotient         -> stored as dDICdALK
"""

import PyCO2SYS as pyco2

# POP tracer units are mmol/m^3 (meq/m^3 for ALK); PyCO2SYS expects umol/kg.
SEAWATER_DENSITY = 1025.0                       # kg/m^3
MMOL_M3_TO_UMOL_KG = 1000.0 / SEAWATER_DENSITY

# output variable names / documentation strings
BETA_NAME = "dDICdCO2"
ETA_NAME = "dDICdALK"
BETA_FORMULA = "(dic - (HCO3 + 2*CO3) / isocapnic_quotient) / CO2"
ETA_FORMULA = "1 / isocapnic_quotient"


def carbonate_sensitivity(alk, dic, salt, temp, po4, sio3):
    """Return ``(beta, eta)`` for a surface carbonate state.

    Parameters
    ----------
    alk, dic, po4, sio3
        mmol/m^3 (meq/m^3 for ``alk``).
    salt
        practical salinity.
    temp
        temperature, degrees C.

    Inputs are array-likes that broadcast together (scalars, numpy arrays,
    xarray DataArrays); the returned ``beta`` and ``eta`` are numpy arrays of
    the broadcast shape.
    """
    csys = pyco2.sys(
        par1=alk * MMOL_M3_TO_UMOL_KG,
        par2=dic * MMOL_M3_TO_UMOL_KG,
        par1_type=1,
        par2_type=2,
        salinity=salt,
        temperature=temp,
        total_silicate=sio3 * MMOL_M3_TO_UMOL_KG,
        total_phosphate=po4 * MMOL_M3_TO_UMOL_KG,
    )

    beta = (
        csys["dic"] - (csys["HCO3"] + 2 * csys["CO3"]) / csys["isocapnic_quotient"]
    ) / csys["CO2"]

    eta = 1.0 / csys["isocapnic_quotient"]

    return beta, eta
