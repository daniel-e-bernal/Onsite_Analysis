"""
This file contains the script to run REopt on the HPC.

This is not meant to run the full optimized techo-economic capabilities of REopt.
We're only doing technical potential so the only real inputs are 
    A) space availability, (constraint)
    B) load profile, (hopefully scaled to the specific industrial sector)
        i. electrical,
        ii. thermal.
    C) lat/long.
"""

using REopt, DataFrames, XLSX, JSON, DeliminatedFiles

#also potentially using PyCall to call on the datasets that are within the HPC storage such as NSRDBX and WINDToolkit
#using PyCall

