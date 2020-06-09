#!/bin/bash

python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingLoose100XV2 --checkBins
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingLoose100XV2 --createBins
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingLoose100XV2 --createHists
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingLoose100XV2 --doFit

python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingLoose100XV2 --doFit --mcSig --altSig
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingLoose100XV2 --doFit --altSig
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingLoose100XV2 --doFit --altBkg

python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingLoose100XV2 --sumUp
