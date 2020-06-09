#!/bin/bash

python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingMVA94XV2wp90 --checkBins
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingMVA94XV2wp90 --createBins
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingMVA94XV2wp90 --createHists
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingMVA94XV2wp90 --doFit

python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingMVA94XV2wp90 --doFit --mcSig --altSig
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingMVA94XV2wp90 --doFit --altSig
python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingMVA94XV2wp90 --doFit --altBkg

python tnpEGM_fitter.py etc/config/settings_pho_UL2018.py --flag passingMVA94XV2wp90 --sumUp
