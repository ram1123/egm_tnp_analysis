#!/bin/bash

python tnpEGM_fitter.py etc/config/settings_ele_UL2018.py --flag passingMVA94Xwp90isoV2 --checkBins
python tnpEGM_fitter.py etc/config/settings_ele_UL2018.py --flag passingMVA94Xwp90isoV2 --createBins
python tnpEGM_fitter.py etc/config/settings_ele_UL2018.py --flag passingMVA94Xwp90isoV2 --createHists
python tnpEGM_fitter.py etc/config/settings_ele_UL2018.py --flag passingMVA94Xwp90isoV2 --doFit

python tnpEGM_fitter.py etc/config/settings_ele_UL2018.py --flag passingMVA94Xwp90isoV2 --doFit --mcSig --altSig
python tnpEGM_fitter.py etc/config/settings_ele_UL2018.py --flag passingMVA94Xwp90isoV2 --doFit --altSig
python tnpEGM_fitter.py etc/config/settings_ele_UL2018.py --flag passingMVA94Xwp90isoV2 --doFit --altBkg

python tnpEGM_fitter.py etc/config/settings_ele_UL2018.py --flag passingMVA94Xwp90isoV2 --sumUp
