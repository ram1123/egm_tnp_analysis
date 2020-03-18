### python specific import
import argparse
import os
import sys
import pickle
import shutil
import ROOT as rt
import math

parser = argparse.ArgumentParser(description='tnp EGM fitter')
parser.add_argument('--checkBins'  , action='store_true'  , help = 'check  bining definition')
parser.add_argument('--createBins' , action='store_true'  , help = 'create bining definition')
parser.add_argument('--createHists', action='store_true'  , help = 'create histograms')
parser.add_argument('--sample'     , default='all'        , help = 'create histograms (per sample, expert only)')
parser.add_argument('--altSig'     , action='store_true'  , help = 'alternate signal model fit')
parser.add_argument('--altBkg'     , action='store_true'  , help = 'alternate background model fit')
parser.add_argument('--doFit'      , action='store_true'  , help = 'fit sample (sample should be defined in settings.py)')
parser.add_argument('--mcSig'      , action='store_true'  , help = 'fit MC nom [to init fit parama]')
parser.add_argument('--doPlot'     , action='store_true'  , help = 'plotting')
parser.add_argument('--sumUp'      , action='store_true'  , help = 'sum up efficiencies')
parser.add_argument('--dumpFits'   , action='store_true'  , help = 'dump fits and efficiencies')
parser.add_argument('--dumpSum'    , action='store_true'  , help = 'dump summary')
parser.add_argument('--iBin'       , dest = 'binNumber'   , type = int,  default=-1, help='bin number (to refit individual bin)')
parser.add_argument('--flag'       , default = None       , help ='WP to test')
parser.add_argument('settings'     , default = None       , help = 'setting file [mandatory]')


args = parser.parse_args()

print '===> settings %s <===' % args.settings
importSetting = 'import %s as tnpConf' % args.settings.replace('/','.').split('.py')[0]
print importSetting
exec(importSetting)

### tnp library
import libPython.binUtils  as tnpBiner
import libPython.rootUtils as tnpRoot


if args.flag is None:
    print '[tnpEGM_fitter] flag is MANDATORY, this is the working point as defined in the settings.py'
    sys.exit(0)
    
if not args.flag in tnpConf.flags.keys() :
    print '[tnpEGM_fitter] flag %s not found in flags definitions' % args.flag
    print '  --> define in settings first'
    print '  In settings I found flags: '
    print tnpConf.flags.keys()
    sys.exit(1)

outputDirectory = '%s/%s/' % (tnpConf.baseOutDir,args.flag)

print '===>  Output directory: '
print outputDirectory

tnpBins = pickle.load( open( '%s/bining.pkl'%(outputDirectory),'rb') )

for s in tnpConf.samplesDef.keys():
    sample =  tnpConf.samplesDef[s]
    if sample is None: continue
    setattr( sample, 'tree'     ,'%s/fitter_tree' % tnpConf.tnpTreeDir )
    setattr( sample, 'histFile' , '%s/%s_%s.root' % ( outputDirectory , sample.name, args.flag ) )

sampleToFit = tnpConf.samplesDef['data']

if sampleToFit is None:
    print '[tnpEGM_fitter, prelim checks]: sample (data or MC) not available... check your settings'
    sys.exit(1)

sampleMC = tnpConf.samplesDef['mcNom']
if sampleMC is None:
    print '[tnpEGM_fitter, prelim checks]: MC sample not available... check your settings'
    sys.exit(1)
for s in tnpConf.samplesDef.keys():
    sample =  tnpConf.samplesDef[s]
    if sample is None: continue
    setattr( sample, 'mcRef'     , sampleMC )
    setattr( sample, 'nominalFit', '%s/%s_%s.nominalFit.root' % ( outputDirectory , sample.name, args.flag ) )
    setattr( sample, 'altSigFit' , '%s/%s_%s.altSigFit.root'  % ( outputDirectory , sample.name, args.flag ) )
    setattr( sample, 'altBkgFit' , '%s/%s_%s.altBkgFit.root'  % ( outputDirectory , sample.name, args.flag ) )

####################################################################
##### dumping fit details as csv file
####################################################################
if args.dumpFits:
    sampleToFit.dump()
    fileName = None
    fitType  = None
    if args.doFit:
        fileName = sampleToFit.nominalFit
        fitType  = 'nominalFit'
    if args.altSig : 
        fileName = sampleToFit.altSigFit
        fitType  = 'altSigFit'
    if args.altBkg : 
        fileName = sampleToFit.altBkgFit
        fitType  = 'altBkgFit'
    
    if fileName is None or fitType is None:
        print '[dumpFits.py, prelim checks]: File is not specified ... Please check your arguments'
        sys.exit(1)

    if not os.path.exists( fileName ):
        print 'Root file %s does not exist.' % fileName
        sys.exit(1)

    print 'opening ', fileName
    rootfile = rt.TFile(fileName,"read")
    
    datafile = rt.TFile( sampleToFit.histFile, "read")

    EffFile = open("%s/%s_%s_dumpFits.csv" % (outputDirectory,sampleToFit.name,fitType) ,"w")
    EffFile.write("Bin; nSigP; nSigPerr; nSigF; nSigFerr; Eff; Efferr; minNllP; minNllF; StatusP; StatusF; devP; devF; dEff; \n")

    for ib in range(len(tnpBins['bins'])):
        if (ib % 10) ==2 or (ib % 10) ==7:
            continue
        tnpBin = tnpBins['bins'][ib]
        if (args.binNumber >= 0 and ib == args.binNumber) or args.binNumber < 0:
            print '  get fit result for bin %02d: %s_resP/F' % (ib, tnpBin['name'])
            resP  = rootfile.Get( '%s_resP' % tnpBin['name'] )
            resF  = rootfile.Get( '%s_resF' % tnpBin['name'] )
            print '  get canvas for bin %02d: %s_Canv' % (ib, tnpBin['name'])
            canv  = rootfile.Get( '%s_Canv' % tnpBin['name'] )
            nSigP = resP.floatParsFinal().find("nSigP")
            nSigF = resF.floatParsFinal().find("nSigF")
            canP = canv.GetPrimitive("c_2")
            pdfP = None
            for prim in canP.GetListOfPrimitives():
                if 'pdf' in prim.GetName() and 'bkg' not in prim.GetName() :
                    pdfP = prim
            canF = canv.GetPrimitive("c_3")
            pdfF = None
            for prim in canF.GetListOfPrimitives():
                if 'pdf' in prim.GetName() and 'bkg' not in prim.GetName() :
                    pdfF = prim

            x=rt.double(0)
            y=rt.double(0)
            peakP = 0
            for i in range(pdfP.GetN()):
                pdfP.GetPoint(i,x,y)
                if y > peakP:
                    peakP = float(y)
                    
            peakF = 0
            for i in range(pdfF.GetN()):
                pdfF.GetPoint(i,x,y)
                if y > peakF:
                    peakF = float(y)
            print 'peakP = %f ; peakF = %f' % (peakP, peakF) 
            hPass = datafile.Get('%s_Pass' % tnpBin['name'] )
            hFail = datafile.Get('%s_Fail' % tnpBin['name'] )
            devP  = 0
            devF  = 0
            for i in range(hPass.GetXaxis().GetNbins()+1) :
                if hPass.GetXaxis().GetBinLowEdge(i)<75 or hPass.GetXaxis().GetBinLowEdge(i) >105 :
                    continue
                if hPass.GetBinContent(i) > peakP : 
                    devP += hPass.GetBinContent(i) - peakP
                    print 'PeakP at %f underestimated for bin %d at %f ' % (peakP, i, hPass.GetBinContent(i))
                if hFail.GetBinContent(i) > peakF : 
                    devF += hFail.GetBinContent(i) - peakF
                    print 'PeakF at %f underestimated for bin %d at %f ' % (peakF, i, hFail.GetBinContent(i))
            nP = nSigP.getVal()
            e_nP = nSigP.getError()
            nF = nSigF.getVal()
            e_nF = nSigF.getError()
            nTot = nP+nF
            eff = nP/nTot
            e_eff = 1/(nTot*nTot)*math.sqrt(nP*nP* e_nF*e_nF + nF*nF * e_nP*e_nP )
            EffFile.write("%02d; %f; %f; %f; %f; %f; %f; %f; %f; %1d; %1d; %f; %f; %f;\n" % (ib, nP, e_nP, nF, e_nF, eff, e_eff, resP.minNll(), resF.minNll(), resP.status(),resF.status(), devP, devF, (nP+devP)/(nTot+devP+devF) - eff ))
    EffFile.close()

    print ' ===> CSV file exported <======='

####################################################################
##### dumping fit details as csv file
####################################################################
if args.dumpSum:
    sampleToFit.dump()
    info = {
        'data'        : sampleToFit.histFile,
        'dataNominal' : sampleToFit.nominalFit,
        'dataAltSig'  : sampleToFit.altSigFit ,
        'dataAltBkg'  : sampleToFit.altBkgFit ,
        'mcNominal'   : sampleToFit.mcRef.histFile,
        'mcAlt'       : None,
        'tagSel'      : None
        }

    if not tnpConf.samplesDef['mcAlt' ] is None:
        info['mcAlt'    ] = tnpConf.samplesDef['mcAlt' ].histFile
    if not tnpConf.samplesDef['tagSel'] is None:
        info['tagSel'   ] = tnpConf.samplesDef['tagSel'].histFile

    effis = None
    SumFileName ='%s/%s_dumpSum.csv' % (outputDirectory, sampleToFit.name)
    fOut = open( SumFileName,'w')
    fOut.write("Bin; EffData; EffMC; systData; systMC; systAltBkg; systAltSig; systAltMC; systAltTagSelec; totalSyst; \n")
    
    for ib in range(len(tnpBins['bins'])):
        if (ib % 10) ==2 or (ib % 10) ==7:
            continue
        effis = tnpRoot.getAllEffi( info, tnpBins['bins'][ib] )

        systAltBkg      = effis['dataAltBkg' ][0] - effis['dataNominal'][0]
        systAltSig      = effis['dataAltSig' ][0] - effis['dataNominal'][0]
        systAltMC       = effis['mcAlt' ][0] - effis['mcNominal'  ][0]
        systAltTagSelec = effis['tagSel'][0] - effis['mcNominal'  ][0]

        totalSyst = 0
        totalSyst += effis['dataNominal'][1]* effis['dataNominal'][1]
        totalSyst += effis['mcNominal'][1]* effis['mcNominal'][1]
        totalSyst += systAltBkg*systAltBkg
        totalSyst += systAltSig*systAltSig
        totalSyst += systAltMC*systAltMC
        totalSyst += systAltTagSelec*systAltTagSelec
        totalSyst = math.sqrt(totalSyst)

        fOut.write( '%02d; %f; %f; %f; %f; %f; %f; %f; %f; %f; \n' %(ib, effis['dataNominal'][0], effis['mcNominal'  ][0], effis['dataNominal'][1], effis['mcNominal'  ][1], systAltBkg, systAltSig, systAltMC, systAltTagSelec, totalSyst))

    fOut.close()



