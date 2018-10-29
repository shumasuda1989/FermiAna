from ROOT import TFile, TH2F, TGraph, TGraphErrors, TGraphAsymmErrors, TMultiGraph, TLegend, TArrow, TList
from ROOT import gStyle, gPad, gROOT
from ROOT import kRed, kBlue, kGreen, kMagenta, kCyan, kPink
import numpy as np

class DataContainer():
    def __init__(self,file,scale=1.,leg='',col=kRed,sty=20,title=''):
        self.file=file
        self.scale=scale
        self.leg=leg
        self.col=col
        self.sty=sty
        self.title=title

def GetData(file, scale=1. ,sed=True, title='SED',barUL=True):
    GetData.Ng+=1
    g = TGraphErrors(file)
    gUL= TGraphErrors()

    if sed:
        for i in range(g.GetN()):
            g.GetY()[i]=pow(g.GetX()[i],2)*g.GetY()[i]*1e-6/scale
            g.GetEY()[i]=pow(g.GetX()[i],2)*g.GetEY()[i]*1e-6/scale
            g.GetX()[i]*=1e-3
  
    idel=0
    nUL=0
    while idel<g.GetN():
        if g.GetEY()[idel]<0:
            gUL.SetPoint(nUL,g.GetX()[idel],g.GetY()[idel])
            if barUL:
                gUL.SetPointError(nUL,0,g.GetY()[idel]*1e-5)
            nUL+=1
            g.RemovePoint(idel)
        else: idel+=1

    if sed:
        g.SetTitle(title+";Energy [GeV];E^{2}dN/dE [TeV cm^{-2} s^{-1}]")
    else:
        g.SetTitle(title+";Energy [MeV];dN/dE [cm^{-2} s^{-1} MeV^{-1}]")

    g.SetLineColor(kRed)
    g.SetMarkerColor(kRed)
    g.SetMarkerStyle(20)
    g.SetName("g%d"%GetData.Ng)
    gUL.SetLineColor(g.GetLineColor())
    gUL.SetName("gUL%d"%GetData.Ng)

    return g, gUL
GetData.Ng=-1


def GetPubData(datafile, errorfile="NONE", f=False):

    xy=np.loadtxt(datafile)
    x=xy[:,0]*1.
    y=xy[:,1]*1.
    if f:  # dN/dE -> E2dN/dE
        y*=pow(x*1e-12,2)

    if errorfile != "NONE":
        ey=np.loadtxt(errorfile)
        eyl=ey[0::2,1]*1.
        eyh=ey[1::2,1]*1.

        for i in range(len(eyh)):
            if eyl[i] > eyh[i]:
                tmp1=eyl[i]
                eyl[i]=eyh[i]
                eyh[i]=tmp1

        if f:  # dN/dE -> E2dN/dE
            eyl*=pow(x*1e-12,2)
            eyh*=pow(x*1e-12,2)

        eyl = y - eyl
        eyh -= y

    # eV -> GeV
    x*=1e-9

    # erg -> TeV
    y/=1.60218

    if errorfile != "NONE":
        eyl/=1.60218
        eyh/=1.60218

    if errorfile != "NONE":
        g = TGraphAsymmErrors(len(x),x,y,np.zeros(len(x)),np.zeros(len(x)),eyl,eyh)
    else:
        g = TGraph(len(x),x,y)

    return g


def DrawUL(gUL):
    alist = TList()
    for i in range(gUL.GetN()):
        a=TArrow(gUL.GetX()[i],gUL.GetY()[i],gUL.GetX()[i],gUL.GetY()[i]/3.,0.01,">")
        a.SetLineColor(gUL.GetLineColor())
        a.Draw()
        alist.Add(a)

    return alist


def SED(flist):

    gStyle.SetOptLogx()
    gStyle.SetOptLogy()
    gStyle.SetPadGridX(True)
    gStyle.SetPadGridY(True)
    gStyle.SetEndErrorSize(10)

    l=TLegend(.65,.6,.96,.95)

    g = []
    gUL = []

    pubdir="/home/smasuda/storage/.Fermi2/Data/DataPoint/GammaCygni/"


    # g.append(GetPubData(pubdir+"GCyg_FermiLande_data.csv",
    #                     pubdir+"GCyg_FermiLande_error.csv"))
    # g[-1].SetName('gLande12')
    # g[-1].SetLineColor(kCyan)
    # g[-1].SetMarkerColor(g[-1].GetLineColor())
    # g[-1].SetMarkerStyle(21)
    # gUL.append(TGraph(0))
    # l.AddEntry(g[-1],'Lande+ \'12','p')

    g.append(GetPubData(pubdir+"GCyg_FermiFraija_data.csv",
                        pubdir+"GCyg_FermiFraija_error.csv"))
    g[-1].SetName('gFraija16')
    g[-1].SetLineColor(kPink+1)
    g[-1].SetLineColor(kCyan)
    g[-1].SetMarkerColor(g[-1].GetLineColor())
    g[-1].SetMarkerStyle(22)
    gUL1=GetPubData(pubdir+"GCyg_FermiFraija_UL.csv")
    gUL1.SetName('gFraija16UL')
    gUL1.SetLineColor(g[-1].GetLineColor())
    gUL1.SetMarkerColor(g[-1].GetLineColor())
    gUL.append(gUL1)
    l.AddEntry(g[-1],'Fraija+ \'16','p')
    # g.append(GetPubData(pubdir+"GCyg_MAGIC_data.csv",
    #                     pubdir+"GCyg_MAGIC_error.csv"))
    # g[-1].SetLineColor(kMagenta)
    # g[-1].SetMarkerColor(kMagenta)
    # g[-1].SetMarkerStyle(29)
    # gUL.append(TGraph(0))

    def ConvertFluxUnit(glist):
        g=glist[-1]
        n=g.GetN()
        x=g.GetX()
        y=g.GetY()
        ey=g.GetEY()
        for i in range(n):
            y[i]*=x[i]*x[i]*1e-6*1e-4 
            if g.ClassName() == "TGraphErrors":
                ey[i]*=x[i]*x[i]*1e-6*1e-4

    rfile=TFile("/home/smasuda/storage/Fermi/Data/GammaCygni2017/VERJ2019p407_spec.root")
    gVER=rfile.Get("gVER")
    gVERUL=rfile.Get("gVERUL")
    fVER=rfile.Get("fVER")
    rfile.Close()
    g.append(gVER)
    gUL.append(gVERUL)
    fVER.SetParameter(0,fVER.GetParameter(0)*1e-4)
    fVER.SetParameter(1,fVER.GetParameter(1)-2.)
    ConvertFluxUnit(g)
    ConvertFluxUnit(gUL)
    colV=kGreen
    gVER.SetMarkerColor(colV)
    gVER.SetLineColor(colV)
    gVERUL.SetMarkerColor(colV)
    gVERUL.SetLineColor(colV)
    fVER.SetLineColor(colV)
    l.AddEntry(g[-1],'VER J2019+407','p')

    rfile=TFile('/home/smasuda/storage/MAGIC/GammaCygni2017_ST7/analysis/flute/CombAll/Unfolding_Output_combunfold_2-Tikhonov.root')
    g.append(rfile.Get("fGraph1E2"))
    rfile.Close()
    g[-1].SetMarkerStyle(21)
    g[-1].SetMarkerColor(kMagenta+1)
    g[-1].SetLineColor(g[-1].GetMarkerColor())
    gUL.append(TGraph(0))
    l.AddEntry(g[-1],'MAGIC this work','lp')

    rfile=TFile("/home/smasuda/storage/Fermi/Data/GammaCygni2017/MarcelFermi.root")
    gtmp=[]
    for key in rfile.GetListOfKeys():
        obj=key.ReadObj()
        if obj.InheritsFrom(TGraph.Class()):
            gtmp.append(obj)

    # g.append(gtmp[0])
    # gUL.append(TGraph(0))
    # g.append(gtmp[1])
    # gUL.append(TGraph(0))
    # l.AddEntry(g[-1],'Marcel disk','lp')
    # g.append(gtmp[2])
    # gUL.append(gtmp[3])
    # l.AddEntry(g[-1],'Marcel gaus','lp')
    # g.append(TGraph(0))
    # gUL.append(gtmp[4])
    # l.AddEntry(gUL[-1],'Marcel arc','lp')
 
    # arr1=rfile.Get("aGausUL")
    # arr2=rfile.Get("aArcUL")

    rfile.Close()



    # XX=np.array([1.])
    # YY=np.array([1.e-10])
    # g.append(TGraph(1,XX,YY))
    # gUL.append(TGraph(0))

    for f in flist:
        tmpg,tmpgUL=GetData(f.file,f.scale,title=f.title)
        tmpg.SetLineColor(f.col)
        tmpg.SetMarkerColor(f.col)
        tmpg.SetMarkerStyle(f.sty)
        tmpgUL.SetLineColor(f.col)
        if f.leg != '':
            l.AddEntry(tmpg,f.leg,'p')
        g.append(tmpg)
        gUL.append(tmpgUL)

    # l.AddEntry('p')

    ng = len(g)

    mg = TMultiGraph()
    for i in range(ng):
        if g[i].GetName() == 'gLP':
            mg.Add(g[i],'3')
        else:
            mg.Add(g[i],'pz')

        if gUL[i].GetN()>0:
            if gUL[i].GetName() == 'gGausUL' or gUL[i].GetName() == 'gArcUL':
                mg.Add(gUL[i],'pz')
            else:
                mg.Add(gUL[i],'p')

    xbin=np.logspace(np.log10(0.9),np.log10(8079.765852892145),101)
    # ybin=np.logspace(np.log10(1.565210892602076e-14),np.log10(1.675199606589398e-10),101)
    ybin=np.logspace(np.log10(3.786556805899183e-14),np.log10(1.05e-10),101)
    frame=TH2F("frame","SED;Energy [GeV];E^{2}dN/dE [TeV cm^{-2} s^{-1}]",100,xbin,100,ybin)
    frame.SetStats(False)
    frame.SetDirectory(0)
    frame.Draw("0")
    mg.SetTitle("SED;Energy [GeV];E^{2}dN/dE [TeV cm^{-2} s^{-1}]")
    mg.Draw()
    fVER.Draw("same")

    arr=[]
    for i in range(len(gUL)):
        if gUL[i].GetN()>0:
            arr.append(DrawUL(gUL[i]))
    # arr.append(arr1)
    # arr.append(arr2)
    # arr1.Draw()
    # arr2.Draw()
    l.Draw()

    gROOT.frame=frame
    gROOT.mg=mg
    gROOT.arr=arr
    gROOT.l=l
    gROOT.fVER=fVER

    gStyle.SetOptLogx(False)
    gStyle.SetOptLogy(False)
