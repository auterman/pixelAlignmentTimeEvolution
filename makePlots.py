#!/usr/bin/env python2

import collections
import glob
import re
import os
import copy
import numpy
import datetime
import string
import shutil
import subprocess
import pickle
import math

import suppressor
with suppressor.suppress_stdout_stderr(): import ROOT
import downloadViaJson
import style

#user specific stuff like workdirs
##import userEnvironment.config 
import imp

from operator import itemgetter
from array import array

userEnvironment = imp.load_source('userEnvironment.config', '.')

class Parameter:
    name = ""
    filename = ""
    label = ""
    cut = 0
    minDraw = 0
    maxDraw = 0
    def __init__(self, n, f, l, c, minDraw, maxDraw):
        self.name = n
	self.filename = f
        self.label = l
        self.cut = c
        self.minDraw = minDraw
        self.maxDraw = maxDraw
parameters = [
    Parameter("Xpos", "Xpos",    "#Deltax (#mum)", 5, -20, 45 ), \
    Parameter("Ypos", "Ypos",    "#Deltay (#mum)", 10, -30, 50 ), \
    Parameter("Zpos", "Zpos",    "#Deltaz (#mum)", 15, -50, 140 ), \
    Parameter("Xrot", "Xrot",    "#Delta#theta_{x} (#murad)", 30, -50, 80 ), \
    Parameter("Yrot", "Yrot",    "#Delta#theta_{y} (#murad)", 30, -50, 80 ), \
    Parameter("Zrot", "Zrot",    "#Delta#theta_{z} (#murad)", 30, -70, 100 ), \
    ]
parDict = collections.OrderedDict( (p.name, p) for p in parameters )
objects = [
    ("FPIX(x+,z-)", ROOT.kBlack,   20),
    ("BPIX(x+)",    ROOT.kBlue,    21),
    ("FPIX(x+,z+)", ROOT.kGreen+2, 22),
    ("FPIX(x-,z-)", ROOT.kRed,     24),
    ("BPIX(x-)",    ROOT.kCyan,    25),
    ("FPIX(x-,z+)", ROOT.kMagenta, 26),
]

plotDir = "/afs/cern.ch/user/a/auterman/public/plots"

class Temp:
    lastTime = ""
    lastRun = 0

def save(name, folder="plots", endings=[".pdf"]):
    for ending in endings:
        ROOT.gPad.GetCanvas().SaveAs(os.path.join(folder,name+ending))

def randomName():
    """
    Generate a random string. This function is useful to give ROOT objects
    different names to avoid overwriting.
    """
    from random import randint
    from sys import maxint
    return "%x"%(randint(0, maxint))

def runFromFilename(filename):
    m = re.match(".*Run(\d+).root", filename)
    if m:
        return int(m.group(1))
    m2 = re.match(".*/Results(\d+)[^\d].*", filename) # for pseudo
    if m2:
        return int(m2.group(1))
    else:
        print "Could not find run number for file", filename
        return 0

def getFromFile(filename, objectname):
    f = ROOT.TFile(filename)
    if f.GetSize()<5000: # DQM files sometimes are empty
        return None
    h = f.Get(objectname)
    h = ROOT.gROOT.CloneObject(h)
    return h

def sortedDict(d):
    return collections.OrderedDict(sorted(d.items(), key=lambda t: t[0]))

def getInputHists(searchPath="root-files/Run*.root"):
    hists = {}
    for filename in glob.glob(searchPath):
        runNr = runFromFilename(filename)
        newHists = {}
        if searchPath.endswith("PCL_SiPixAl_DQM.root"):
            c = getFromFile(filename, "PCL_SiPixAl_Expert")
            for pad in c.GetListOfPrimitives():
                pad.cd()
                for x in pad.GetListOfPrimitives():
                    if isinstance(x, ROOT.TH1F):
                        newHists[x.GetName()] = x.Clone()
                        break
        else: # dqm plots
            for p in parameters:
                h = getFromFile(filename, p.name)
                if h:
                    newHists[p.name] = h
        if newHists: hists[runNr] = newHists
    return sortedDict(hists)

def exceedsCuts(h, cutDict=False):
    maxErrCut = 10
    sigCut = 2.5
    maxCut = 200
    var = h.GetName().split("_")[0]
    cut = parDict[var].cut

    binInfos = []
    for bin in range(1,h.GetNbinsX()+1):
        c = abs(h.GetBinContent(bin))
        e = h.GetBinError(bin)
        if c > maxCut or e > maxErrCut:
            binInfos.append("fail")
        elif c > cut and e and c/e > sigCut:
            binInfos.append("update")
        else:
            binInfos.append("good")
    if "fail" in binInfos:
        return "fail"
    elif "update" in binInfos:
        return "update"
    else:
        return "good"

def getRunEndTime(run):
    #returs a string similar to 2016-06-16 23:30:32
    return subprocess.check_output(["das_client.py --limit=0 --query=\"run={} | grep run.end_time\"".format(run)], shell=True)

def getValidRunBefore(run):
    foundRun = False
    run -= 1
    while not foundRun:
        try:
            out = subprocess.check_output(["das_client.py --limit=0 --query=\"run={} | grep run.end_time\"".format(run)], shell=True)
            foundRun = True
        except:
            run -=1
    return run


def getFieldFromDB(run):
    bfield = subprocess.check_output(["das_client", "--limit", "0",  "--query", "run={} | grep run.bfield".format(run)])
    bfieldvals=re.findall("\d+\.\d+", bfield)
    return -1 if len(bfieldvals) == 0 else bfieldvals[0]

def getField(run, dbName="runField.pkl"):
    db = {}
    if os.path.exists(dbName):
        with open(dbName) as f:
            db = pickle.load(f)
    if run not in db or db[run] == "\n":
    	db[run] = getFieldFromDB(run)
	print "Got B-Field for run {}: {}".format(run, db[run])
    with open(dbName, "wb") as f:
        pickle.dump(db, f)
    return db[run]


def getLuminosity(minRun):
    """Expects something like
    +-------+------+--------+--------+-------------------+------------------+
    | nfill | nrun | nls    | ncms   | totdelivered(/fb) | totrecorded(/fb) |
    +-------+------+--------+--------+-------------------+------------------+
    | 73    | 327  | 142418 | 138935 | 19.562            | 18.036           |
    +-------+------+--------+--------+-------------------+------------------+
    And extracts the total recorded luminosity (/fb).
    """
    output = subprocess.check_output(["/afs/cern.ch/user/a/auterman/.local/bin/brilcalc", "lumi", "-b", "STABLE BEAMS", "--normtag=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_BRIL.json", "-u", "/fb", "--begin", str(minRun)])
    return float(output.split("\n")[-4].split("|")[-2])

def getTime(run, dbName="runTime.pkl"):
    db = {}
    if os.path.exists(dbName):
        with open(dbName) as f:
            db = pickle.load(f)
    if not run in db:
        db[run] = getRunEndTime(run)
	db[run] = db[run].replace('"','')
    #print "run ",run," >>>",db[run],"<<<"
    lastRun=run
    while db[run] == "\n": 
        lastRun = getValidRunBefore(lastRun)
        db[run] = getRunEndTime(lastRun)
        db[run] = db[run].replace('"','')
        print "run ",run," ->Get Time for previous run ",lastRun,": db=",db[run]
    with open(dbName, "wb") as f:
        pickle.dump(db, f)
    return db[run]

def sendMail(adress, subject="", body=""):
    os.system("echo \"{}\" | mail -s \"{}\" {}".format(body, subject, adress))

def drawHists(hmap, savename, run):
    hnames = ["Xpos", "Ypos","Zpos", "Xrot", "Yrot", "Zrot"]
    line = ROOT.TLine()
    line.SetLineColor(ROOT.kRed)
    c = ROOT.TCanvas(randomName(),"",1200,600)
    c.Divide(3,2)
    dbUpdated = False
    for ih, hname in enumerate(hnames):
        c.cd(ih+1)
        h = hmap[hname]
        h.SetLineColor(ROOT.kBlack)
        h.SetFillColor(ROOT.kGreen-7)
        cutStatus = exceedsCuts(h)
        if cutStatus == "update":
            h.SetFillColor(ROOT.kOrange-9)
            dbUpdated = True
        elif cutStatus == "fail":
            h.SetFillColor(ROOT.kRed)
        for bin in range(1,7):
            h.GetXaxis().SetBinLabel(bin,objects[bin-1][0])
        h.GetXaxis().SetRange(1,6)
        h.GetYaxis().SetRangeUser(-50,50)
        h.SetTitle("")
        h.GetYaxis().SetTitle(parameters[ih].label)
        h.Draw("histe")
        cut = h.GetBinContent(8)
        if not cut:
            cuts = {"Xpos":5, "Ypos":10, "Zpos":15, "Xrot":30, "Yrot":30, "Zrot":30}
            cut = cuts[h.GetName().split("_")[0]]
        line.DrawLine(0,-cut,6,-cut)
        line.DrawLine(0,+cut,6,+cut)
    c.cd(0)
    text = ROOT.TLatex()
    text.SetTextSize(.75*text.GetTextSize())
    text.DrawLatexNDC(.05, .967, "#scale[1.2]{#font[61]{CMS}} #font[52]{Preliminary}")
    text.DrawLatexNDC(.82, .967, "Run {} (13TeV)".format(run))
    #text.DrawLatexNDC(.15, .899, "#scale[0.85]{Tracker alignment in 2016 data-taking used as reference}")
    save(savename, plotDir, [".pdf",".png", ".root"])
#    if dbUpdated:
#        sendMail("kiesel@cern.ch auterman@cern.ch", "[PCL] Cuts exceeded", "Run: {}\nSee cern.ch/test-cmsPixelAlignmentSurveillance".format(run))


def drawPublicStyleHists(hmap, savename, run):
    hnames = ["Xpos", "Ypos","Zpos", "Xrot", "Yrot", "Zrot"]
    line = ROOT.TLine()
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(2)
    c = ROOT.TCanvas(randomName(),"",600,600)
    text = ROOT.TLatex()
    text.SetTextSize(.66*text.GetTextSize())
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.12)
    for ih, hname in enumerate(hnames):
        h = hmap[hname]
        h.SetLineColor(ROOT.kBlack)
        h.SetFillColor(ROOT.kGreen-7)
        cutStatus = exceedsCuts(h)
        if cutStatus == "update":
            h.SetFillColor(ROOT.kOrange-9)
            dbUpdated = True
        elif cutStatus == "fail":
            h.SetFillColor(ROOT.kRed)
        for bin in range(1,7):
            h.GetXaxis().SetBinLabel(bin,objects[bin-1][0])
        h.GetXaxis().SetRange(1,6)
        h.GetYaxis().SetRangeUser(-50,50)
        h.SetTitle("")
        h.GetYaxis().SetTitle(parameters[ih].label)
	h.GetYaxis().SetTitleOffset(1.2)
        h.Draw("histe")
        cut = h.GetBinContent(8)
        if not cut:
            cuts = {"Xpos":5, "Ypos":10, "Zpos":15, "Xrot":30, "Yrot":30, "Zrot":30}
            cut = cuts[h.GetName().split("_")[0]]
        line.DrawLine(0,-cut,6,-cut)
        line.DrawLine(0,+cut,6,+cut)
        text.DrawLatexNDC(.15, .95, "#scale[1.2]{#font[61]{CMS}} #font[52]{Preliminary}")
        text.DrawLatexNDC(.5, .95, "Run {} (13TeV, 2016)".format(run))
        #text.DrawLatexNDC(.15, .899, "#scale[0.85]{Tracker alignment in 2016 data-taking used as reference}")
        save(savename+hname, plotDir, [".pdf",".png", ".root"])



def drawGraphsVsX(gmap, xaxis, savename, magnetGraph, params, objcts, specialRuns=[]):
    """ Options for xaxis: time, run"""
    #lumi = getLuminosity(273000)
    lumi = getLuminosity(278888)
    print "Int. lumi. of plotted run range is: ",lumi
    if not gmap: return
    line = ROOT.TLine()
    line.SetLineColor(ROOT.kBlack)
    blackline = ROOT.TLine()
    blackline.SetLineColor(ROOT.kBlack)
    updateLine = ROOT.TLine()
    updateLine.SetLineStyle(2)
    updateLine.SetLineColor(18)
    leg = ROOT.TLegend(.2, .65, .8, .928)
    leg.SetHeader("Tracker alignment in 2016 data-taking used as reference")
    leg.SetNColumns(3)
    leg.AddEntry(line, "Update threshold", "l")
    leg.AddEntry(updateLine, "Alignment update", "l")
    leg.AddEntry(magnetGraph, "Magnet < 3.8 T", "f")
    for ip, p in enumerate(params):
        c = ROOT.TCanvas(randomName(),"",1200,600)
	for io,o in enumerate(objcts):
	    g = gmap[p.name][o[0]]
            if xaxis == "time":
                g.SetTitle(";Time;{}".format(p.label))
                g.GetXaxis().SetTimeDisplay(1)
                g.GetXaxis().SetTimeOffset(25)
                g.GetXaxis().SetTimeFormat("%Y-%m-%d")
                g.GetXaxis().SetNdivisions(6,0,0)
            elif xaxis == "run":
                g.SetTitle(";Run;{}".format(p.label))
                g.GetXaxis().SetNoExponent()
                g.GetXaxis().SetNdivisions(7,0,0)
            else:
                print "No idea what to do with x-axis", xaxis
            g.GetYaxis().SetRangeUser(p.minDraw, p.maxDraw)
            g.SetMarkerColor( o[1])
            g.SetLineColor(   o[1])
	    g.SetMarkerStyle( o[2])
	    g.SetMarkerSize(1)
	    g.SetLineWidth(2)
            g.Draw( "same pl" if io>0 else "ap")
	    if not io:
		magnetGraph.Draw("same, e3")  
		xax = g.GetXaxis()
                xmin, xmax = xax.GetXmin(), xax.GetXmax()
		blackline.DrawLine(xmin, 0, xmax, 0)  
                g.Draw( "same pl" )
            if not ip:
	        leg.AddEntry(g, o[0], "lp")
	    ##end: for g,ig
        line.DrawLine(xmin, -p.cut, xmax, -p.cut)
        line.DrawLine(xmin, +p.cut, xmax, +p.cut)
        for r in specialRuns:
            updateLine.DrawLine(r, p.minDraw, r, p.maxDraw)
        text = ROOT.TLatex()
        text.DrawLatexNDC(.08, .945, "#scale[1.2]{#font[61]{CMS}} #font[52]{Preliminary}")
        text.DrawLatexNDC(.56, .945, "13 TeV data  (Aug. 16 - Dec. 5, 2016)")
        #text.DrawLatexNDC(.79, .945, "{:.1f} fb^{{-1}} (13TeV)".format(lumi))
        #if ip == 0: 
	leg.Draw()
        ROOT.gPad.RedrawAxis()
	save(savename+"_"+p.filename, plotDir, endings=[".pdf",".png", ".root"])


def string2Time(timeStr):
    #print timeStr
    return ROOT.TDatime(timeStr).Convert(0)

def getGraphsVsRun(inputHists, minRun=-1, convertToTime=False):
    inputHists = sortedDict(dict((key,value) for key, value in inputHists.iteritems() if key >= minRun))
    ##gdefault = ROOT.TGraphErrors()
    gdefault = ROOT.TGraph()
    graphsVsRun = {}
    #showFirstAfter281000 = False
    for iRun, (runNr, hmap) in enumerate(inputHists.iteritems()):
        xVar = string2Time(getTime(runNr)) if convertToTime else runNr
        for hname, h in hmap.iteritems():
            if hname not in graphsVsRun: graphsVsRun[hname] = dict([(obj[0],gdefault.Clone()) for obj in objects ])
            for bin in range(1,7):
                c = h.GetBinContent(bin)
		#if c == 0: continue
                e = h.GetBinError(bin)
                if abs(c) < 1e-15 or abs(e) > 50: continue
                n = graphsVsRun[hname][objects[bin-1][0]].GetN()
		##e.g. graphsVsRun[Xpos][FPIX(x+,z-)]
                graphsVsRun[hname][objects[bin-1][0]].SetPoint(n, xVar, c)
		#if not showFirstAfter281000 and not convertToTime and xVar>284658:
		#   print "First Run After Magnet back at 3.8T: ",xVar
		#   showFirstAfter281000 = 1
    return graphsVsRun

def updateFile(source, dest, changes={}):
    with open(source) as f:
        x = string.Template(f.read())
    with open(dest, "w") as f:
        f.write(x.safe_substitute(changes))

def getNthLastRun(inputHists, N):
    sortedRuns = sorted(inputHists.keys())
    return sortedRuns[-min(N, len(sortedRuns))]

def isFilledRun(hmap):
    globalMax = 0
    for k, v in hmap.iteritems():
        for bin in range(1,7):
            globalMax = max(globalMax, v.GetBinContent(bin))
    return abs(globalMax) > 1e-6

def getTableString(inputHists, maxPlots=5):
    inputHists = collections.OrderedDict(reversed(list(inputHists.items())))
    outString = "<table>\n<tr> <td> Run </td> <td> End time </td> <td>Parameters</td> </tr>"
    for run, hmap in inputHists.iteritems():
        link = "No results"
        if isFilledRun(hmap):
            if maxPlots > 0:
                link = "<a href=plots/Run{0}.pdf><img src='plots/Run{0}.png' border='0'/></a>".format(run)
                maxPlots -= 1
            else:
                link = "<a href=plots/Run{0}.pdf>pdf</a>".format(run)
        outString += "\n<tr> <td>{0}</td> <td>{1}</td> <td>{2}</td> </tr>".format(run, getTime(run), link)
    outString += "\n</table>"
    return outString

def getUpdateRuns(tag):
    out = subprocess.check_output(["conddb", "list", tag])
    return [int(x.split()[0]) for x in out.split("\n")[2:] if x]


def stringToSqlTimeString(timeStr):
    return timeStr.replace(".", "-")

def getRunFromTime(inputTime, dbFile ="runTime.pkl"):
    #time is ROOT::TDatime
    #sqlTimeStr = stringToSqlTimeString(time)
    with open(dbFile) as f:
        db = pickle.load(f)
    timeConvertedDb = dict( [(run, string2Time(time)) for run,time in db.iteritems()])
    for run, time in sorted(timeConvertedDb.iteritems(), key=itemgetter(1)):
        if inputTime < time: return run
    print "No matching run found"
    return 0
    
    
def ReadMagnetFieldHistory(filename, convertToTime=False):
    mgnt = ROOT.TGraphErrors()
    with open(filename) as f:
        content = f.readlines()
        for i, line in enumerate(content):
            timeStr, fieldStr = line.split(",")
	    
            time = string2Time(stringToSqlTimeString(timeStr))
            field = float(fieldStr)
	    if not convertToTime:
	        run = getRunFromTime(time)
		if run!=0:
		   print time, run
		   mgnt.SetPoint(i, run, 0)
	    else:
                mgnt.SetPoint(i, time, 0)
            mgnt.SetPointError(i, 0, 50 if field>3.7 else 0)
    return mgnt

def GetFieldHistoryByHand(convertToTime=False):
   tstr=["2016.05.04 13:35:26","2016.06.24 12:33:31","2016.06.24 12:33:31","2016.06.24 16:25:35","2016.06.24 16:26:24","2016.06.30 02:20:45","2016.06.30 02:20:45","2016.07.02 13:36:46","2016.07.07 09:40:23","2016.07.26 09:32:26","2016.07.26 09:32:26","2016.07.26 19:06:44","2016.07.26 19:06:44","2016.09.01 09:09:43","2016.09.01 09:09:43","2016.09.01 21:49:36","2016.09.01 21:49:36","2016.09.10 04:25:10","2016.09.10 04:25:10","2016.09.16 19:00:09","2016.09.16 19:00:09","2016.10.31 04:57:21","2016.10.31 04:57:21","2016.11.04 18:44:15","2016.11.04 18:44:15","2016.12.05 07:08:57"]
   x  = [ 273150, 275673, 275674, 275674, 275752, 275960, 276044, 276217, 276235, 277420, 277421, 277784, 277785, 279844, 279845, 279865, 279887, 280384, 280385, 281009, 281010, 284044, 284045, 284658, 284659, 286520, 286521 ]
   y  = [      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0]
   ex = [      0,      0,      1,      1,      0,      0,      1,      1,      0,      0,      1,      1,      0,      0,      1,      1,      0,      0,      1,      1,      0,      0,      1,      1,      0,      0,      1]
   ey = [      0,      0,     80,     80,      0,      0,     80,     80,      0,      0,     80,     80,      0,      0,     80,     80,      0,      0,     80,     80,      0,      0,     80,     80,      0,      0,     80]
   time=[]
   for t in tstr: 
      time.append(string2Time(stringToSqlTimeString(t)))
      
   if not convertToTime:
       ge = ROOT.TGraphErrors(len(x), array('d', x), array('d', y), array('d', ex), array('d', ey))
   else:
       ge = ROOT.TGraphErrors(len(time), array('d', time), array('d', y), array('d', ex), array('d', ey))       
   ge.SetFillColor(17)
   ge.SetLineColor(ROOT.kWhite)
   return ge;


if __name__ == "__main__":
    #downloadViaJson.getGridCertificat()
    #downloadViaJson.downloadViaJson()
    inputHists = getInputHists()

    # draw new runs:
    alreadyPlotted = [ int(x[3:9]) for x in os.listdir(plotDir) if x.endswith(".pdf") and x.startswith("Run")]
    for run, hmap in inputHists.iteritems():
        if run not in alreadyPlotted:
            drawHists(hmap, "Run{}".format(run), run)
        #if run == 285090: drawPublicStyleHists(hmap, "Run285090", 285090)
        #if run == 285216: drawPublicStyleHists(hmap, "Run285216", 285216)

    filename = "MagnetHistory.txt"
    #magnetGraphvsTime = ReadMagnetFieldHistory(filename, convertToTime=True)
    #magnetGraphvsRun  = ReadMagnetFieldHistory(filename, convertToTime=False)
    magnetGraphvsRun   = GetFieldHistoryByHand(convertToTime=False)
    magnetGraphvsTime  = GetFieldHistoryByHand(convertToTime=True)


    simple_parameters = [
        Parameter("Xpos", "Xsimple", "#Deltax (#mum)", 5, -20, 45 ), \
        Parameter("Ypos", "Ysimple", "#Deltay (#mum)", 10, -30, 50 ), \
        Parameter("Zpos", "Zsimple", "#Deltaz (#mum)", 15, -50, 140 )
        ]
    parDict = collections.OrderedDict( (p.name, p) for p in simple_parameters )
    simple_objects = [
        ("BPIX(x+)",    ROOT.kBlue,    21),
        ("BPIX(x-)",    ROOT.kCyan,    25),
    ]


    # vs run
    #updateRuns = [x for x in getUpdateRuns("TrackerAlignment_PCL_byRun_v0_express") if x >= 273000]
    updateRuns = [x for x in getUpdateRuns("TrackerAlignment_PCL_byRun_v0_express") if x >= 278888]
    updateFields = [getField(x) for x in updateRuns]
    graphsVsRun = getGraphsVsRun(inputHists, 278887)
    #drawGraphsVsX(graphsVsRun, "run", "vsRun", magnetGraphvsRun, updateRuns)
    drawGraphsVsX(graphsVsRun, "run", "vsRun", magnetGraphvsRun, parameters, objects, updateRuns)
    drawGraphsVsX(graphsVsRun, "run", "vsRun", magnetGraphvsRun, simple_parameters, simple_objects, updateRuns)

    # vs time
    updateTimes = [string2Time(getTime(x)) for x in updateRuns]
    graphsVsTime = getGraphsVsRun(inputHists, 278887, convertToTime=True)
    drawGraphsVsX(graphsVsTime, "time", "vsTime", magnetGraphvsTime, parameters, objects, updateTimes)
    drawGraphsVsX(graphsVsTime, "time", "vsTime", magnetGraphvsTime, simple_parameters, simple_objects, updateTimes)
#    updateFile("indexTemplate.html", "/afs/cern.ch/user/a/auterman/public/index.html",
#        {
#            "date": datetime.datetime.today().isoformat(' '),
#            "table": getTableString(inputHists)
#        })


