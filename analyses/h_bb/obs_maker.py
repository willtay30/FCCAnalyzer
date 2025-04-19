import uproot
import subprocess
import os
import json
import ROOT 
import math

#####
# Warning!
# Run using python, source with only jsetup.sh (ROOT)

# Purpose is to take the combined root file (after running rootmerge.py) and format the file to allow for datacard creation (COMBINE)
# Creates histograms of the total events per HZ subprocess (data_obs) and new histograms for the last bin of the cutflow(_last)
# Also copies root file and and puts a version in the datacard folder

# please write output of this script to a file as information needs to be referred to for cardwriter.py
####

output_file = f"end_of_h_ssApril7.root" #root file

# below if some processes we are ignoring
bad = ["wzp6_ee_mumu_ecm240", "wzp6_ee_tautau_ecm240", "wzp6_egamma_eZ_Zmumu_ecm240", "wzp6_gammae_eZ_Zmumu_ecm240", "wzp6_gaga_mumu_60_ecm240", "wzp6_gaga_tautau_60_ecm240", "wzp6_ee_nuenueZ_ecm240"]

hists = {}
maxNumCuts = 100

fOut = ROOT.TFile(output_file, "READ")

for key in fOut.GetListOfKeys():
    #print(key)
    if key.GetName() in bad:
        #print("         BAD")
        continue

    process = fOut.Get(key.GetName())
    for jey in process.GetListOfKeys():

        #print("j:", jey)
        event = fOut.Get(f"{key.GetName()}/{jey.GetName()}")
        #print(event.GetName())
        if event.Class_Name() == "TH1D":
            name = event.GetName().split("_")
            
            if name[0] == "cutFlow":
                #print("got cutflow: ", name[1])
                if name[1] in hists:
                    #print("in hist")
                    for i in range(maxNumCuts):
                        hists[name[1]].SetBinContent(i, hists[name[1]].GetBinContent(i) + event.GetBinContent(i))
                    
                    #print(hists[name[1]].GetBinContent(1))
                else:
                    #print("not in hist")
                    h_new = ROOT.TH1D("data_obs", "", event.GetNbinsX(), event.GetBinLowEdge(1), event.GetBinLowEdge(event.GetNbinsX() + 1))
                    for i in range(maxNumCuts):
                        h_new.SetBinContent(i, event.GetBinContent(i))
                    
                    hists[name[1]] = h_new
                    
                    
    
binmap = {}
print("****************************************")
print(hists)
for hist in hists.keys():
    print("-----", hist, "-----")
    print(type(hists[hist]))
    h_new = ROOT.TH1D("data_obs", "", hists[hist].GetNbinsX(), hists[hist].GetBinLowEdge(1), hists[hist].GetBinLowEdge(hists[hist].GetNbinsX() + 1))

    for i in range(maxNumCuts):
        hists[hist].SetBinError(i, math.sqrt(hists[hist].GetBinContent(i)))
        print(i, hists[hist].GetBinContent(i), hists[hist].GetBinError(i))
        if hists[hist].GetBinContent(i) == 0 and i != 0:
             print("laster: ", i-1, " | ", hists[hist].GetBinContent(i-1))
             h_new.SetBinContent(1, hists[hist].GetBinContent(i - 1))
             h_new.SetBinError(1, math.sqrt(hists[hist].GetBinContent(i - 1)))
             binmap[hist] = i - 1
             break
    print("new obhist: ", h_new.GetBinContent(1), h_new.GetBinError(1))
    cOut = ROOT.TFile(output_file, "UPDATE")
    cOut.cd()
    cOut.mkdir(hist)
    cOut.cd(hist)
    h_new.Write()
    cOut.cd()
    cOut.Close() 

print("****************************************")
#print(hists)
print(binmap)
print("****************************************")
notzero = {}
for key in fOut.GetListOfKeys():
    #print(key)
    if key.GetName() in bad:
        #print("     BAD")
        continue

    process = fOut.Get(key.GetName())
    for jey in process.GetListOfKeys():

        #print("j:", jey)
        event = fOut.Get(f"{key.GetName()}/{jey.GetName()}")
        #print(event.GetName())
        if event.Class_Name() == "TH1D":
            name = event.GetName().split("_")
            
            if name[0] == "cutFlow":
                
                h_new = ROOT.TH1D(f"{event.GetName()}_last", "", event.GetNbinsX(), event.GetBinLowEdge(1), event.GetBinLowEdge(event.GetNbinsX() + 1))
                h_new.SetBinContent(1, event.GetBinContent(binmap[name[1]]))
                h_new.SetBinError(1, event.GetBinError(binmap[name[1]]))
                
                if h_new.GetBinContent(1) != 0:
                    print("got process: ", key.GetName(), "| cutflow procduct: ", name[1])
                    print("                 new obhist: ", h_new.GetName(), h_new.GetBinContent(1), h_new.GetBinError(1))
                    if event.GetName() in notzero:
                        notzero[event.GetName()].append(f"{key.GetName()}")
                    else:
                        notzero[event.GetName()] = [f"{key.GetName()}"]
                    

                cOut = ROOT.TFile(output_file, "UPDATE")
                cOut.cd()
                cOut.cd(key.GetName())
                h_new.Write()
                cOut.cd()
                cOut.Close() 


print("****************************************")
print("not zeros")
for key in notzero.keys():
    print("--------------", key, "-----------")
    print(notzero[key])

os.system(f"cp {output_file} datacards/{output_file}")