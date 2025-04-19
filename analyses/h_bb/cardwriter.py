import sys
import os
##################################################################
# Run with python, must source both env_standalone.sh (COMBINE) and jsetup.sh (ROOT) in that order
# sourcing for COMBINE must be done in the HiggsAnalysis/CombinedLimit directory

#   Options
# For each subprocess, check:
# 1) rootfile (made using obs_maker.py)
# 2) hdecay (should match root file)
# 3) cut (Z decay product)
# 4) totalevents (either top of obs_maker.py output or from h_bb.ipynb)
# 5) lumi, acc, and eff (%, from h_bb.ipynb)
# 6) not zero processes (bottom of obs_maker.py output)

# At the bottom, we run a merging script to format the datacard nicely. Change to your path

# all other options are set, but might be changed if creating a different kind of datacard



rootfile = "end_of_h_ssApril7.root"
hdecay = "s" # 1 lowercase
cut = "ss"  # 2 lowercase
totalevents = "76490.2309307271"
# 4 dec places
lumi = "1.0050"
acc = "1.0031"
eff = "1.0054"

notzero = ['wzp6_ee_tautauH_Hbb_ecm240', 'wzp6_ee_qqH_Hbb_ecm240', 'wzp6_ee_ssH_Hbb_ecm240', 'wzp6_ee_ccH_Hbb_ecm240', 'wzp6_ee_eeH_Hcc_ecm240', 'wzp6_ee_tautauH_Hcc_ecm240', 'wzp6_ee_qqH_Hcc_ecm240', 'wzp6_ee_ssH_Hcc_ecm240', 'wzp6_ee_ccH_Hcc_ecm240', 'wzp6_ee_eeH_Hgg_ecm240', 'wzp6_ee_tautauH_Hgg_ecm240', 'wzp6_ee_qqH_Hgg_ecm240', 'wzp6_ee_ssH_Hgg_ecm240', 'wzp6_ee_ccH_Hgg_ecm240', 'wzp6_ee_bbH_Hgg_ecm240', 'wzp6_ee_eeH_Hss_ecm240', 'wzp6_ee_mumuH_Hss_ecm240', 'wzp6_ee_tautauH_Hss_ecm240', 'wzp6_ee_qqH_Hss_ecm240', 'wzp6_ee_ssH_Hss_ecm240', 'wzp6_ee_ccH_Hss_ecm240', 'wzp6_ee_bbH_Hss_ecm240', 'p8_ee_WW_ecm240', 'p8_ee_ZZ_ecm240']

Zprods = ['ee', 'mumu', 'tautau', 'nunu', 'qq', 'ss', 'cc', 'bb']
bb_sig = [f'wzp6_ee_{i}H_Hbb_ecm240' for i in Zprods]
cc_sig = [f'wzp6_ee_{i}H_Hcc_ecm240' for i in Zprods]
gg_sig = [f'wzp6_ee_{i}H_Hgg_ecm240' for i in Zprods]
ss_sig = [f'wzp6_ee_{i}H_Hss_ecm240' for i in Zprods]


datasets_bkg = ["p8_ee_WW_ecm240", "p8_ee_ZZ_ecm240"]

processes = bb_sig + cc_sig + gg_sig + ss_sig + datasets_bkg

signal = [f"wzp6_ee_{cut}H_H{hdecay}{hdecay}_ecm240"]

outfile = f"datacards/{hdecay}_{cut}.txt"
#######################################################


tempfile = f"datacards/{hdecay}_{cut}_temp.txt"
f = open(tempfile, "w")

# writing i j k

f.write("imax *\njmax *\nkmax *\n")

f.write("---------------------------------------------\n") # Writing shapes input

f.write(f"shapes * cutFlow_{cut} {rootfile} $PROCESS/$CHANNEL_last\n")
f.write(f"shapes data_obs cutFlow_{cut}  {rootfile} {cut}/data_obs\n")

f.write("---------------------------------------------\n") # Writing bins and obs

f.write(f"bin         cutFlow_{cut}\n")
f.write(f"observation {totalevents}\n")

f.write("---------------------------------------------\n") # writing processes

f.write("bin                       ")
for i in range(len(processes)):
    f.write(f"cutFlow_{cut}                     ")


f.write("\nprocess              ")
for proc in processes:
    f.write(f"{proc}      ")


f.write("\nprocess                        ")
k = 1
l = 0
for proc in processes:
    if proc in signal:
        f.write(f"{l}                               ")
        
    else:
        f.write(f"{k}                               ")
        k = k + 1


f.write("\nrate                          ")
for proc in processes:
    if proc in notzero:
        f.write(f"-1                              ")
    else:
        f.write(f"0                               ")

f.write("\n---------------------------------------------\n") # writing nuisance params

f.write("lumi         lnN          ")
for proc in processes:
        f.write(f"{lumi}                         ")

f.write(f"\nacc_{hdecay}_{cut}     lnN          ")
for proc in processes:
        f.write(f"{acc}                         ")

f.write(f"\neff_{hdecay}_{cut}     lnN          ")
for proc in processes:
        f.write(f"{eff}                         ")

f.write("\n---------------------------------------------\n") # Writing autoMCstats

f.write(f"cutFlow_{cut} autoMCStats 0")

f.close()

# change below python file to your file path
os.system(f"python /home/submit/aniketkg/FCCAnalyzer/HiggsAnalysis/CombinedLimit/scripts/combineCards.py {tempfile} > {outfile}")
os.remove(tempfile)