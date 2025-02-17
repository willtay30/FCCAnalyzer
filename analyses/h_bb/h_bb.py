# import python.functions as functions      # must be commented out for running fccanalyis
import python.helpers as helpers
import ROOT
import argparse
import logging

# import helper_jetclustering
# import helper_flavourtagger

from addons.ONNXRuntime.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.jetClusteringHelper import ExclusiveJetClusteringHelper
from examples.FCCee.weaver.config import collections, njets

#***************************************************
# Notice:
# -Currently only the use of fccanalysis run [] (with proper sourcing) works
# -Much of the code is commented to allow for further testing, file may not work when uncommented
# -Code takes a long time to run to completion (>1hr)
# -Output is stored as seperate root files @ output/hbb_tagging/
#****************************************************

#modifying names 

logger = logging.getLogger("fcclogger")

# parser = functions.make_def_argparser()
# args = parser.parse_args()
# functions.set_threads(args)

# functions.add_include_file("analyses/higgs_mass_xsec/functions.h")
# functions.add_include_file("analyses/higgs_mass_xsec/functions_gen.h")

flavour = "C"  # Change to B or C or S as needed

# list of all processes
fraction = 0.05
processList = {    #Hbb sigs
    'wzp6_ee_eeH_Hbb_ecm240':          {'fraction':fraction},
    'wzp6_ee_mumuH_Hbb_ecm240':          {'fraction':fraction},
    'wzp6_ee_tautauH_Hbb_ecm240':          {'fraction':fraction},
    'wzp6_ee_nunuH_Hbb_ecm240':          {'fraction':fraction},
    'wzp6_ee_qqH_Hbb_ecm240':          {'fraction':fraction},
    'wzp6_ee_ssH_Hbb_ecm240':          {'fraction':fraction},
    'wzp6_ee_ccH_Hbb_ecm240':          {'fraction':fraction},
    'wzp6_ee_bbH_Hbb_ecm240':          {'fraction':fraction},

    #Hcc sigs
    'wzp6_ee_eeH_Hcc_ecm240':          {'fraction':fraction},
    'wzp6_ee_mumuH_Hcc_ecm240':          {'fraction':fraction},
    'wzp6_ee_tautauH_Hcc_ecm240':          {'fraction':fraction},
    'wzp6_ee_nunuH_Hcc_ecm240':          {'fraction':fraction},
    'wzp6_ee_qqH_Hcc_ecm240':          {'fraction':fraction},
    'wzp6_ee_ssH_Hcc_ecm240':          {'fraction':fraction},
    'wzp6_ee_ccH_Hcc_ecm240':          {'fraction':fraction},
    'wzp6_ee_bbH_Hcc_ecm240':          {'fraction':fraction},

    #Hgg sigs
    'wzp6_ee_eeH_Hgg_ecm240':          {'fraction':fraction},
    'wzp6_ee_mumuH_Hgg_ecm240':          {'fraction':fraction},
    'wzp6_ee_tautauH_Hgg_ecm240':          {'fraction':fraction},
    'wzp6_ee_nunuH_Hgg_ecm240':          {'fraction':fraction},
    'wzp6_ee_qqH_Hgg_ecm240':          {'fraction':fraction},
    'wzp6_ee_ssH_Hgg_ecm240':          {'fraction':fraction},
    'wzp6_ee_ccH_Hgg_ecm240':          {'fraction':fraction},
    'wzp6_ee_bbH_Hgg_ecm240':          {'fraction':fraction},

    #bkgs
    'p8_ee_WW_ecm240':          {'fraction':fraction},
    'p8_ee_ZZ_ecm240':          {'fraction':fraction},
    'wzp6_ee_mumu_ecm240':          {'fraction':fraction},
    'wzp6_ee_tautau_ecm240':          {'fraction':fraction},
    'wzp6_egamma_eZ_Zmumu_ecm240':          {'fraction':fraction},
    'wzp6_gammae_eZ_Zmumu_ecm240':          {'fraction':fraction},
    'wzp6_gaga_mumu_60_ecm240':          {'fraction':fraction},
    'wzp6_gaga_tautau_60_ecm240':          {'fraction':fraction},
    'wzp6_ee_nuenueZ_ecm240':          {'fraction':fraction},

}

inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../higgs_mass_xsec/functions.h", "../higgs_mass_xsec/functions_gen.h", "otherfuncs.h", "otherfunc_gen.h"]


# output directory
outputDir   = f"output2/h{flavour}{flavour}_tagging/histmaker/"


# define histograms

bins_m = (250, 0, 250)
bins_p = (200, 0, 200)
bins_m_zoom = (200, 110, 130) # 100 MeV


bins_theta = (500, 0, 5)
bins_phi = (400, -4, 4)

bins_count = (100, 0, 100)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)

bins_resolution = (10000, 0.95, 1.05)
bins_resolution_1 = (20000, 0, 2)

jet_energy = (1000, 0, 100) # 100 MeV bins
dijet_m = (2000, 0, 200) # 100 MeV bins
visMass = (2000, 0, 200) # 100 MeV bins
missEnergy  = (2000, 0, 200) # 100 MeV bins

dijet_m_final = (500, 50, 100) # 100 MeV bins

bins_cos = (100, -1, 1)
bins_aco = (1000,-360,360)
bins_cosThetaMiss = (10000, 0, 1)

bins_prob = (400, 0, 2)
bins_pfcand = (200, -10, 10)

# setup clustering and flavour taggingv helper_jetclustering. helper_flavourtagger.
# 2 jets

# jet2Cluster = helper_jetclustering.ExclusiveJetClusteringHelper(2, "rps_no_leps")
# jet2Flavour = helper_flavourtagger.JetFlavourHelper(collections, jet2Cluster.jets, jet2Cluster.constituents)

njets = 2
jet2Cluster = ExclusiveJetClusteringHelper("rps_no_leps", njets)
jet2Flavour = JetFlavourHelper(collections, jet2Cluster.jets, jet2Cluster.constituents, "")

# 4 jets

# jet4Cluster = helper_jetclustering.ExclusiveJetClusteringHelper(4, "ReconstructedParticles")
# jet4Flavour = helper_flavourtagger.JetFlavourHelper(collections, jet4Cluster.jets, jet4Cluster.constituents)

njets = 4
jet4Cluster = ExclusiveJetClusteringHelper("ReconstructedParticles", njets)
jet4Flavour = JetFlavourHelper(collections, jet4Cluster.jets, jet4Cluster.constituents, "")

# path = "/home/submit/aniketkg/FCCAnalyzer_ag/data/flavourtagger/fccee_flavtagging_edm4hep_wc_v1"
# urlpath = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/fccee_flavtagging_edm4hep_wc_v1"
# jet2Flavour.load(f"{path}.json", f"{path}.onnx")
# jet4Flavour.load(f"{path}.json", f"{path}.onnx")

model_name = "fccee_flavtagging_edm4hep_wc_v1"

# model files needed for unit testing in CI
url_model_dir = "https://fccsw.web.cern.ch/fccsw/testsamples/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
url_preproc = "{}/{}.json".format(url_model_dir, model_name)
url_model = "{}/{}.onnx".format(url_model_dir, model_name)

# model files locally stored on /eos
model_dir = "/eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_13_01_2022/"
local_preproc = "{}/{}.json".format(model_dir, model_name)
local_model = "{}/{}.onnx".format(model_dir, model_name)

def get_file_path(url, filename):
    import os
    if os.path.exists(filename):
        return os.path.abspath(filename)
    else:
        import urllib.request
        urllib.request.urlretrieve(url, os.path.basename(url))
        return os.path.basename(url)

weaver_preproc = get_file_path(url_preproc, local_preproc)
weaver_model = get_file_path(url_model, local_model)

def build_graph(df, dataset):

    #logging.info(f"build graph {dataset.name}")
    results, cols = [], []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")
    df = helpers.defineCutFlowVars(df) # make the cutX=X variables
    
    # define collections
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")


    # muons
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df = df.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df = df.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df = df.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(25)(muons_all)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    
    # electrons
    df = df.Alias("Electron0", "Electron#0.index")
    df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    df = df.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
    df = df.Define("electrons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons_all)")
    df = df.Define("electrons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons_all)")
    df = df.Define("electrons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons_all)")
    df = df.Define("electrons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons_all)")

    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(25)(electrons_all)")
    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df = df.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df = df.Define("electrons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons)")
    df = df.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
    df = df.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons)")

    # lepton kinematic histograms
    results.append(df.Histo1D(("muons_all_p_cut0", "", *bins_p), "muons_all_p"))
    results.append(df.Histo1D(("muons_all_theta_cut0", "", *bins_theta), "muons_all_theta"))
    results.append(df.Histo1D(("muons_all_phi_cut0", "", *bins_phi), "muons_all_phi"))
    results.append(df.Histo1D(("muons_all_q_cut0", "", *bins_charge), "muons_all_q"))
    results.append(df.Histo1D(("muons_all_no_cut0", "", *bins_count), "muons_all_no"))

    results.append(df.Histo1D(("electrons_all_p_cut0", "", *bins_p), "electrons_all_p"))
    results.append(df.Histo1D(("electrons_all_theta_cut0", "", *bins_theta), "electrons_all_theta"))
    results.append(df.Histo1D(("electrons_all_phi_cut0", "", *bins_phi), "electrons_all_phi"))
    results.append(df.Histo1D(("electrons_all_q_cut0", "", *bins_charge), "electrons_all_q"))
    results.append(df.Histo1D(("electrons_all_no_cut0", "", *bins_count), "electrons_all_no"))


    #########
    ### CUT 0: all events
    #########
    results.append(df.Histo1D(("cutFlow_mumu", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_ee", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_nunu", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_qq", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_ss", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_cc", "", *bins_count), "cut0"))
    results.append(df.Histo1D(("cutFlow_bb", "", *bins_count), "cut0"))
    

    #########
    ### CUT 1: select Z decay product
    #########
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(240., ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    results.append(df.Histo1D(("missingEnergy_nOne", "", *missEnergy), "missingEnergy"))
    
    select_mumu = "muons_no == 2 && electrons_no == 0 && missingEnergy < 30"
    select_ee   = "muons_no == 0 && electrons_no == 2 && missingEnergy < 30"
    select_nunu = "muons_no == 0 && electrons_no == 0 && missingEnergy > 102 && missingEnergy < 110"
    select_qq   = "muons_no == 0 && electrons_no == 0 && missingEnergy < 35"
    
    df_mumu   = df.Filter(select_mumu)
    df_ee     = df.Filter(select_ee)
    df_nunu   = df.Filter(select_nunu)
    df_quarks = df.Filter(select_qq)
    
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut1"))
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut1"))
    results.append(df_nunu.Histo1D(("cutFlow_nunu", "", *bins_count), "cut1"))
    results.append(df_quarks.Histo1D(("cutFlow_bb", "", *bins_count), "cut1"))
    results.append(df_quarks.Histo1D(("cutFlow_cc", "", *bins_count), "cut1"))
    results.append(df_quarks.Histo1D(("cutFlow_ss", "", *bins_count), "cut1"))
    results.append(df_quarks.Histo1D(("cutFlow_qq", "", *bins_count), "cut1"))


    #########
    ### CUT 2: we want to detect Z->mumu/ee (so we don't have Z->qq interfering with measurement of H->bb)
    ###        Z->qq and Z->nunu do not get the next few cuts
    #########

    # build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system (Z), index and 2 the leptons of the pair
    
    # muons
    df_mumu = df_mumu.Define("hbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0, 240, false)(muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df_mumu = df_mumu.Filter("hbuilder_result.size() > 0")
    
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut2"))

    df_mumu = df_mumu.Define("zmumu", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[0]}") # the Z
    df_mumu = df_mumu.Define("zmumu_tlv", "FCCAnalyses::makeLorentzVectors(zmumu)") # the muons
    df_mumu = df_mumu.Define("zmumu_leps", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[1],hbuilder_result[2]}")
    df_mumu = df_mumu.Define("zmumu_leps_tlv", "FCCAnalyses::makeLorentzVectors(zmumu_leps)")

    df_mumu = df_mumu.Define("zmumu_m", "FCCAnalyses::ReconstructedParticle::get_mass(zmumu)[0]")
    df_mumu = df_mumu.Define("zmumu_p", "FCCAnalyses::ReconstructedParticle::get_p(zmumu)[0]")
    df_mumu = df_mumu.Define("zmumu_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(zmumu)")
    df_mumu = df_mumu.Define("zmumu_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(zmumu_recoil)[0]")

    # electrons
    df_ee = df_ee.Define("hbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0, 240, false)(electrons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df_ee = df_ee.Filter("hbuilder_result.size() > 0")
    
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut2"))

    df_ee = df_ee.Define("zee", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[0]}") # the Z
    df_ee = df_ee.Define("zee_tlv", "FCCAnalyses::makeLorentzVectors(zee)") # the electrons
    df_ee = df_ee.Define("zee_leps", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{hbuilder_result[1],hbuilder_result[2]}")
    df_ee = df_ee.Define("zee_leps_tlv", "FCCAnalyses::makeLorentzVectors(zee_leps)")

    df_ee = df_ee.Define("zee_m", "FCCAnalyses::ReconstructedParticle::get_mass(zee)[0]")
    df_ee = df_ee.Define("zee_p", "FCCAnalyses::ReconstructedParticle::get_p(zee)[0]")
    df_ee = df_ee.Define("zee_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(240)(zee)")
    df_ee = df_ee.Define("zee_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(zee_recoil)[0]")



    #########
    ### CUT 3: recoil cut (H mass)
    #########
    results.append(df_mumu.Histo1D(("mumu_recoil_m_nOne", "", *bins_m), "zmumu_recoil_m"))
    #make a histogram that plots what's in zmumu_recoil. make another histogram right after the filter
    #do that for all cuts
    #look thru histograms for each cut
    #make slides: signal before and after each cut (ZH->mumucc)
    df_mumu = df_mumu.Filter("zmumu_recoil_m > 123 && zmumu_recoil_m < 132")
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut3"))
    
    results.append(df_ee.Histo1D(("ee_recoil_m_nOne", "", *bins_m), "zee_recoil_m"))
    df_ee = df_ee.Filter("zee_recoil_m > 123 && zee_recoil_m < 132")
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut3"))
    
    # graphs for inspecting the WW background
    #results.append(df_mumu.Graph("missingEnergy", "zmumu_recoil_m"))
    #results.append(df_mumu.Graph("zmumu_recoil_m", "missingEnergy"))
    

    #########
    ### CUT 4: momentum
    #########
    results.append(df_mumu.Histo1D(("mumu_p_nOne", "", *bins_p), "zmumu_p"))
    df_mumu = df_mumu.Filter("zmumu_p > 45 && zmumu_p < 55")
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut4"))

    results.append(df_ee.Histo1D(("ee_p_nOne", "", *bins_p), "zee_p"))
    df_ee = df_ee.Filter("zee_p > 45 && zee_p < 55")
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut4"))


    #########
    ### CUT 5: cut on Z mass
    #########
    results.append(df_mumu.Histo1D(("zmumu_m_nOne", "", *bins_m), "zmumu_m"))
    df_mumu = df_mumu.Filter("zmumu_m > 85 && zmumu_m < 95")
    results.append(df_mumu.Histo1D(("cutFlow_mumu", "", *bins_count), "cut5"))

    results.append(df_ee.Histo1D(("zee_m_nOne", "", *bins_m), "zee_m"))
    df_ee = df_ee.Filter("zee_m > 85 && zee_m < 95")
    results.append(df_ee.Histo1D(("cutFlow_ee", "", *bins_count), "cut5"))


    for leps, df in [("muons", df_mumu), ("electrons", df_ee), ("neutrinos", df_nunu)]:
        # define PF candidates collection by removing the leptons
        if leps != "neutrinos":
            df = df.Define("rps_no_leps", f"FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, {leps})")
        else:
            df = df.Alias("rps_no_leps", "ReconstructedParticles")
        
        # clustering
        df = jet2Cluster.define(df)
        #df = df.Define("jet_tlv", "FCCAnalyses::makeLorentzVectors(jet_px, jet_py, jet_pz, jet_e)")
        df = df.Define("jet_tlv", f"JetConstituentsUtils::compute_tlv_jets({jet2Cluster.jets})")

        # calculate dijet m and p
        df = df.Define("dijet", "jet_tlv[0] + jet_tlv[1]")
        df = df.Define("dijet_m", "dijet.M()")
        df = df.Define("dijet_p", "dijet.P()")
        
        # for neutrinos, cut on Higgs mass and momentum
        if leps == "neutrinos":
            results.append(df.Histo1D((f"z{leps}_h_m_nOne", "", *bins_m), "dijet_m"))
            results.append(df.Histo1D((f"z{leps}_h_p_nOne", "", *bins_m), "dijet_p"))
            
            df = df.Filter("dijet_p > 40 && dijet_p < 58")
            results.append(df.Histo1D(("cutFlow_nunu", "", *bins_count), "cut2"))
            
            df = df.Filter("dijet_m > 115 && dijet_m < 128")
            results.append(df.Histo1D(("cutFlow_nunu", "", *bins_count), "cut3"))
        
        # flavour tagging

        #df = jet2Flavour.define_and_inference(df)
        df = jet2Flavour.define(df) # define variables
        df = jet2Flavour.inference(weaver_preproc, weaver_model, df) # run inference
        
        # cut on jet confidence
        #LOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOK HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
        df = df.Filter(f"recojet_is{flavour}[0] > 0.5 && recojet_is{flavour}[1] > 0.5")
        
        if leps != "neutrinos":
            results.append(df.Histo1D((f"cutFlow_{'mumu' if leps == 'muons' else 'ee'}", "", *bins_count), "cut6"))
            
            # store the final recoil mass of ee and mumu for a fit
            results.append(df.Histo1D((f"z{leps}_final_recoil_m", "", *bins_m), f"z{'mumu' if leps == 'muons' else 'ee'}_recoil_m"))
        else:
            results.append(df.Histo1D(("cutFlow_nunu", "", *bins_count), "cut4"))
        
        # store final dijet mass and momentum
        results.append(df.Histo1D((f"z{leps}_h_m", "", *bins_m), "dijet_m"))
        results.append(df.Histo1D((f"z{leps}_h_p", "", *bins_m), "dijet_p"))



    # jet analysis for the case of 2 jets (Z -> leps) --------------



    # Z->qq analyses
    # clustering
    df_quarks = jet4Cluster.define(df_quarks)
    #df_quarks = df_quarks.Define("jet_tlv", "FCCAnalyses::makeLorentzVectors(jet_px, jet_py, jet_pz, jet_e)")
    df_quarks = df_quarks.Define("jet_tlv", f"JetConstituentsUtils::compute_tlv_jets({jet4Cluster.jets})")

    # pair jets based on distance to Z and H masses
    df_quarks = df_quarks.Define("zh_min_idx", """
            FCCAnalyses::Vec_i min{0, 0, 0, 0};
            float distm = INFINITY;
            for (int i = 0; i < 3; i++)
                for (int j = i + 1; j < 4; j++)
                    for (int k = 0; k < 3; k++) {
                        if (i == k || j == k) continue;
                        for (int l = k + 1; l < 4; l++) {
                            if (i == l || j == l) continue;
                            float distz = (jet_tlv[i] + jet_tlv[j]).M() - 91.2;
                            float disth = (jet_tlv[k] + jet_tlv[l]).M() - 125;
                            if (distz*distz/91.2 + disth*disth/125 < distm) {
                                distm = distz*distz/91.2 + disth*disth/125;
                                min[0] = i; min[1] = j; min[2] = k; min[3] = l;
                            }
                        }
                    }
            return min;""")

    # compute Z and H masses and momenta
    df_quarks = df_quarks.Define("z_dijet", "jet_tlv[zh_min_idx[0]] + jet_tlv[zh_min_idx[1]]")
    df_quarks = df_quarks.Define("h_dijet", "jet_tlv[zh_min_idx[2]] + jet_tlv[zh_min_idx[3]]")

    df_quarks = df_quarks.Define("z_dijet_m", "z_dijet.M()")
    df_quarks = df_quarks.Define("z_dijet_p", "z_dijet.P()")
    df_quarks = df_quarks.Define("h_dijet_m", "h_dijet.M()")
    df_quarks = df_quarks.Define("h_dijet_p", "h_dijet.P()")

    results.append(df_quarks.Histo1D(("quarks_z_m_nOne", "", *bins_m), "z_dijet_m"))
    results.append(df_quarks.Histo1D(("quarks_z_p_nOne", "", *bins_m), "z_dijet_p"))
    results.append(df_quarks.Histo1D(("quarks_h_m_nOne", "", *bins_m), "h_dijet_m"))
    results.append(df_quarks.Histo1D(("quarks_h_p_nOne", "", *bins_m), "h_dijet_p"))

    # filter on Z momentum
    df_quarks = df_quarks.Filter("z_dijet_p > 45 && z_dijet_p < 56")
    results.append(df_quarks.Histo1D(("cutFlow_bb", "", *bins_count), "cut2"))
    results.append(df_quarks.Histo1D(("cutFlow_cc", "", *bins_count), "cut2"))
    results.append(df_quarks.Histo1D(("cutFlow_ss", "", *bins_count), "cut2"))
    results.append(df_quarks.Histo1D(("cutFlow_qq", "", *bins_count), "cut2"))

    # filter on H mass
    df_quarks = df_quarks.Filter("h_dijet_m > 122 && h_dijet_m < 128")
    results.append(df_quarks.Histo1D(("cutFlow_bb", "", *bins_count), "cut3"))
    results.append(df_quarks.Histo1D(("cutFlow_cc", "", *bins_count), "cut3"))
    results.append(df_quarks.Histo1D(("cutFlow_ss", "", *bins_count), "cut3"))
    results.append(df_quarks.Histo1D(("cutFlow_qq", "", *bins_count), "cut3"))


    # flavour tagging
    #df_quarks = jet4Flavour.define_and_inference(df_quarks)
    df_quarks = jet4Flavour.define(df_quarks) # define variables
    df_quarks = jet4Flavour.inference(weaver_preproc, weaver_model, df_quarks) # run inference


    # get tag confidence
    df_quarks = df_quarks.Define(f"H{flavour}{flavour}_prob", f"std::min(recojet_is{flavour}[zh_min_idx[2]], recojet_is{flavour}[zh_min_idx[3]])")
    results.append(df_quarks.Histo1D((f"H{flavour}{flavour}_prob_nOne", "", *bins_prob), f"H{flavour}{flavour}_prob"))

    df_quarks = df_quarks.Define("Zbb_prob", "std::min(recojet_isB[zh_min_idx[0]], recojet_isB[zh_min_idx[1]])")
    df_quarks = df_quarks.Define("Zcc_prob", "std::min(recojet_isC[zh_min_idx[0]], recojet_isC[zh_min_idx[1]])")
    df_quarks = df_quarks.Define("Zss_prob", "std::min(recojet_isS[zh_min_idx[0]], recojet_isS[zh_min_idx[1]])")
    df_quarks = df_quarks.Define("Zqq_prob", "std::min(recojet_isQ[zh_min_idx[0]], recojet_isQ[zh_min_idx[1]])")

    results.append(df_quarks.Histo1D(("Zbb_prob_nOne", "", *bins_prob), "Zbb_prob"))
    results.append(df_quarks.Histo1D(("Zcc_prob_nOne", "", *bins_prob), "Zcc_prob"))
    results.append(df_quarks.Histo1D(("Zss_prob_nOne", "", *bins_prob), "Zss_prob"))
    results.append(df_quarks.Histo1D(("Zqq_prob_nOne", "", *bins_prob), "Zqq_prob"))

    # sort by most likely tag
    df_quarks = df_quarks.Define("Zbb_like", "recojet_isB[zh_min_idx[0]] + recojet_isB[zh_min_idx[1]]")
    df_quarks = df_quarks.Define("Zcc_like", "recojet_isC[zh_min_idx[0]] + recojet_isC[zh_min_idx[1]]")
    df_quarks = df_quarks.Define("Zss_like", "recojet_isS[zh_min_idx[0]] + recojet_isS[zh_min_idx[1]]")
    df_quarks = df_quarks.Define("Zqq_like", "recojet_isQ[zh_min_idx[0]] + recojet_isQ[zh_min_idx[1]]")

    df_quarks = df_quarks.Define("best_tag", """
            if (Zbb_like > Zcc_like && Zbb_like > Zss_like && Zbb_like > Zqq_like) {
                return 0;
            } else if (Zcc_like > Zss_like && Zcc_like > Zqq_like) {
                return 1;
            } else if (Zss_like > Zqq_like) {
                return 2;
            } else {
                return 3;
            } """)

    # sort by maximum likelihood tag
    df_bb = df_quarks.Filter("best_tag == 0")
    df_cc = df_quarks.Filter("best_tag == 1")
    df_ss = df_quarks.Filter("best_tag == 2")
    df_qq = df_quarks.Filter("best_tag == 3")

    results.append(df_bb.Histo1D(("cutFlow_bb", "", *bins_count), "cut4"))
    results.append(df_cc.Histo1D(("cutFlow_cc", "", *bins_count), "cut4"))
    results.append(df_ss.Histo1D(("cutFlow_ss", "", *bins_count), "cut4"))
    results.append(df_qq.Histo1D(("cutFlow_qq", "", *bins_count), "cut4"))

    results.append(df_bb.Graph(f"H{flavour}{flavour}_prob", "Zbb_prob"))
    results.append(df_cc.Graph(f"H{flavour}{flavour}_prob", "Zcc_prob"))
    results.append(df_ss.Graph(f"H{flavour}{flavour}_prob", "Zss_prob"))
    results.append(df_qq.Graph(f"H{flavour}{flavour}_prob", "Zqq_prob"))

    # make sure there are two b (or c or s) jets
    # switch case based on flavour so that you don't have to change the probabilities each time you change flavour
    # 16 possibilities
    # probabilities obtained thru optimization
    df_qq = df_qq.Filter(f"H{flavour}{flavour}_prob > 0.032") #bb: .032
    df_ss = df_ss.Filter(f"H{flavour}{flavour}_prob > 0.032") #bb: .032
    df_cc = df_cc.Filter(f"H{flavour}{flavour}_prob > 0.029") #bb: .029
    df_bb = df_bb.Filter(f"H{flavour}{flavour}_prob > 0.011") #bb: .011

    results.append(df_bb.Histo1D(("cutFlow_bb", "", *bins_count), "cut5"))
    results.append(df_cc.Histo1D(("cutFlow_cc", "", *bins_count), "cut5"))
    results.append(df_ss.Histo1D(("cutFlow_ss", "", *bins_count), "cut5"))
    results.append(df_qq.Histo1D(("cutFlow_qq", "", *bins_count), "cut5"))

    # check that the Z jets are the right type
    df_bb = df_bb.Filter("Zbb_prob > 0.042")
    df_cc = df_cc.Filter("Zcc_prob > 0.134")
    df_ss = df_ss.Filter("Zss_prob > 0.095")
    df_qq = df_qq.Filter("Zqq_prob > 0.053")

    results.append(df_bb.Histo1D(("cutFlow_bb", "", *bins_count), "cut6"))
    results.append(df_cc.Histo1D(("cutFlow_cc", "", *bins_count), "cut6"))
    results.append(df_ss.Histo1D(("cutFlow_ss", "", *bins_count), "cut6"))
    results.append(df_qq.Histo1D(("cutFlow_qq", "", *bins_count), "cut6"))

    # make final mass and momentum histograms
    for q, df in [("bb", df_bb), ("cc", df_cc), ("ss", df_ss), ("qq", df_qq)]:
        results.append(df.Histo1D((f"z{q}_z_m", "", *bins_m), "z_dijet_m"))
        results.append(df.Histo1D((f"z{q}_h_m", "", *bins_m), "h_dijet_m"))
        results.append(df.Histo1D((f"z{q}_z_p", "", *bins_m), "z_dijet_p"))
        results.append(df.Histo1D((f"z{q}_h_p", "", *bins_m), "h_dijet_p"))


    
    return results, weightsum


# Main function must be uncommented when running file using python, else keep commented

# if __name__ == "__main__":

#     datadict = functions.get_datadicts() # get default datasets

#     Zprods = ["ee", "mumu", "tautau", "nunu", "qq", "ss", "cc", "bb"] 
#     bb_sig = [f"wzp6_ee_{i}H_Hbb_ecm240" for i in Zprods]
#     cc_sig = [f"wzp6_ee_{i}H_Hcc_ecm240" for i in Zprods]
#     gg_sig = [f"wzp6_ee_{i}H_Hgg_ecm240" for i in Zprods]
    
#     quark_test = ["wzp6_ee_qqH_Hbb_ecm240", "wzp6_ee_ssH_Hbb_ecm240", "wzp6_ee_ccH_Hbb_ecm240", "wzp6_ee_bbH_Hbb_ecm240"]
    
#     datasets_bkg = ["p8_ee_WW_ecm240", "p8_ee_ZZ_ecm240", "wzp6_ee_mumu_ecm240", "wzp6_ee_tautau_ecm240", "wzp6_egamma_eZ_Zmumu_ecm240", "wzp6_gammae_eZ_Zmumu_ecm240", "wzp6_gaga_mumu_60_ecm240", "wzp6_gaga_tautau_60_ecm240", "wzp6_ee_nuenueZ_ecm240"]

#     datasets_to_run = bb_sig + cc_sig + gg_sig + datasets_bkg[:2]
    
#     result = functions.build_and_run(datadict, datasets_to_run, build_graph, f"end_of_h_bb.root", args, norm=True, lumi=7200000)
