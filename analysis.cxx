#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

static std::vector<std::string> myVarSet = {"run", "lumi", "evt",
                                            "npv", "rho",

                                            "tau_pt", "tau_eta", "tau_phi", "tau_mass",
                                            "tau_charge", "tau_decayMode",
                                            "tau_byDeepTau2017v2p1VSjetraw", "tau_byDeepTau2017v2p1VSjet",
                                            "tau_byDeepTau2017v2p1VSeraw", "tau_byDeepTau2017v2p1VSe",
                                            "tau_byDeepTau2017v2p1VSmuraw", "tau_byDeepTau2017v2p1VSmu",

                                            "jet_pt", "jet_eta", "jet_phi", "jet_mass",
                                            "jet_partonFlavour", "jet_hadronFlavour",

                                            "ele_pt", "ele_eta", "ele_phi", "ele_mass",

                                            "muon_pt", "muon_eta", "muon_phi", "muon_mass", "muon_type",

                                            "pfCand_pt", "pfCand_eta", "pfCand_phi", "pfCand_mass",
                                            "pfCand_charge", "pfCand_tauIso", "pfCand_tauSignal", "pfCand_particleType", "pfCand_fromPV",
                                            "pfCand_puppiWeight", "pfCand_puppiWeightNoLep"
};

namespace basefunctions {

void makeSnapshot(auto df, const std::string treename, const std::string outputname, const std::vector<std::string> varSet) {

    df.template Snapshot<unsigned int, unsigned int, unsigned long long, // "run", "lumi", "evt",
                int, float, // "npv", "rho",
                float, float, float, float, // "tau_pt", "tau_eta", "tau_phi", "tau_mass",
                int, int, // "tau_charge", "tau_decayMode",
                float, unsigned short, // "tau_byDeepTau2017v2p1VSjetraw", "tau_byDeepTau2017v2p1VSjet",
                float, unsigned short, // "tau_byDeepTau2017v2p1VSeraw", "tau_byDeepTau2017v2p1VSe",
                float, unsigned short, // "tau_byDeepTau2017v2p1VSmuraw", "tau_byDeepTau2017v2p1VSmu",
                float, float, float, float, // "jet_pt", "jet_eta", "jet_phi", "jet_mass",
                int, int, // "jet_partonFlavour", "jet_hadronFlavour",
                std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, // "ele_pt", "ele_eta", "ele_phi", "ele_mass",
                std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<unsigned int>, // "muon_pt", "muon_eta", "muon_phi", "muon_mass", "muon_type"
                std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, // "pfCand_pt", "pfCand_eta", "pfCand_phi", "pfCand_mass",
                std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, // "pfCand_charge", "pfCand_tauIso", "pfCand_tauSignal", "pfCand_particleType", "pfCand_fromPV"
                std::vector<float>, std::vector<float>, // "pfCand_puppiWeight", "pfCand_puppiWeightNoLep"
                int, // gen_match
                ROOT::VecOps::RVec<float>, // ele_hcalOverEcal
                ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>, // ele_deltaR, ele_deltaPhi, ele_deltaEta,
                ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>, // muon_deltaR, muon_deltaPhi, muon_deltaEta,
                ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float> // pfCand_deltaR, pfCand_deltaPhi, pfCand_deltaEta,
                >(treename, outputname, varSet);
}

auto GetGenMatch() {
    return [](const float &tau_pt, const float &tau_eta, const float &tau_phi,
              const float &gen_pt, const float &gen_eta, const float &gen_phi,
              const int &gen_index, const int &gen_kind, const int &jet_index){
        auto deltaR = ROOT::VecOps::DeltaR(tau_eta, gen_eta, tau_phi, gen_phi);
        auto deltaR_match = deltaR <= 0.2;
        auto gen_lepton_valid = gen_index >= 0;
        auto gen_jet_valid = jet_index >= 0;

        auto gen_jet = (!gen_lepton_valid && gen_jet_valid);

        auto good_gen_lepton = deltaR_match && gen_lepton_valid && (gen_pt >= 8.0);

        auto gen_prompt_electron = good_gen_lepton && (gen_kind == 1);
        auto gen_prompt_muon = good_gen_lepton && (gen_kind == 2);

        auto gen_tau_electron = good_gen_lepton && (gen_kind == 3);
        auto gen_tau_muon = good_gen_lepton && (gen_kind == 4);

        auto gen_had_tau = good_gen_lepton && (gen_pt >= 15.0) && (gen_kind == 5);

        auto gen_match = 1 * gen_prompt_electron +
                         2 * gen_prompt_muon +
                         3 * gen_tau_electron +
                         4 * gen_tau_muon +
                         5 * gen_had_tau +
                         6 * gen_jet;
        return gen_match;
    };
}

auto GetDeltaR(){
    return [](const float &tau_eta, const float &tau_phi,
              const ROOT::RVec<float> &cands_eta, const ROOT::RVec<float> &cands_phi){
        ROOT::RVec<float> taus_eta;
        ROOT::RVec<float> taus_phi;
        for(unsigned int i =0; i < cands_eta.size(); i++){
            taus_eta.push_back(tau_eta);
            taus_phi.push_back(tau_phi);
        }
        return ROOT::VecOps::DeltaR(taus_eta, cands_eta, taus_phi, cands_phi);
    };
}

auto GetDeltaPhi(){
    return [](const float &tau_phi,
              const ROOT::RVec<float> &cands_phi){
        ROOT::RVec<float> taus_phi;
        for(unsigned int i =0; i < cands_phi.size(); i++){
            taus_phi.push_back(tau_phi);
        }
        return ROOT::VecOps::DeltaPhi(taus_phi, cands_phi);
    };
}

auto GetDeltaEta(){
    return [](const float &tau_eta,
              const ROOT::RVec<float> &cands_eta){
        ROOT::RVec<float> deltaEta;
        for(unsigned int i =0; i < cands_eta.size(); i++){
            deltaEta.push_back(std::abs(tau_eta - cands_eta.at(i)));
        }
        return deltaEta;
    };
}

} // end namespace basefunctions

int main(int argc, char* argv[]) {

    if (argc != 2)
    {
        std::cout << "Please provide the name of the sample as argument\n";
        std::cout << "Usage: ./analysis <sample-name>\n";
        exit(1);
    }

    std::string sample_name(argv[1]);
    std::cout << "Using files for sample: " << sample_name << std::endl;
    std::string folder = "/ceph/akhmet/forAndrewIsaac/prod_2018_v2/" + sample_name + "/";
    std::string out = "/ceph/akhmet/forAndrewIsaac/prod_2018_v2_processed_v3/";
    //std::string out = "./";

    gErrorIgnoreLevel = kError;
    ROOT::EnableImplicitMT(20);

    ROOT::RDataFrame df(
        "taus",
        folder+"*.root");
        //{folder+"eventTuple_1.root", folder+"eventTuple_2.root"});

    auto df_taupt_gteq20 = df.Filter([](const ROOT::RVec<float> &pt){return Any(pt >= 20.0);}, {"tau_pt"});
    auto df_taudm_good = df_taupt_gteq20.Filter([](const ROOT::RVec<int> &dm){return Any(dm < 5 || dm > 9);}, {"tau_decayMode"});

    auto df_gen_match = df_taudm_good.Define("gen_match", basefunctions::GetGenMatch(),
                                  {"tau_pt", "tau_eta", "tau_phi",
                                   "genLepton_vis_pt", "genLepton_vis_eta", "genLepton_vis_phi",
                                   "genLepton_index", "genLepton_kind", "genJet_index"});
    myVarSet.push_back("gen_match");

    auto df_hOe = df_gen_match.Define("ele_hcalOverEcal", [](const ROOT::RVec<float> hOe1, const ROOT::RVec<float> hOe2){return hOe1 + hOe2;}, {"ele_hcalDepth1OverEcal", "ele_hcalDepth2OverEcal"});
    myVarSet.push_back("ele_hcalOverEcal");

    auto df_ele_deltaR = df_hOe.Define("ele_deltaR",basefunctions::GetDeltaR(), {"tau_eta", "tau_phi", "ele_eta", "ele_phi"});
    myVarSet.push_back("ele_deltaR");

    auto df_ele_deltaPhi = df_ele_deltaR.Define("ele_deltaPhi",basefunctions::GetDeltaPhi(), {"tau_phi", "ele_phi"});
    myVarSet.push_back("ele_deltaPhi");

    auto df_ele_deltaEta = df_ele_deltaPhi.Define("ele_deltaEta",basefunctions::GetDeltaEta(), {"tau_eta", "ele_eta"});
    myVarSet.push_back("ele_deltaEta");

    auto df_muon_deltaR = df_ele_deltaEta.Define("muon_deltaR",basefunctions::GetDeltaR(), {"tau_eta", "tau_phi", "muon_eta", "muon_phi"});
    myVarSet.push_back("muon_deltaR");

    auto df_muon_deltaPhi = df_muon_deltaR.Define("muon_deltaPhi",basefunctions::GetDeltaPhi(), {"tau_phi", "muon_phi"});
    myVarSet.push_back("muon_deltaPhi");

    auto df_muon_deltaEta = df_muon_deltaPhi.Define("muon_deltaEta",basefunctions::GetDeltaEta(), {"tau_eta", "muon_eta"});
    myVarSet.push_back("muon_deltaEta");

    auto df_pfCand_deltaR = df_muon_deltaEta.Define("pfCand_deltaR",basefunctions::GetDeltaR(), {"tau_eta", "tau_phi", "pfCand_eta", "pfCand_phi"});
    myVarSet.push_back("pfCand_deltaR");

    auto df_pfCand_deltaPhi = df_pfCand_deltaR.Define("pfCand_deltaPhi",basefunctions::GetDeltaPhi(), {"tau_phi", "pfCand_phi"});
    myVarSet.push_back("pfCand_deltaPhi");

    auto df_pfCand_deltaEta = df_pfCand_deltaPhi.Define("pfCand_deltaEta",basefunctions::GetDeltaEta(), {"tau_eta", "pfCand_eta"});
    myVarSet.push_back("pfCand_deltaEta");

    auto df_final = df_pfCand_deltaEta;

    auto df_electrons = df_final.Filter([](const int &match){return (match == 1 || match == 3);}, {"gen_match"});
    auto df_muons = df_final.Filter([](const int &match){return (match == 2 || match == 4);}, {"gen_match"});
    auto df_taus = df_final.Filter([](const int &match){return (match == 5);}, {"gen_match"});
    auto df_jets = df_final.Filter([](const int &match){return (match == 6);}, {"gen_match"});
    auto df_invalid = df_final.Filter([](const int &match){return (match == 0);}, {"gen_match"});

    basefunctions::makeSnapshot(df_electrons, "taus", out + sample_name + "_electrons.root", myVarSet);
    basefunctions::makeSnapshot(df_muons, "taus", out + sample_name + "_muons.root", myVarSet);
    basefunctions::makeSnapshot(df_taus, "taus", out + sample_name + "_taus.root", myVarSet);
    basefunctions::makeSnapshot(df_jets, "taus", out + sample_name + "_jets.root", myVarSet);
    basefunctions::makeSnapshot(df_invalid, "taus", out + sample_name + "_invalid.root", myVarSet);
}
