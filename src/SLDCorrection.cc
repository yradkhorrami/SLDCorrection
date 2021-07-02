#include "SLDCorrection.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include "EVENT/MCParticle.h"
#include "EVENT/Vertex.h"
#include "EVENT/ReconstructedParticle.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "IMPL/ParticleIDImpl.h"
#include "UTIL/PIDHandler.h"
#include "marlin/VerbosityLevels.h"
#include <GeometryUtil.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TF1.h"
#include "TTree.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TRatioPlot.h"
#include "TAxis.h"
#include "TLine.h"

using namespace lcio ;
using namespace marlin ;

SLDCorrection aSLDCorrection;

SLDCorrection::SLDCorrection() :

	Processor("SLDCorrection"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	m_Bfield(0.f),
	c(0.),
	mm2m(0.),
	eV2GeV(0.),
	eB(0.),
	foundFlightDirection(true),
	m_nTauSLDecay(0),
	m_nTauNeutrino(0),
	m_nNeutrino(0),
	m_nChargedPFOwoTrack(0),
	n_NuPxResidual(0),
	n_NuPyResidual(0),
	n_NuPzResidual(0),
	n_NuEResidual(0),
	n_NuPxNormalizedResidual(0),
	n_NuPyNormalizedResidual(0),
	n_NuPzNormalizedResidual(0),
	n_NuENormalizedResidual(0)
{
	_description = "SLDCorrection finds semi-leptonic decays within jets and performs a correction to 4-momentum of the jet due to the missing neutrino(s)";

	registerInputCollection(	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					m_mcParticleCollection,
					std::string("MCParticle")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"PfoCollection",
					"Name of input pfo collection",
					m_inputPfoCollection,
					std::string("PandoraPFOs")
				);

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"JetCollection",
					"Name of input jet collection",
					m_inputJetCollection,
					std::string("Durham_nJets")
				);

	registerInputCollection(	LCIO::MCPARTICLE,
					"jetFlavour" ,
					"Name of input Jet Flavour collection",
					m_jetFlavour ,
					std::string("jetFlavour")
				);

	registerInputCollection(	LCIO::VERTEX,
					"PrimaryVertex",
					"Name of Primary Vertex Collection",
					m_inputPrimaryVertex,
					std::string("PrimaryVertex")
				);

	registerInputCollection(	LCIO::VERTEX,
					"BuildUpVertex",
					"Name of BuildUp Vertex Collection",
					m_inputBuildUpVertex,
					std::string("BuildUpVertex")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"RecoMCTruthLinkCollection",
					"Name of input RecoMCTruthLink Collection",
					m_RecoMCTruthLinkCollection,
					std::string("RecoMCTruthLinkCollection")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthRecoLinkCollection",
					"Name of input MCTruthRecoLink Collection",
					m_MCTruthRecoLinkCollection,
					std::string("MCTruthRecoLinkCollection")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"TrackMCTruthLinkCollection",
					"Name of input TrackMCTruthLink Collection",
					m_TrackMCTruthLinkCollection,
					std::string("MarlinTrkTracksMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthTrackLinkCollection",
					"Name of input MCTruthTrackLink Collection",
					m_MCTruthTrackLinkCollection,
					std::string("MCTruthMarlinTrkTracksLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"ClusterMCTruthLinkCollection",
					"Name of input m_ClusterMCTruthLink Collection",
					m_ClusterMCTruthLinkCollection,
					std::string("ClusterMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthClusterLinkCollection",
					"Name of input MCTruthClusterLink Collection",
					m_MCTruthClusterLinkCollection,
					std::string("MCTruthClusterLink")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"SLDNeutrinoCollection",
					"Name of semi-leptonic Neutrino collection",
					m_SLDNuCollection,
					std::string("SLDNeutrinoCollection")
				);

	registerProcessorParameter(	"includbJets",
					"Include b-jets for semi-decay correction",
					m_includbJets,
					bool(false)
				);

	registerProcessorParameter(	"includcJets",
					"Include c-jets for semi-decay correction",
					m_includcJets,
					bool(false)
				);

	registerProcessorParameter(	"includgJets",
					"Include g-jets for semi-decay correction",
					m_includgJets,
					bool(false)
				);

	registerProcessorParameter(	"includOthers",
					"Include Other final states for semi-decay correction",
					m_includOthers,
					bool(false)
				);

	registerProcessorParameter(	"includeBSLD",
					"do correction for semi-leptonic decays of B-Hadrons",
					m_includeBSLD,
					bool(true)
				);

	registerProcessorParameter(	"includeCSLD",
					"do correction for semi-leptonic decays of C-Hadrons",
					m_includeCSLD,
					bool(true)
				);

	registerProcessorParameter(	"includeTSLD",
					"do correction for semi-leptonic decays of Tau-Leptons",
					m_includeTSLD,
					bool(true)
				);

	registerProcessorParameter(	"cheatSLDLeptons",
					"Cheat semi-leptonic decays lepton from MCTruth",
					m_cheatSLDLeptons,
					bool(true)
				);

	registerProcessorParameter(	"cheatFlightDirection",
					"Cheat Flight direction of mother hadron",
					m_cheatFlightDirection,
					bool(true)
				);

	registerProcessorParameter(	"useJetAxisAsFlightDirection",
					"use jet axis as flight direction of mother hadron",
					m_useJetAxisAsFlightDirection,
					bool(true)
				);

	registerProcessorParameter(	"considerParentCharge",
					"Consider charge of parent hadron for calculating flight direction",
					m_considerParentCharge,
					bool(true)
				);

	registerProcessorParameter(	"cheatVertices",
					"Cheat vertices of mother hadron for calculating flight direction",
					m_cheatVertices,
					bool(true)
				);

	registerProcessorParameter(	"cheatLepton4momentum",
					"Cheat FourMomentum of lepton in semi-leptonic decays",
					m_cheatLepton4momentum,
					bool(true)
				);

	registerProcessorParameter(	"cheatCharged4momentum",
					"Cheat FourMomentum of charged visibles in semi-leptonic decays",
					m_cheatCharged4momentum,
					bool(true)
				);

	registerProcessorParameter(	"cheatNeutral4momentum",
					"Cheat FourMomentum of neutral visibles in semi-leptonic decays",
					m_cheatNeutral4momentum,
					bool(true)
				);

	registerProcessorParameter(	"nIterFlightDirCorrection",
					"Number of iterations for correcting flight direction of CHARGED parent hadron",
					m_nIterFlightDirCorrection,
					int(1)
				);

	registerProcessorParameter(	"recoFourMomentumOfVisibles",
					"0: get 4p from linked track/cluster to MCP, 1: get 4p from PFO with linked track/cluster, 2: get 4p from linked PFO",
					m_recoFourMomentumOfVisibles,
					int(0)
				);

	registerProcessorParameter(	"fillRootTree",
					"Fill root tree to check processor performance",
					m_fillRootTree,
					bool(true)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("Output.root")
				);

}

void SLDCorrection::init()
{
	streamlog_out(DEBUG) << "	init called  " << std::endl;
	m_Bfield = 3.5;
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;
	printParameters();
	if ( m_fillRootTree )
	{
		m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
		m_pTTree = new TTree("SLDCorrection", "SLDCorrection");
		m_pTTree->SetDirectory(m_pTFile);
		m_pTTree->Branch("event", &m_nEvt, "event/I");
		m_pTTree->Branch("nTauSLDecay",&m_nTauSLDecay,"nTauSLDecay/I");
		m_pTTree->Branch("nTauNeutrino",&m_nTauNeutrino,"nTauNeutrino/I");
		m_pTTree->Branch("nNeutrino",&m_nNeutrino,"nNeutrino/I");
		m_pTTree->Branch("jetFlavourPDG",&m_jetFlavourPDG);
		m_pTTree->Branch("nSLD_chargedMCPwoTrack",&m_nSLD_chargedMCPwoTrack);
		m_pTTree->Branch("GenStatParentHadron",&m_GenStatParentHadron);
		m_pTTree->Branch("ChargeParentHadron",&m_ChargeParentHadron);
		m_pTTree->Branch("foundRecoLepton",&m_foundRecoLepton);
		m_pTTree->Branch("foundBuildUpVertex",&m_foundBuildUpVertex);
		m_pTTree->Branch("foundRecoLeptonInBuildUpVertex",&m_foundRecoLeptonInBuildUpVertex);
		m_pTTree->Branch("foundRecoLeptonInPrimaryVertex",&m_foundRecoLeptonInPrimaryVertex);
		m_pTTree->Branch("lostChargedMCP_CosTheta",&m_lostChargedMCP_CosTheta);
		m_pTTree->Branch("lostChargedMCP_Energy",&m_lostChargedMCP_Energy);
		m_pTTree->Branch("lostChargedMCP_Pt",&m_lostChargedMCP_Pt);
		m_pTTree->Branch("SLDecayXi", &m_SLDecayXi);
		m_pTTree->Branch("SLDecayYi", &m_SLDecayYi);
		m_pTTree->Branch("SLDecayZi", &m_SLDecayZi);
		m_pTTree->Branch("SLDecayRi", &m_SLDecayRi);
		m_pTTree->Branch("SLDecayXf", &m_SLDecayXf);
		m_pTTree->Branch("SLDecayYf", &m_SLDecayYf);
		m_pTTree->Branch("SLDecayZf", &m_SLDecayZf);
		m_pTTree->Branch("SLDecayRf", &m_SLDecayRf);
		m_pTTree->Branch("trueNuPx", &m_trueNuPx);
		m_pTTree->Branch("trueNuPy", &m_trueNuPy);
		m_pTTree->Branch("trueNuPz", &m_trueNuPz);
		m_pTTree->Branch("trueNuE", &m_trueNuE);
		m_pTTree->Branch("recoNuCloseInitialPx", &m_recoNuCloseInitialPx);
		m_pTTree->Branch("recoNuCloseInitialPy", &m_recoNuCloseInitialPy);
		m_pTTree->Branch("recoNuCloseInitialPz", &m_recoNuCloseInitialPz);
		m_pTTree->Branch("recoNuCloseInitialE", &m_recoNuCloseInitialE);
		m_pTTree->Branch("recoNuClosePx", &m_recoNuClosePx);
		m_pTTree->Branch("recoNuClosePy", &m_recoNuClosePy);
		m_pTTree->Branch("recoNuClosePz", &m_recoNuClosePz);
		m_pTTree->Branch("recoNuCloseE", &m_recoNuCloseE);
		m_pTTree->Branch("recoNuPosPx", &m_recoNuPosPx);
		m_pTTree->Branch("recoNuPosPy", &m_recoNuPosPy);
		m_pTTree->Branch("recoNuPosPz", &m_recoNuPosPz);
		m_pTTree->Branch("recoNuPosE", &m_recoNuPosE);
		m_pTTree->Branch("recoNuNegPx", &m_recoNuNegPx);
		m_pTTree->Branch("recoNuNegPy", &m_recoNuNegPy);
		m_pTTree->Branch("recoNuNegPz", &m_recoNuNegPz);
		m_pTTree->Branch("recoNuNegE", &m_recoNuNegE);
		m_pTTree->Branch("NuPxResidual", &m_NuPxResidual);
		m_pTTree->Branch("NuPyResidual", &m_NuPyResidual);
		m_pTTree->Branch("NuPzResidual", &m_NuPzResidual);
		m_pTTree->Branch("NuEResidual", &m_NuEResidual);
		m_pTTree->Branch("NuPxNormalizedResidual", &m_NuPxNormalizedResidual);
		m_pTTree->Branch("NuPyNormalizedResidual", &m_NuPyNormalizedResidual);
		m_pTTree->Branch("NuPzNormalizedResidual", &m_NuPzNormalizedResidual);
		m_pTTree->Branch("NuENormalizedResidual", &m_NuENormalizedResidual);
		m_pTTree->Branch("flightDirectionError", &m_FlightDirectionError);
		h_NuPxResidual = new TH1F( "PxResidual" , "; _{}p_{x,#nu}^{REC} - p_{x,#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 2000 , -10.0 , 10.0 ); n_NuPxResidual = 0;
		h_NuPyResidual = new TH1F( "PyResidual" , "; _{}p_{y,#nu}^{REC} - p_{y,#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 2000 , -10.0 , 10.0 ); n_NuPyResidual = 0;
		h_NuPzResidual = new TH1F( "PzResidual" , "; _{}p_{z,#nu}^{REC} - p_{z,#nu}^{MC}  [GeV]; Normalized Entries / 0.1" , 2000 , -10.0 , 10.0 ); n_NuPzResidual = 0;
		h_NuEResidual = new TH1F( "EResidual" , "; _{}E_{#nu}^{REC} - E_{#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 2000 , -10.0 , 10.0 ); n_NuEResidual = 0;
		h_NuPxNormalizedResidual = new TH1F( "PxNormalizedResidual" , "; ( _{}p_{x,#nu}^{REC} - p_{x,#nu}^{MC} ) / #sigma_{p_{x,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPxNormalizedResidual = 0;
		h_NuPyNormalizedResidual = new TH1F( "PyNormalizedResidual" , "; ( _{}p_{y,#nu}^{REC} - p_{y,#nu}^{MC} ) / #sigma_{p_{y,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPyNormalizedResidual = 0;
		h_NuPzNormalizedResidual = new TH1F( "PzNormalizedResidual" , "; ( _{}p_{z,#nu}^{REC} - p_{z,#nu}^{MC} ) / #sigma_{p_{z,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPzNormalizedResidual = 0;
		h_NuENormalizedResidual = new TH1F( "ENormalizedResidual" , "; ( _{}E_{#nu}^{REC} - E_{#nu}^{MC} ) / #sigma_{E_{#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuENormalizedResidual = 0;
		h_recoNuPx_mcNuPx = new TH2F( "p_{x}^{#nu}" , "; _{}p_{x,#nu}^{MC} [GeV] ;  _{}p_{x,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuPy_mcNuPy = new TH2F( "p_{y}^{#nu}" , "; _{}p_{y,#nu}^{MC} [GeV] ;  _{}p_{y,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuPz_mcNuPz = new TH2F( "p_{z}^{#nu}" , "; _{}p_{z,#nu}^{MC} [GeV] ;  _{}p_{z,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
		h_recoNuE_mcNuE = new TH2F( "E^{#nu}" , "; _{}E_{#nu}^{MC} [GeV] ;  _{}E_{#nu}^{REC} [GeV]" , 100 , 0.0 , 100.0 , 100 , 0.0 , 100.0 );
		h_parentPx_daughtersPx = new TH2F( "p_{x} conservation" , "; _{}Px_{parentHadron}^{MC} [GeV] ;  _{}Px_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentPy_daughtersPy = new TH2F( "p_{y} conservation" , "; _{}Py_{parentHadron}^{MC} [GeV] ;  _{}Py_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentPz_daughtersPz = new TH2F( "p_{z} conservation" , "; _{}Pz_{parentHadron}^{MC} [GeV] ;  _{}Pz_{daughters}^{MC} [GeV]" , 200 , -100.0 , 100.0 , 200 , -100.0 , 100.0 );
		h_parentE_daughtersE = new TH2F( "E conservation" , "; _{}E_{parentHadron}^{MC} [GeV] ;  _{}E_{daughters}^{MC} [GeV]" , 100 , 0.0 , 100.0 , 100 , 0.0 , 100.0 );
		h_recoPFOLinkedToElectron_Type = new TH1I( "PFOTypeofTrueElectron" , "; PFO Type" , 8 , 0 , 8 );
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(1,"e^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(3,"#pi^{#pm}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(4,"#gamma");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(5,"K^{0}_{S}");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(6,"#Lambda");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(7,"Other");
		h_recoPFOLinkedToElectron_Type->GetXaxis()->SetBinLabel(8,"Not Found");
		h_recoPFOLinkedToMuon_Type = new TH1I( "PFOTypeofTrueMuon" , "; PFO Type" , 8 , 0 , 8 );
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(1,"e^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(3,"#pi^{#pm}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(4,"#gamma");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(5,"K^{0}_{S}");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(6,"#Lambda");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(7,"Other");
		h_recoPFOLinkedToMuon_Type->GetXaxis()->SetBinLabel(8,"Not Found");
		h_SLDecayOrder = new TH1I( "SLDecayOrder" , "; SLDecay Type" , 3 , 0 , 3 );
		h_SLDecayOrder->GetXaxis()->SetBinLabel(1,"upStream");
		h_SLDecayOrder->GetXaxis()->SetBinLabel(2,"primary");
		h_SLDecayOrder->GetXaxis()->SetBinLabel(3,"downStream");
		h_foundVertex = new TH2I( "Viertices" , "; primary vertex ; secondary vertex" , 2 , 0 , 2 , 2 , 0 , 2 );
		h_foundVertex->GetXaxis()->SetBinLabel(1,"vertex not found");
		h_foundVertex->GetXaxis()->SetBinLabel(2,"vertex found");
		h_foundVertex->GetYaxis()->SetBinLabel(1,"vertex not found");
		h_foundVertex->GetYaxis()->SetBinLabel(2,"vertex found");
		h_secondaryVertex = new TH1I( "secondary vertices" , ";" , 6 , 0 , 6 );
		h_secondaryVertex->GetXaxis()->SetBinLabel(1,"lep in BUp vtx");
		h_secondaryVertex->GetXaxis()->SetBinLabel(2,"lep in Prim vtx");
		h_secondaryVertex->GetXaxis()->SetBinLabel(3,"SLD with downStream vtx");
		h_secondaryVertex->GetXaxis()->SetBinLabel(4,"Sec. vtx not found");
		h_secondaryVertex->GetXaxis()->SetBinLabel(5,"reco lep not found");
		h_secondaryVertex->GetXaxis()->SetBinLabel(6,"other(?)");
		h_parentHadronCharge = new TH1I( "parentHadronCharge" , "; Parent Hadron Charge" , 5 , 0 , 5 );
		h_parentHadronCharge->GetXaxis()->SetBinLabel(1,"-2");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(2,"-1");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(3,"0");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(4,"1");
		h_parentHadronCharge->GetXaxis()->SetBinLabel(5,"2");
		h_MCPTracks = new TH1I( "chargedMCPTracks" , ";" , 2 , 0 , 2 );
		h_MCPTracks->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks->GetXaxis()->SetBinLabel(2,"track lost");
		h_MCPTracks_Eweighted = new TH1I( "chargedMCPTracks" , "Energy weighted;" , 2 , 0 , 2 );
		h_MCPTracks_Eweighted->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks_Eweighted->GetXaxis()->SetBinLabel(2,"track lost");
		h_MCPTracks_Ptweighted = new TH1I( "chargedMCPTracks" , "p_{T} weighted;" , 2 , 0 , 2 );
		h_MCPTracks_Ptweighted->GetXaxis()->SetBinLabel(1,"track found");
		h_MCPTracks_Ptweighted->GetXaxis()->SetBinLabel(2,"track lost");
	}
}

void SLDCorrection::Clear()
{
	m_jetFlavourPDG.clear();	//4: c-jet, 5: b-jet, 21: g-jet, 0:Other
	m_nSLD_chargedMCPwoTrack.clear();
	m_GenStatParentHadron.clear();
	m_ChargeParentHadron.clear();
	m_foundRecoLepton.clear();
	m_foundBuildUpVertex.clear();
	m_foundRecoLeptonInBuildUpVertex.clear();
	m_foundRecoLeptonInPrimaryVertex.clear();
	m_lostChargedMCP_CosTheta.clear();
	m_lostChargedMCP_Energy.clear();
	m_lostChargedMCP_Pt.clear();
	m_nTauSLDecay = 0;
	m_nTauNeutrino = 0;
	m_nNeutrino = 0;
	m_nChargedPFOwoTrack = 0;
	m_SLDecayXi.clear();
	m_SLDecayYi.clear();
	m_SLDecayZi.clear();
	m_SLDecayRi.clear();
	m_SLDecayXf.clear();
	m_SLDecayYf.clear();
	m_SLDecayZf.clear();
	m_SLDecayRf.clear();
	m_trueNuPx.clear();
	m_trueNuPy.clear();
	m_trueNuPz.clear();
	m_trueNuE.clear();
	m_recoNuCloseInitialPx.clear();
	m_recoNuCloseInitialPy.clear();
	m_recoNuCloseInitialPz.clear();
	m_recoNuCloseInitialE.clear();
	m_recoNuClosePx.clear();
	m_recoNuClosePy.clear();
	m_recoNuClosePz.clear();
	m_recoNuCloseE.clear();
	m_recoNuPosPx.clear();
	m_recoNuPosPy.clear();
	m_recoNuPosPz.clear();
	m_recoNuPosE.clear();
	m_recoNuNegPx.clear();
	m_recoNuNegPy.clear();
	m_recoNuNegPz.clear();
	m_recoNuNegE.clear();
	m_NuPxResidual.clear();
	m_NuPyResidual.clear();
	m_NuPzResidual.clear();
	m_NuEResidual.clear();
	m_NuPxNormalizedResidual.clear();
	m_NuPyNormalizedResidual.clear();
	m_NuPzNormalizedResidual.clear();
	m_NuENormalizedResidual.clear();
	m_FlightDirectionError.clear();
}

void SLDCorrection::processRunHeader()
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}

void SLDCorrection::processEvent( EVENT::LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
//	LCCollection *RecoJetCollection{};
	LCCollection *MCParticleCollection{};
	LCCollection *jetFlacourCollection{};
	int nTauNeutrino = 0;
	int m_bJet = 0;
	int m_cJet = 0;
	int m_gJet = 0;
	int m_other = 0;
	int jetFlavourPDG = 0;
	++m_nEvtSum;
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

        try
        {
		jetFlacourCollection = pLCEvent->getCollection( m_jetFlavour );
		m_bJet = jetFlacourCollection->getParameters().getIntVal("isDecayedTob"); if ( m_bJet == 1 ) jetFlavourPDG = 5;
		m_cJet = jetFlacourCollection->getParameters().getIntVal("isDecayedToc"); if ( m_cJet == 1 ) jetFlavourPDG = 4;
		m_gJet = jetFlacourCollection->getParameters().getIntVal("isDecayedTog"); if ( m_gJet == 1 ) jetFlavourPDG = 21;
		m_other = jetFlacourCollection->getParameters().getIntVal("isDecayedToother"); if ( m_other == 1 ) jetFlavourPDG = 0;
		if ( m_bJet == 1 && !m_includbJets ) return;
		if ( m_cJet == 1 && !m_includcJets ) return;
		if ( m_gJet == 1 && !m_includgJets ) return;
		if ( m_other == 1 && !m_includOthers ) return;

		MCParticleCollection = pLCEvent->getCollection( m_mcParticleCollection );
		int nMCP = MCParticleCollection->getNumberOfElements();
		for ( int i_mcp = 0 ; i_mcp < nMCP ; ++i_mcp )
		{
			MCParticle *testLepton = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_mcp ) );
			if ( abs( testLepton->getPDG() ) == 16 && ( testLepton->getGeneratorStatus() ) == 1 ) ++nTauNeutrino;
			bool primarySLDecay = false;
			bool downStreamSLDecay = false;
			bool upStreamSLDecay = false;
			bool isBHadronSLDecay = false;
			bool isCHadronSLDecay = false;
			bool isTauLeptonSLDecay = false;
			if ( ( abs( testLepton->getPDG() ) == 11 || abs( testLepton->getPDG() ) == 13 || abs( testLepton->getPDG() ) == 15 ) && ( testLepton->getGeneratorStatus() ) == 1 )
			{
				for ( long unsigned int i_parent = 0 ; i_parent < ( testLepton->getParents() ).size() ; ++i_parent )
				{
					MCParticle *parent = testLepton->getParents()[ i_parent ];
					primarySLDecay = hasPrimarySLDecay( parent );
					if ( primarySLDecay ) downStreamSLDecay = hasDownStreamSLDecay( parent );
					if ( primarySLDecay ) upStreamSLDecay = hasUpStreamSLDecay( parent );
				}
				if ( primarySLDecay )
				{
					std::vector< TLorentzVector > recoNeutrinoFourMomentum;
					std::vector< std::vector< float > > recoNeutrinoCovMat;
					isBHadronSLDecay = checkBHadronSLDecay( testLepton );
					isCHadronSLDecay = checkCHadronSLDecay( testLepton );
					isTauLeptonSLDecay = checkTauLeptonSLDecay( testLepton );
					if ( isBHadronSLDecay && !m_includeBSLD ) continue;
					if ( isCHadronSLDecay && !m_includeCSLD ) continue;
					if ( isTauLeptonSLDecay && !m_includeTSLD ) continue;
					streamlog_out(DEBUG0) << "" << std::endl;
					streamlog_out(DEBUG0) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					streamlog_out(DEBUG0) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<< Found a primary semi-leptonic decay >>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
					if ( downStreamSLDecay )
					{
						streamlog_out(DEBUG4) << "	There is/are downstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
						h_SLDecayOrder->Fill( 2.5 );
					}
					if ( upStreamSLDecay )
					{
						streamlog_out(DEBUG4) << "	There is/are upstream semi-leptonic(s) decay in primary semi-leptonic decay products" << std::endl;
						h_SLDecayOrder->Fill( 0.5 );
					}
					if ( !downStreamSLDecay && !upStreamSLDecay )
					{
						streamlog_out(DEBUG0) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						streamlog_out(DEBUG0) << "	<<<<<<<<<<<<<<<< There are no upstream and downstream semi-leptonic decay >>>>>>>>>>>>>>>>>" << std::endl;
						h_SLDecayOrder->Fill( 1.5 );
						m_jetFlavourPDG.push_back( jetFlavourPDG );
						doSLDCorrection( pLCEvent , testLepton );
					}
				}
			}
		}
		m_nTauNeutrino = nTauNeutrino;

		m_pTTree->Fill();
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
        }

}


bool SLDCorrection::hasPrimarySLDecay( MCParticle *parentHadron )
{
	bool hasSLDecay = false;
	if ( parentHadron->getGeneratorStatus() == 2 && ( floor( abs( parentHadron->getPDG() ) / 100 ) == 5 || ( floor( abs( parentHadron->getPDG() ) / 1000 ) == 5 ) || floor( abs( parentHadron->getPDG() ) / 100 ) == 4 || ( floor( abs( parentHadron->getPDG() ) / 1000 ) == 4 ) || ( abs( parentHadron->getPDG() ) == 15 ) ) )
	{
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
		{
			MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
			if ( daughter->getGeneratorStatus() == 1 )
			{
				if ( abs( daughter->getPDG() ) == 11 || abs( daughter->getPDG() ) == 13 || abs( daughter->getPDG() ) == 15 )
				{
					for ( long unsigned int i_NuCandidate = 0 ; i_NuCandidate < ( parentHadron->getDaughters() ).size() ; ++i_NuCandidate )
					{
						MCParticle *secondDaughter = parentHadron->getDaughters()[ i_NuCandidate ];
						if ( ( abs( secondDaughter->getPDG() ) == abs( daughter->getPDG() ) + 1 ) && secondDaughter->getGeneratorStatus() == 1 )
						{
							hasSLDecay = true;
						}
					}
				}
			}
			else if( abs( daughter->getPDG() ) == 15 )
			{
				for ( long unsigned int i_NuCandidate = 0 ; i_NuCandidate < ( parentHadron->getDaughters() ).size() ; ++i_NuCandidate )
				{
					MCParticle *secondDaughter = parentHadron->getDaughters()[ i_NuCandidate ];
					if ( ( abs( secondDaughter->getPDG() ) == abs( daughter->getPDG() ) + 1 ) && secondDaughter->getGeneratorStatus() == 1 )
					{
						hasSLDecay = true;
						m_nTauSLDecay += 1;
					}
				}
			}
		}
	}
	return hasSLDecay;
}

bool SLDCorrection::hasDownStreamSLDecay( MCParticle *parentHadron )
{
	bool hasSLDecay = false;
	bool primarySLDecay = false;
	bool downStreamSLDecay = false;
	for ( long unsigned int i_primaryDaughter = 0 ; i_primaryDaughter < ( parentHadron->getDaughters() ).size() ; ++i_primaryDaughter )
	{
		MCParticle *primaryDaughter = parentHadron->getDaughters()[ i_primaryDaughter ];
		primarySLDecay = primarySLDecay || hasPrimarySLDecay( primaryDaughter );
		downStreamSLDecay = downStreamSLDecay || hasDownStreamSLDecay( primaryDaughter );
	}
	hasSLDecay = primarySLDecay || downStreamSLDecay ;
	return hasSLDecay;
}

bool SLDCorrection::hasUpStreamSLDecay( MCParticle *parentHadron )
{
	bool hasSLDecay = false;
	bool primarySLDecay = false;
	bool upStreamSLDecay = false;
	for ( long unsigned int i_upperParent = 0 ; i_upperParent < ( parentHadron->getParents() ).size() ; ++i_upperParent )
	{
		MCParticle *upperParent = parentHadron->getParents()[ i_upperParent ];
		primarySLDecay = primarySLDecay || hasPrimarySLDecay( upperParent );
		upStreamSLDecay = upStreamSLDecay || hasUpStreamSLDecay( upperParent );
	}
	hasSLDecay = primarySLDecay || upStreamSLDecay ;
	return hasSLDecay;
}

bool SLDCorrection::checkBHadronSLDecay( MCParticle *SLDLepton )
{
	bool isBHadronSLDecay = false;
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	if ( floor( abs( parentHadron->getPDG() ) / 100 ) == 5 || floor( abs( parentHadron->getPDG() ) / 1000 ) == 5 ) isBHadronSLDecay = true;
	return isBHadronSLDecay;
}

bool SLDCorrection::checkCHadronSLDecay( MCParticle *SLDLepton )
{
	bool isCHadronSLDecay = false;
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	if ( floor( abs( parentHadron->getPDG() ) / 100 ) == 4 || floor( abs( parentHadron->getPDG() ) / 1000 ) == 4 ) isCHadronSLDecay = true;
	return isCHadronSLDecay;
}

bool SLDCorrection::checkTauLeptonSLDecay( MCParticle *SLDLepton )
{
	bool TauLeptonSLDecay = false;
	MCParticle *parent = SLDLepton->getParents()[ 0 ];
	if ( abs( parent->getPDG() ) == 15 ) TauLeptonSLDecay = true;
	return TauLeptonSLDecay;
}

TLorentzVector SLDCorrection::getParentHadron4mom( MCParticle *SLDLepton )
{
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	TLorentzVector parentHadron4mom( parentHadron->getMomentum()[ 0 ] , parentHadron->getMomentum()[ 1 ] , parentHadron->getMomentum()[ 2 ] , parentHadron->getEnergy() );
	return parentHadron4mom;
}

TLorentzVector SLDCorrection::getDaughters4mom( MCParticle *SLDLepton )
{
	TLorentzVector daughters4mom( 0.0 , 0.0 , 0.0 , 0.0 );
	MCParticle *parentHadron = SLDLepton->getParents()[ 0 ];
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
	{
		MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
		daughters4mom += TLorentzVector( daughter->getMomentum()[ 0 ] , daughter->getMomentum()[ 1 ] , daughter->getMomentum()[ 2 ] , daughter->getEnergy() );
	}
	return daughters4mom;
}

void SLDCorrection::doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	TLorentzVector recoNeutrinoFourMomentumPos( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumNeg( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumClose( 0.0 , 0.0 , 0.0 , 0.0 );
	streamlog_out(DEBUG0) << "			      ( PDG	, Mass		, Px		, Py		, Pz		, E		, Charge		)" << std::endl;
	TLorentzVector fourMomentumLepton = getLeptonFourMomentum( pLCEvent , SLDLepton );
	TLorentzVector visibleFourMomentumCharged = getVisibleFourMomentum( pLCEvent , SLDLepton , SLDLepton->getParents()[ 0 ] , true , false );
	TLorentzVector visibleFourMomentumNeutral = getVisibleFourMomentum( pLCEvent , SLDLepton , SLDLepton->getParents()[ 0 ] , false , true );
	TLorentzVector visibleFourMomentum = fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral;
	TLorentzVector trueNeutrinoFourMomentum = getTrueNeutrinoFourMomentum( SLDLepton );
	TLorentzVector totalFourMomentum = visibleFourMomentum + trueNeutrinoFourMomentum;

	m_nSLD_chargedMCPwoTrack.push_back( m_nChargedPFOwoTrack );
	m_GenStatParentHadron.push_back( ( SLDLepton->getParents()[ 0 ] )->getGeneratorStatus() );
	m_ChargeParentHadron.push_back( ( SLDLepton->getParents()[ 0 ] )->getCharge() );
	m_nChargedPFOwoTrack = 0;

//	Test for Energy-momentum conservation
	TLorentzVector parentHadron4mom = getParentHadron4mom( SLDLepton );
	TLorentzVector daughters4mom = getDaughters4mom( SLDLepton );
	daughters4mom = totalFourMomentum;
	h_parentPx_daughtersPx->Fill( parentHadron4mom.Px() , daughters4mom.Px() );
	h_parentPy_daughtersPy->Fill( parentHadron4mom.Py() , daughters4mom.Py() );
	h_parentPz_daughtersPz->Fill( parentHadron4mom.Pz() , daughters4mom.Pz() );
	h_parentE_daughtersE->Fill( parentHadron4mom.E() , daughters4mom.E() );

	double sigmaParentHadronMass = 0.0;
	std::vector< float > flightDirectionCovMat( 6 , 0.0 );
	std::vector< float > LeptonCovMat( 10 , 0.0 );
	std::vector< float > visibleChargedCovMat( 10 , 0.0 );
	std::vector< float > visibleNeutralCovMat( 10 , 0.0 );
	std::vector< float > NeutrinoCovMat( 10 , 0.0 );

	streamlog_out(DEBUG0) << "	-----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG4) << "	Visible Lepton 4-Momentum:			( " << fourMomentumLepton.Px() << "	, " << fourMomentumLepton.Py() << "	, " << fourMomentumLepton.Pz() << "	, " << fourMomentumLepton.E() << " )" << std::endl;
	streamlog_out(DEBUG4) << "	Visible Charged 4-Momentum:			( " << visibleFourMomentumCharged.Px() << "	, " << visibleFourMomentumCharged.Py() << "	, " << visibleFourMomentumCharged.Pz() << "	, " << visibleFourMomentumCharged.E() << " )" << std::endl;
	streamlog_out(DEBUG4) << "	Visible Neutral 4-Momentum:			( " << visibleFourMomentumNeutral.Px() << "	, " << visibleFourMomentumNeutral.Py() << "	, " << visibleFourMomentumNeutral.Pz() << "	, " << visibleFourMomentumNeutral.E() << " )" << std::endl;
	streamlog_out(DEBUG0) << "" << std::endl;
	streamlog_out(DEBUG4) << "	Visible Total 4-Momentum:			( " << visibleFourMomentum.Px() << "	, " << visibleFourMomentum.Py() << "	, " << visibleFourMomentum.Pz() << "	, " << visibleFourMomentum.E() << " )" << std::endl;
	streamlog_out(DEBUG4) << "	Invisible Neutrino 4-Momentum:		( " << trueNeutrinoFourMomentum.Px() << "	, " << trueNeutrinoFourMomentum.Py() << "	, " << trueNeutrinoFourMomentum.Pz() << "	, " << trueNeutrinoFourMomentum.E() << " )" << std::endl;
	streamlog_out(DEBUG4) << "	Visible + Invisible 4-Momentum:		( " << totalFourMomentum.Px() << "	, " << totalFourMomentum.Py() << "	, " << totalFourMomentum.Pz() << "	, " << totalFourMomentum.E() << " )" << std::endl;
	streamlog_out(DEBUG0) << "	-----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	double parentHadronMass = getParentHadronMass( SLDLepton );
	foundFlightDirection = true;
	TVector3 flightDirection = getFlightDirection( pLCEvent , SLDLepton );
	if ( !foundFlightDirection ) return;
	m_FlightDirectionError.push_back( flightDirection.Dot( cheatFlightDirection( SLDLepton ) ) );

	recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , +1 );
	recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , -1 );

	recoNeutrinoFourMomentumClose = ( fabs( recoNeutrinoFourMomentumPos.E() - trueNeutrinoFourMomentum.E() ) < fabs( recoNeutrinoFourMomentumNeg.E() - trueNeutrinoFourMomentum.E() ) ? recoNeutrinoFourMomentumPos : recoNeutrinoFourMomentumNeg );
	streamlog_out(DEBUG4) << "" << std::endl;
	streamlog_out(DEBUG4) << "	Closest Neutrino 4-Momentum:			( " << recoNeutrinoFourMomentumClose.Px() << "	, " << recoNeutrinoFourMomentumClose.Py() << "	, " << recoNeutrinoFourMomentumClose.Pz() << "	, " << recoNeutrinoFourMomentumClose.E() << " )" << std::endl;
	streamlog_out(DEBUG4) << "	True Neutrino 4-Momentum:			( " << trueNeutrinoFourMomentum.Px() << "	, " << trueNeutrinoFourMomentum.Py() << "	, " << trueNeutrinoFourMomentum.Pz() << "	, " << trueNeutrinoFourMomentum.E() << " )" << std::endl;
	if ( fabs( recoNeutrinoFourMomentumClose.E() - trueNeutrinoFourMomentum.E() ) > 10.0 )
	{
		streamlog_out(DEBUG4) << "	!!! Big Difference between true and reco neutrino Energy : " << recoNeutrinoFourMomentumClose.E() - trueNeutrinoFourMomentum.E() << "  GeV" << std::endl;
	}
	streamlog_out(DEBUG4) << "	----------------------------------------------------------------------------------------------------" << std::endl;
	NeutrinoCovMat = getNeutrinoCovMat( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , flightDirectionCovMat , LeptonCovMat , visibleChargedCovMat , visibleNeutralCovMat , sigmaParentHadronMass );
	streamlog_out(DEBUG4) << "	Reconstructed Neutrino CovMat:" << std::endl;
	streamlog_out(DEBUG4) << "						" << NeutrinoCovMat[ 0 ] << std::endl;
	streamlog_out(DEBUG4) << "						" << NeutrinoCovMat[ 1 ] << "	, " << NeutrinoCovMat[ 2 ] << std::endl;
	streamlog_out(DEBUG4) << "						" << NeutrinoCovMat[ 3 ] << "	, " << NeutrinoCovMat[ 4 ]  << "	, " << NeutrinoCovMat[ 5 ] << std::endl;
	streamlog_out(DEBUG4) << "						" << NeutrinoCovMat[ 6 ] << "	, " << NeutrinoCovMat[ 7 ]  << "	, " << NeutrinoCovMat[ 8 ] << "	, " << NeutrinoCovMat[ 9 ]  << std::endl;
	streamlog_out(DEBUG4) << "	----------------------------------------------------------------------------------------------------" << std::endl;
	plotHistograms( trueNeutrinoFourMomentum , recoNeutrinoFourMomentumClose , NeutrinoCovMat );
	m_trueNuPx.push_back( trueNeutrinoFourMomentum.Px() );
	m_trueNuPy.push_back( trueNeutrinoFourMomentum.Py() );
	m_trueNuPz.push_back( trueNeutrinoFourMomentum.Pz() );
	m_trueNuE.push_back( trueNeutrinoFourMomentum.E() );
	m_recoNuClosePx.push_back( recoNeutrinoFourMomentumClose.Px() );
	m_recoNuClosePy.push_back( recoNeutrinoFourMomentumClose.Py() );
	m_recoNuClosePz.push_back( recoNeutrinoFourMomentumClose.Pz() );
	m_recoNuCloseE.push_back( recoNeutrinoFourMomentumClose.E() );
	m_recoNuPosPx.push_back( recoNeutrinoFourMomentumPos.Px() );
	m_recoNuPosPy.push_back( recoNeutrinoFourMomentumPos.Py() );
	m_recoNuPosPz.push_back( recoNeutrinoFourMomentumPos.Pz() );
	m_recoNuPosE.push_back( recoNeutrinoFourMomentumPos.E() );
	m_recoNuNegPx.push_back( recoNeutrinoFourMomentumNeg.Px() );
	m_recoNuNegPy.push_back( recoNeutrinoFourMomentumNeg.Py() );
	m_recoNuNegPz.push_back( recoNeutrinoFourMomentumNeg.Pz() );
	m_recoNuNegE.push_back( recoNeutrinoFourMomentumNeg.E() );
}

TVector3 SLDCorrection::getFlightDirection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	TVector3 flightDirection( 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumPos( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumNeg( 0.0 , 0.0 , 0.0 , 0.0 );
	TLorentzVector recoNeutrinoFourMomentumClose( 0.0 , 0.0 , 0.0 , 0.0 );
	double parentHadronMass = getParentHadronMass( SLDLepton );
	TLorentzVector trueNeutrinoFourMomentum = getTrueNeutrinoFourMomentum( SLDLepton );
	TLorentzVector fourMomentumLepton = getLeptonFourMomentum( pLCEvent , SLDLepton );
	TLorentzVector visibleFourMomentumCharged = getVisibleFourMomentum( pLCEvent , SLDLepton , SLDLepton->getParents()[ 0 ] , true , false );
	TLorentzVector visibleFourMomentumNeutral = getVisibleFourMomentum( pLCEvent , SLDLepton , SLDLepton->getParents()[ 0 ] , false , true );
	TLorentzVector parentFourMomentumPos = fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral;
	TLorentzVector parentFourMomentumNeg = fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral;
	std::vector< double > primaryVertex = getPrimaryVertex( pLCEvent , SLDLepton );
	std::vector< double > secondaryVertex = getSecondaryVertex( pLCEvent , SLDLepton );
	bool foundPrimaryVertex = ( primaryVertex[ 3 ] > 0 ? true : false );
	bool foundSecondaryVertex = ( secondaryVertex[ 3 ] > 0 ? true : false );
	h_foundVertex->Fill( ( foundPrimaryVertex ? 1.5 : 0.5 ) , ( foundSecondaryVertex ? 1.5 : 0.5 ) );
	int parentCharge = getParentCharge( SLDLepton );
	if ( parentCharge == -2 )
	{
		h_parentHadronCharge->Fill( 0.5 );
	}
	else if ( parentCharge == -1 )
	{
		h_parentHadronCharge->Fill( 1.5 );
	}
	else if ( parentCharge == 0 )
	{
		h_parentHadronCharge->Fill( 2.5 );
	}
	else if ( parentCharge == 1 )
	{
		h_parentHadronCharge->Fill( 3.5 );
	}
	else if ( parentCharge == 2 )
	{
		h_parentHadronCharge->Fill( 4.5 );
	}
	streamlog_out(DEBUG0) << "	-----------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	if ( m_cheatFlightDirection )
	{
		foundFlightDirection = true;
		flightDirection = cheatFlightDirection( SLDLepton );
	}
	else if ( m_useJetAxisAsFlightDirection )
	{
		flightDirection = getJetAxis( pLCEvent , SLDLepton );
		if ( flightDirection.Mag() > 0.5 )
		{
			foundFlightDirection = true;
		}
		else
		{
			foundFlightDirection = false;
		}
	}
	else if ( foundPrimaryVertex )
	{
		if ( foundSecondaryVertex )
		{
			foundFlightDirection = true;
			if ( m_considerParentCharge )
			{
				flightDirection = getParentHadronUnitMomentum( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral , primaryVertex , secondaryVertex , 0 );
				float flightDistance = sqrt( pow( secondaryVertex[ 0 ] - primaryVertex[ 0 ] , 2 ) + pow( secondaryVertex[ 1 ] - primaryVertex[ 1 ] , 2 ) + pow( secondaryVertex[ 2 ] - primaryVertex[ 2 ] , 2 ) );
				if ( flightDistance == 0.0 )
				{
					flightDirection = ( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral ).Vect();
					flightDirection.SetMag( 1.0 );
				}
				recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , +1 );
				recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , -1 );
				recoNeutrinoFourMomentumClose = ( fabs( recoNeutrinoFourMomentumPos.E() - trueNeutrinoFourMomentum.E() ) < fabs( recoNeutrinoFourMomentumNeg.E() - trueNeutrinoFourMomentum.E() ) ? recoNeutrinoFourMomentumPos : recoNeutrinoFourMomentumNeg );
				m_recoNuCloseInitialPx.push_back( recoNeutrinoFourMomentumClose.Px() );
				m_recoNuCloseInitialPy.push_back( recoNeutrinoFourMomentumClose.Py() );
				m_recoNuCloseInitialPz.push_back( recoNeutrinoFourMomentumClose.Pz() );
				m_recoNuCloseInitialE.push_back( recoNeutrinoFourMomentumClose.E() );
				for ( int i_iter = 0 ; i_iter < m_nIterFlightDirCorrection ; ++i_iter )
				{
					flightDirection = getParentHadronUnitMomentum( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral + recoNeutrinoFourMomentumPos , primaryVertex , secondaryVertex , parentCharge );
					if ( flightDistance == 0.0 )
					{
						flightDirection = ( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral + recoNeutrinoFourMomentumPos ).Vect();
						flightDirection.SetMag( 1.0 );
					}
					recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , +1 );
					flightDirection = getParentHadronUnitMomentum( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral + recoNeutrinoFourMomentumNeg , primaryVertex , secondaryVertex , parentCharge );
					if ( flightDistance == 0.0 )
					{
						flightDirection = ( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral + recoNeutrinoFourMomentumNeg ).Vect();
						flightDirection.SetMag( 1.0 );
					}
					recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , -1 );
				}
			}
			else
			{
				flightDirection = getParentHadronUnitMomentum( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral , primaryVertex , secondaryVertex , 0 );
				float flightDistance = sqrt( pow( secondaryVertex[ 0 ] - primaryVertex[ 0 ] , 2 ) + pow( secondaryVertex[ 1 ] - primaryVertex[ 1 ] , 2 ) + pow( secondaryVertex[ 2 ] - primaryVertex[ 2 ] , 2 ) );
				if ( flightDistance == 0.0 )
				{
					flightDirection = ( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral ).Vect();
					flightDirection.SetMag( 1.0 );
				}
				recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , +1 );
				recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , -1 );
				recoNeutrinoFourMomentumClose = ( fabs( recoNeutrinoFourMomentumPos.E() - trueNeutrinoFourMomentum.E() ) < fabs( recoNeutrinoFourMomentumNeg.E() - trueNeutrinoFourMomentum.E() ) ? recoNeutrinoFourMomentumPos : recoNeutrinoFourMomentumNeg );
				m_recoNuCloseInitialPx.push_back( recoNeutrinoFourMomentumClose.Px() );
				m_recoNuCloseInitialPy.push_back( recoNeutrinoFourMomentumClose.Py() );
				m_recoNuCloseInitialPz.push_back( recoNeutrinoFourMomentumClose.Pz() );
				m_recoNuCloseInitialE.push_back( recoNeutrinoFourMomentumClose.E() );
//					flightDirection = getParentHadronUnitMomentum( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral + recoNeutrinoFourMomentumPos , primaryVertex , secondaryVertex , 0 );
				if ( flightDistance == 0.0 )
				{
					flightDirection = ( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral + recoNeutrinoFourMomentumPos ).Vect();
					flightDirection.SetMag( 1.0 );
				}
				recoNeutrinoFourMomentumPos = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , +1 );
//				flightDirection = getParentHadronUnitMomentum( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral + recoNeutrinoFourMomentumNeg , primaryVertex , secondaryVertex , 0 );
				if ( flightDistance == 0.0 )
				{
					flightDirection = ( fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral + recoNeutrinoFourMomentumNeg ).Vect();
					flightDirection.SetMag( 1.0 );
				}
				recoNeutrinoFourMomentumNeg = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , -1 );
			}
		}
		else
		{

		}
	}
	else
	{
		foundFlightDirection = false;
	}
	return flightDirection;
}

std::vector< double > SLDCorrection::getPrimaryVertex( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	std::vector< double > primaryVertex( 4 , 0.0 );
	primaryVertex[ 3 ] = -1.0;// to clarify whether vertex is found or not?
	if ( m_cheatVertices )
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		primaryVertex[ 0 ] = MotherHadron->getVertex()[ 0 ];
		primaryVertex[ 1 ] = MotherHadron->getVertex()[ 1 ];
		primaryVertex[ 2 ] = MotherHadron->getVertex()[ 2 ];
		primaryVertex[ 3 ] = 1.0;
		streamlog_out(DEBUG4) << "		true primary Vertex (x,y,z): 		" << primaryVertex[ 0 ] << "	, " << primaryVertex[ 1 ] << "	, " << primaryVertex[ 2 ] << std::endl;
	}
	else
	{
		LCCollection *primaryVertexCollection = pLCEvent->getCollection( m_inputPrimaryVertex );
		Vertex* primaryVtx = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
		primaryVertex[ 0 ] = primaryVtx->getPosition()[ 0 ];
		primaryVertex[ 1 ] = primaryVtx->getPosition()[ 1 ];
		primaryVertex[ 2 ] = primaryVtx->getPosition()[ 2 ];
		primaryVertex[ 3 ] = 1.0;
		streamlog_out(DEBUG4) << "		reco primary Vertex (x,y,z): 	" << primaryVertex[ 0 ] << "	, " << primaryVertex[ 1 ] << "	, " << primaryVertex[ 2 ] << std::endl;
	}
	m_SLDecayXi.push_back( primaryVertex[ 0 ] );
	m_SLDecayYi.push_back( primaryVertex[ 1 ] );
	m_SLDecayZi.push_back( primaryVertex[ 2 ] );
	m_SLDecayRi.push_back( sqrt( pow( primaryVertex[ 0 ] , 2 ) + pow( primaryVertex[ 1 ] , 2 ) + pow( primaryVertex[ 2 ] , 2 ) ) );
	return primaryVertex;
}

std::vector< double > SLDCorrection::getSecondaryVertex( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	std::vector< double > secondaryVertex( 4 , 0.0 );
	bool foundRecoLeptonInBuildUpVertex = false;
	bool foundRecoLeptonInPrimaryVertex = false;
	secondaryVertex[ 3 ] = -1.0;
	// to clarify whether vertex is found or not?
	// -1: not found Sec.Vtx
	//  0: cheated Sec.Vtx
	//  1: found RecoLepton in Sec.Vtx
	//  2: use downStream Vtx and intersect with RecoLepton as Sec.Vtx
	//  3: Reco Lepton in Prim. Vtx, use jet asix for Sec.Vtx
	//  4: No Sec.Vtx, jet axis is used
	if ( m_cheatVertices )
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		secondaryVertex[ 0 ] = MotherHadron->getEndpoint()[ 0 ];
		secondaryVertex[ 1 ] = MotherHadron->getEndpoint()[ 1 ];
		secondaryVertex[ 2 ] = MotherHadron->getEndpoint()[ 2 ];
		secondaryVertex[ 3 ] = 0.0;
		streamlog_out(DEBUG4) << "		true secondary Vertex (x,y,z): 	" << secondaryVertex[ 0 ] << "	, " << secondaryVertex[ 1 ] << "	, " << secondaryVertex[ 2 ] << std::endl;
	}
	else
	{
		try
		{
			LCCollection *BuildUpVertexCollection = pLCEvent->getCollection( m_inputBuildUpVertex );
			streamlog_out(DEBUG2) << "	There are " << BuildUpVertexCollection->getNumberOfElements() << " BuildUp Vertices" << std::endl;
		}
		catch(DataNotAvailableException &e)
	        {
	        	streamlog_out(MESSAGE) << "	BuildUp Vertex collection not found" << std::endl;
			return secondaryVertex;
	        }
		ReconstructedParticle* linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , true , false );
		m_foundRecoLepton.push_back( ( linkedRecoLepton == NULL ? 0 : 1 ) );
		LCCollection *BuildUpVertexCollection = pLCEvent->getCollection( m_inputBuildUpVertex );
		int n_VTX = BuildUpVertexCollection->getNumberOfElements();
		m_foundBuildUpVertex.push_back( ( n_VTX == 0 ? 0 : 1 ) );
		if ( linkedRecoLepton != NULL )
		{
			for ( int i_vtx = 0 ; i_vtx < n_VTX ; ++i_vtx )
			{
				Vertex* secondaryVtx = dynamic_cast<Vertex*>( BuildUpVertexCollection->getElementAt( i_vtx ) );
				ReconstructedParticle* parent = secondaryVtx->getAssociatedParticle();
				int n_daughters = ( parent->getParticles() ).size();
				for ( int i_daughter = 0 ; i_daughter < n_daughters ; ++i_daughter )
				{
					ReconstructedParticle* daughter = parent->getParticles()[ i_daughter ];
					if ( linkedRecoLepton == daughter )
					{
						secondaryVertex[ 0 ] = secondaryVtx->getPosition()[ 0 ];
						secondaryVertex[ 1 ] = secondaryVtx->getPosition()[ 1 ];
						secondaryVertex[ 2 ] = secondaryVtx->getPosition()[ 2 ];
						secondaryVertex[ 3 ] = 1.0;
						foundRecoLeptonInBuildUpVertex = true;
						streamlog_out(DEBUG0) << "	Found Reco Lepton in BuildUp Vertex" << std::endl;
					}
				}
			}
			m_foundRecoLeptonInBuildUpVertex.push_back( ( foundRecoLeptonInBuildUpVertex ? 1 : 0 ) );
			if ( !foundRecoLeptonInBuildUpVertex )
			{
				if ( hasDownStreamVertex( SLDLepton->getParents()[ 0 ] ) )
				{
					secondaryVertex = getDownStreamVertex( pLCEvent , SLDLepton->getParents()[ 0 ] );
					if ( secondaryVertex[ 3 ] == 1.0 )
					{
						secondaryVertex[ 3 ] = 2.0;
					}
					else
					{
						secondaryVertex[ 3 ] = -1.0;
					}
				}
				if ( hasDownStreamVertex( SLDLepton->getParents()[ 0 ] ) || secondaryVertex[ 3 ] == -1.0 )
				{
					secondaryVertex[ 0 ] = getJetAxis( pLCEvent , SLDLepton ).X() + getPrimaryVertex( pLCEvent , SLDLepton )[ 0 ];
					secondaryVertex[ 1 ] = getJetAxis( pLCEvent , SLDLepton ).Y() + getPrimaryVertex( pLCEvent , SLDLepton )[ 1 ];
					secondaryVertex[ 2 ] = getJetAxis( pLCEvent , SLDLepton ).Z() + getPrimaryVertex( pLCEvent , SLDLepton )[ 2 ];

					LCCollection *primaryVertexCollection = pLCEvent->getCollection( m_inputPrimaryVertex );
					Vertex* primaryVtx = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
					ReconstructedParticle* parent = primaryVtx->getAssociatedParticle();
					int n_daughters = ( parent->getParticles() ).size();
					for ( int i_daughter = 0 ; i_daughter < n_daughters ; ++i_daughter )
					{
						ReconstructedParticle* daughter = parent->getParticles()[ i_daughter ];
						if ( linkedRecoLepton == daughter )
						{
							foundRecoLeptonInPrimaryVertex = true;
							secondaryVertex[ 3 ] = 3.0;
							streamlog_out(DEBUG0) << "	Found Reco Lepton in Primary Vertex" << std::endl;
						}
					}
					if ( !foundRecoLeptonInPrimaryVertex ) secondaryVertex[ 3 ] = 4.0;
				}
			}
			m_foundRecoLeptonInPrimaryVertex.push_back( ( foundRecoLeptonInPrimaryVertex ? 1 : 0 ) );

		}

		if ( foundRecoLeptonInBuildUpVertex )
		{
			h_secondaryVertex->Fill( 0.5 );
		}
		else if ( foundRecoLeptonInPrimaryVertex )
		{
			h_secondaryVertex->Fill( 1.5 );
		}
		else if ( linkedRecoLepton != NULL )
		{
			streamlog_out(DEBUG0) << "	Not Found Reco Lepton in BuildUp/Primary Vertex" << std::endl;
			if ( linkedRecoLepton->getStartVertex() != NULL ) streamlog_out(DEBUG0) << "		StartVertex of Reco Lepton [" << linkedRecoLepton->getStartVertex() << "] at (x,y,z) = (" << ( linkedRecoLepton->getStartVertex() )->getPosition()[ 0 ] << " , " << ( linkedRecoLepton->getStartVertex() )->getPosition()[ 1 ] << " , " << ( linkedRecoLepton->getStartVertex() )->getPosition()[ 2 ] << " )" << std::endl;
			if ( hasDownStreamVertex( SLDLepton->getParents()[ 0 ] ) )
			{
				h_secondaryVertex->Fill( 2.5 );
			}
			else
			{
				h_secondaryVertex->Fill( 3.5 );
			}
		}
		else if ( linkedRecoLepton == NULL )
		{
			streamlog_out(DEBUG0) << "	Not Found Reco Lepton linked to mcLepton" << std::endl;
			h_secondaryVertex->Fill( 4.5 );
		}
		else
		{
			h_secondaryVertex->Fill( 5.5 );
		}
	}
	m_SLDecayXf.push_back( secondaryVertex[ 0 ] );
	m_SLDecayYf.push_back( secondaryVertex[ 1 ] );
	m_SLDecayZf.push_back( secondaryVertex[ 2 ] );
	m_SLDecayRf.push_back( sqrt( pow( secondaryVertex[ 0 ] , 2 ) + pow( secondaryVertex[ 1 ] , 2 ) + pow( secondaryVertex[ 2 ] , 2 ) ) );
	return secondaryVertex;
}

bool SLDCorrection::hasDownStreamVertex( MCParticle *neutralParticle )
{
	bool hasPrimVertex = false;
	bool hasDSVertex = false;
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( neutralParticle->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = neutralParticle->getDaughters()[ i_daughter ];
		hasPrimVertex = hasPrimVertex || hasVertex( daughter );
		hasDSVertex = hasDSVertex || hasDownStreamVertex( daughter );
	}
	return hasPrimVertex || hasDSVertex;
}

bool SLDCorrection::hasVertex( MCParticle *neutralParticle )
{
	bool hasChargedVertex = false;
	for ( long unsigned int i_daughter1 = 0 ; i_daughter1 < ( neutralParticle->getDaughters() ).size() ; ++i_daughter1 )
	{
		EVENT::MCParticle *daughter1 = neutralParticle->getDaughters()[ i_daughter1 ];
		if ( daughter1->getCharge() < -0.5 && daughter1->getGeneratorStatus() == 1 )
		{
			for ( long unsigned int i_daughter2 = 0 ; i_daughter2 < ( neutralParticle->getDaughters() ).size() ; ++i_daughter2 )
			{
				EVENT::MCParticle *daughter2 = neutralParticle->getDaughters()[ i_daughter2 ];
				if ( daughter2->getCharge() > 0.5 && daughter2->getGeneratorStatus() == 1 )
				{
					hasChargedVertex = true;
				}
			}
		}
	}
	return hasChargedVertex;
}

std::vector<double> SLDCorrection::getDownStreamVertex( EVENT::LCEvent *pLCEvent , MCParticle *neutralParticle )
{
	std::vector< double > downStreamVertex( 4 , 0.0 );
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( neutralParticle->getDaughters() ).size() ; ++i_daughter )
	{
		EVENT::MCParticle *daughter = neutralParticle->getDaughters()[ i_daughter ];
		if ( hasVertex( daughter ) )
		{
			downStreamVertex = getVertex( pLCEvent , daughter );
		}
		else
		{
			downStreamVertex = getDownStreamVertex( pLCEvent , daughter );
		}
	}
	return downStreamVertex;
}

std::vector<double> SLDCorrection::getVertex( EVENT::LCEvent *pLCEvent , MCParticle *neutralParticle )
{
	std::vector< double > downStreamVertex( 4 , 0.0 );
	downStreamVertex[ 3 ] = -1.0;
	for ( long unsigned int i_daughter1 = 0 ; i_daughter1 < ( neutralParticle->getDaughters() ).size() ; ++i_daughter1 )
	{
		EVENT::MCParticle *daughter1 = neutralParticle->getDaughters()[ i_daughter1 ];
		if ( daughter1->getCharge() < -0.5 && daughter1->getGeneratorStatus() == 1 )
		{
			for ( long unsigned int i_daughter2 = 0 ; i_daughter2 < ( neutralParticle->getDaughters() ).size() ; ++i_daughter2 )
			{
				EVENT::MCParticle *daughter2 = neutralParticle->getDaughters()[ i_daughter2 ];
				if ( daughter2->getCharge() > 0.5 && daughter2->getGeneratorStatus() == 1 )
				{
					ReconstructedParticle* linkedChargedPFO1 = getLinkedPFO( pLCEvent , daughter1 , true , false );
					ReconstructedParticle* linkedChargedPFO2 = getLinkedPFO( pLCEvent , daughter2 , true , false );
					LCCollection *BuildUpVertexCollection = pLCEvent->getCollection( m_inputBuildUpVertex );
					int n_VTX = BuildUpVertexCollection->getNumberOfElements();
					if ( linkedChargedPFO1 != NULL && linkedChargedPFO2 != NULL )
					{
						for ( int i_vtx = 0 ; i_vtx < n_VTX ; ++i_vtx )
						{
							Vertex* downStreamVtx = dynamic_cast<Vertex*>( BuildUpVertexCollection->getElementAt( i_vtx ) );
							ReconstructedParticle* parent = downStreamVtx->getAssociatedParticle();
							int n_daughters = ( parent->getParticles() ).size();
							for ( int i_recoDaughter1 = 0 ; i_recoDaughter1 < n_daughters ; ++i_recoDaughter1 )
							{
								ReconstructedParticle* recoDaughter1 = parent->getParticles()[ i_recoDaughter1 ];
								if ( linkedChargedPFO1 == recoDaughter1 )
								{
									for ( int i_recoDaughter2 = 0 ; i_recoDaughter2 < n_daughters ; ++i_recoDaughter2 )
									{
										ReconstructedParticle* recoDaughter2 = parent->getParticles()[ i_recoDaughter2 ];
										if ( linkedChargedPFO2 == recoDaughter2 )
										{
											downStreamVertex[ 0 ] = downStreamVtx->getPosition()[ 0 ];
											downStreamVertex[ 1 ] = downStreamVtx->getPosition()[ 1 ];
											downStreamVertex[ 2 ] = downStreamVtx->getPosition()[ 2 ];
											downStreamVertex[ 3 ] = 1.0;
											streamlog_out(DEBUG0) << "	Found Down Stream Vertex in BuildUp Vertex" << std::endl;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return downStreamVertex;
}

/*
std::vector< double > SLDCorrection::getDownStreamVertex( MCParticle *parentHadron , bool foundDownStreamVertex )
{
	std::vector< double > downStreamVertex( 4 , 0.0 );
	downStreamVertex[ 3 ] = -1.0;// to clarify whether vertex is found or not?
	for ( int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
	{
		const EVENT::MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
		if ( fabs( daughter->getCharge() ) < 0.5 )
		{

		}
	}




	if ( m_cheatVertices )
	{
		primaryVertex[ 0 ] = MotherHadron->getVertex()[ 0 ];
		primaryVertex[ 1 ] = MotherHadron->getVertex()[ 1 ];
		primaryVertex[ 2 ] = MotherHadron->getVertex()[ 2 ];
		primaryVertex[ 3 ] = 1.0;
		streamlog_out(DEBUG4) << "		true primary Vertex (x,y,z): 		" << primaryVertex[ 0 ] << "	, " << primaryVertex[ 1 ] << "	, " << primaryVertex[ 2 ] << std::endl;
	}
	else
	{
		LCCollection *primaryVertexCollection = pLCEvent->getCollection( m_inputPrimaryVertex );
		Vertex* primaryVtx = dynamic_cast<Vertex*>( primaryVertexCollection->getElementAt( 0 ) );
		primaryVertex[ 0 ] = primaryVtx->getPosition()[ 0 ];
		primaryVertex[ 1 ] = primaryVtx->getPosition()[ 1 ];
		primaryVertex[ 2 ] = primaryVtx->getPosition()[ 2 ];
		primaryVertex[ 3 ] = 1.0;
		streamlog_out(DEBUG4) << "		reco primary Vertex (x,y,z): 	" << primaryVertex[ 0 ] << "	, " << primaryVertex[ 1 ] << "	, " << primaryVertex[ 2 ] << std::endl;
	}
	m_SLDecayXi.push_back( primaryVertex[ 0 ] );
	m_SLDecayYi.push_back( primaryVertex[ 1 ] );
	m_SLDecayZi.push_back( primaryVertex[ 2 ] );
	m_SLDecayRi.push_back( sqrt( pow( primaryVertex[ 0 ] , 2 ) + pow( primaryVertex[ 1 ] , 2 ) + pow( primaryVertex[ 2 ] , 2 ) ) );
	return primaryVertex;
}
*/
int SLDCorrection::getParentCharge( MCParticle *SLDLepton )
{
	const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	int parentCharge = MotherHadron->getCharge();
	streamlog_out(DEBUG0) << "		parent Hadron Charge:		" << parentCharge << std::endl;
	return parentCharge;
}

TVector3 SLDCorrection::cheatFlightDirection( MCParticle *SLDLepton )
{
	const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
	TVector3 flightDirection( MotherHadron->getMomentumAtEndpoint()[ 0 ] , MotherHadron->getMomentumAtEndpoint()[ 1 ] , MotherHadron->getMomentumAtEndpoint()[ 2 ] );
	streamlog_out(DEBUG0) << "		Parent Hadron Flight Direction	( x		, y		, z		)" << std::endl;
	flightDirection.SetMag( 1.0 );
	return flightDirection;
}

TVector3 SLDCorrection::getJetAxis( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	TVector3 jetAxis( 0.0 , 0.0 , 0.0 );
	ReconstructedParticle* linkedRecoLepton{};
	LCCollection *RecoJetCollection{};
	bool foundRecoLeptonInJet = false;
	try
	{
		linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , true , false );
		RecoJetCollection = pLCEvent->getCollection( m_inputJetCollection );
		streamlog_out(DEBUG2) << "	Found " << RecoJetCollection->getNumberOfElements() << " jet(s) in event " << m_nEvt << std::endl;
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "	Input jet collection not found in event " << m_nEvt << std::endl;
		return jetAxis;
	}

	linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , true , false );
	for ( int i_jet = 0 ; i_jet < RecoJetCollection->getNumberOfElements() ; ++i_jet )
	{
		ReconstructedParticle *recoJet = dynamic_cast<EVENT::ReconstructedParticle*>( RecoJetCollection->getElementAt( i_jet ) );
		ReconstructedParticleVec jetPFOs  = recoJet->getParticles();
		for ( long unsigned int i_pfo = 0 ; i_pfo < jetPFOs.size() ; ++i_pfo )
		{
			ReconstructedParticle *recoPFO = jetPFOs[ i_pfo ];
			if ( recoPFO == linkedRecoLepton )
			{
				TVector3 jetMomentum( recoJet->getMomentum()[ 0 ] , recoJet->getMomentum()[ 1 ] , recoJet->getMomentum()[ 2 ] );
				jetAxis = jetMomentum;
				foundRecoLeptonInJet = true;
			}
		}
	}
	if ( foundRecoLeptonInJet )
	{
		jetAxis.SetMag(1.0);
		return jetAxis;
	}
	else
	{
		streamlog_out(DEBUG2) << "	Couldn't find recoLepton in jets for assigning jet axis to flight direction of mother hadron" << m_nEvt << std::endl;
		return jetAxis;
	}
}

TVector3 SLDCorrection::getParentHadronUnitMomentum( TLorentzVector parentFourMomentum , std::vector< double > primaryVertex , std::vector< double > secondaryVertex , int parentCharge )
{
	TVector3 initialFlightDirection( secondaryVertex[ 0 ] - primaryVertex[ 0 ] , secondaryVertex[ 1 ] - primaryVertex[ 1 ] , secondaryVertex[ 2 ] - primaryVertex[ 2 ] );
	TVector3 parentPt( parentFourMomentum.Px() , parentFourMomentum.Py() , 0.0 );
	double pT = parentPt.Mag();
	double flightDistance2D = sqrt( pow( secondaryVertex[ 0 ] - primaryVertex[ 0 ] , 2 ) + pow( secondaryVertex[ 1 ] - primaryVertex[ 1 ] , 2 ) );
	double sinCorrectionAngle = parentCharge * eB * flightDistance2D / ( 2.0 * pT );
	double cosCorrectionAngle = sqrt( 1.0 - pow( sinCorrectionAngle , 2 ) );
	double correctFlightDirectionX = cosCorrectionAngle * initialFlightDirection.X() + sinCorrectionAngle * initialFlightDirection.Y();
	double correctFlightDirectionY = -sinCorrectionAngle * initialFlightDirection.X() + cosCorrectionAngle * initialFlightDirection.Y();
	TVector3 correctFlightDirection( correctFlightDirectionX , correctFlightDirectionY , initialFlightDirection.Z() );
	correctFlightDirection.SetMag( 1.0 );
	return correctFlightDirection;
}

double SLDCorrection::getParentHadronMass( MCParticle *SLDLepton )
{
	double TrueParentHadronMass = 0.0;
	try
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		TrueParentHadronMass = MotherHadron->getMass();
		streamlog_out(DEBUG4) << "		Parent Hadron (PDG: " << MotherHadron->getPDG() << "):" << std::endl;
		streamlog_out(DEBUG0) << "			Energy:									  " << MotherHadron->getEnergy() << std::endl;
		streamlog_out(DEBUG0) << "			Mass:		" << MotherHadron->getMass() << std::endl;
		streamlog_out(DEBUG0) << "			Momentum (px,py,pz): 		" << MotherHadron->getMomentum()[ 0 ] << "	, " << MotherHadron->getMomentum()[ 1 ] << "	, " << MotherHadron->getMomentum()[ 2 ] << std::endl;
		streamlog_out(DEBUG0) << "			MomentumEP (px,py,pz): 	" << MotherHadron->getMomentumAtEndpoint()[ 0 ] << "	, " << MotherHadron->getMomentumAtEndpoint()[ 1 ] << "	, " << MotherHadron->getMomentumAtEndpoint()[ 2 ] << std::endl;
		streamlog_out(DEBUG0) << "			Vertex (x,y,z): 		" << MotherHadron->getVertex()[ 0 ] << "	, " << MotherHadron->getVertex()[ 1 ] << "	, " << MotherHadron->getVertex()[ 2 ] << std::endl;
		streamlog_out(DEBUG0) << "			EndPoint (x,y,z): 		" << MotherHadron->getEndpoint()[ 0 ] << "	, " << MotherHadron->getEndpoint()[ 1 ] << "	, " << MotherHadron->getEndpoint()[ 2 ] << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Lepton for semi-leptonic decay not found" << std::endl;
        }
	return TrueParentHadronMass;
}

TLorentzVector SLDCorrection::getLeptonFourMomentum( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	TLorentzVector leptonFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	ReconstructedParticle* linkedRecoLepton = NULL;
	if ( m_cheatLepton4momentum )
	{
		leptonFourMomentum = TLorentzVector( SLDLepton->getMomentum()[ 0 ] , SLDLepton->getMomentum()[ 1 ] , SLDLepton->getMomentum()[ 2 ] , SLDLepton->getEnergy() );
	}
	else
	{
		if ( m_recoFourMomentumOfVisibles == 0 )
		{
			linkedRecoLepton = getLinkedTrack4MomOfPFO( pLCEvent , SLDLepton );
		}
		else if ( m_recoFourMomentumOfVisibles == 1 )
		{
			linkedRecoLepton = getLinkedChargedPFO( pLCEvent , SLDLepton );
		}
		else
		{
			linkedRecoLepton = getLinkedPFO( pLCEvent , SLDLepton , true , false );
		}

		if ( linkedRecoLepton != NULL )
		{
			leptonFourMomentum = TLorentzVector( linkedRecoLepton->getMomentum()[ 0 ] , linkedRecoLepton->getMomentum()[ 1 ] , linkedRecoLepton->getMomentum()[ 2 ] , linkedRecoLepton->getEnergy() );
		}
		else
		{
			leptonFourMomentum = TLorentzVector( SLDLepton->getMomentum()[ 0 ] , SLDLepton->getMomentum()[ 1 ] , SLDLepton->getMomentum()[ 2 ] , SLDLepton->getEnergy() );
		}
	}
	streamlog_out(DEBUG0) << "		Lepton" << std::endl;
	streamlog_out(DEBUG0) << "			True:(	" << SLDLepton->getPDG() << "	, " << SLDLepton->getMass() << "	, " << SLDLepton->getMomentum()[ 0 ] << "	, " << SLDLepton->getMomentum()[ 1 ] << "	, " << SLDLepton->getMomentum()[ 2 ] << "	, " << SLDLepton->getEnergy() << "	, " << SLDLepton->getCharge() << "	)" << std::endl;
	if ( !m_cheatLepton4momentum && linkedRecoLepton != NULL ) streamlog_out(DEBUG0) << "			Reco:(	" << linkedRecoLepton->getType() << "	, " << linkedRecoLepton->getMass() << "	, " << linkedRecoLepton->getMomentum()[ 0 ] << "	, " << linkedRecoLepton->getMomentum()[ 1 ] << "	, " << linkedRecoLepton->getMomentum()[ 2 ] << "	, " << linkedRecoLepton->getEnergy() << "	, " << linkedRecoLepton->getCharge() << "	)" << std::endl;
	return leptonFourMomentum;
}

TLorentzVector SLDCorrection::getVisibleFourMomentum( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton , MCParticle *parentHadron , bool getChargedTLV , bool getNeutralTLV )
{
	TLorentzVector VisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	try
	{
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( parentHadron->getDaughters() ).size() ; ++i_daughter )
		{
			EVENT::MCParticle *daughter = parentHadron->getDaughters()[ i_daughter ];
			ReconstructedParticle* linkedPFO = NULL;
			streamlog_out(DEBUG2) << "	Daughter[ " << i_daughter <<" ]: genStatus = " << daughter->getGeneratorStatus() << " , PDG = " << daughter->getPDG() << " , Charge = " << daughter->getCharge() << std::endl;
			if ( daughter->getGeneratorStatus() == 1 )
			{
				if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && daughter != SLDLepton )
				{
					if ( getChargedTLV && fabs( daughter->getCharge() ) >= 0.1 )
					{
						streamlog_out(DEBUG2) << "		Charged:" << std::endl;
						streamlog_out(DEBUG2) << "			True:(	" << daughter->getPDG() << "	, " << daughter->getMass() << "	, " << daughter->getMomentum()[ 0 ] << "	, " << daughter->getMomentum()[ 1 ] << "	, " << daughter->getMomentum()[ 2 ] << "	, " << daughter->getEnergy() << "	, " << daughter->getCharge() << "	)" << std::endl;
						if ( m_cheatCharged4momentum )// || linkedPFO == NULL )
						{
							VisibleFourMomentum += TLorentzVector( daughter->getMomentum()[ 0 ] , daughter->getMomentum()[ 1 ] , daughter->getMomentum()[ 2 ] , daughter->getEnergy() );
						}
						else
						{
							if ( m_recoFourMomentumOfVisibles == 0 )
							{
								linkedPFO = getLinkedTrack4MomOfPFO( pLCEvent , daughter );
							}
							else if ( m_recoFourMomentumOfVisibles == 1 )
							{
								linkedPFO = getLinkedChargedPFO( pLCEvent , daughter );
							}
							else
							{
								linkedPFO = getLinkedPFO( pLCEvent , SLDLepton , true , false );
							}
							if ( linkedPFO == NULL )
							{
//								VisibleFourMomentum += TLorentzVector( daughter->getMomentum()[ 0 ] , daughter->getMomentum()[ 1 ] , daughter->getMomentum()[ 2 ] , daughter->getEnergy() );
								++m_nChargedPFOwoTrack;
								h_MCPTracks->Fill( 1.5 );
								h_MCPTracks_Eweighted->Fill( 1.5 , daughter->getEnergy() );
								h_MCPTracks_Ptweighted->Fill( 1.5 , sqrt( pow( daughter->getMomentum()[ 0 ] , 2 ) + pow( daughter->getMomentum()[ 1 ] , 2 ) ) );
								m_lostChargedMCP_CosTheta.push_back( sqrt( pow( daughter->getMomentum()[ 0 ] , 2 ) + pow( daughter->getMomentum()[ 1 ] , 2 ) + pow( daughter->getMomentum()[ 2 ] , 2 ) ) / daughter->getMomentum()[ 2 ] );
								m_lostChargedMCP_Energy.push_back( daughter->getEnergy() );
								m_lostChargedMCP_Pt.push_back( sqrt( pow( daughter->getMomentum()[ 0 ] , 2 ) + pow( daughter->getMomentum()[ 1 ] , 2 ) ) );
							}
							else
							{
								h_MCPTracks->Fill( 0.5 );
								h_MCPTracks_Eweighted->Fill( 0.5 , daughter->getEnergy() );
								h_MCPTracks_Ptweighted->Fill( 0.5 , sqrt( pow( daughter->getMomentum()[ 0 ] , 2 ) + pow( daughter->getMomentum()[ 1 ] , 2 ) ) );
								VisibleFourMomentum += TLorentzVector( linkedPFO->getMomentum()[ 0 ] , linkedPFO->getMomentum()[ 1 ] , linkedPFO->getMomentum()[ 2 ] , linkedPFO->getEnergy() );
								streamlog_out(DEBUG2) << "			Reco:(	" << linkedPFO->getType() << "	, " << linkedPFO->getMass() << "	, " << linkedPFO->getMomentum()[ 0 ] << "	, " << linkedPFO->getMomentum()[ 1 ] << "	, " << linkedPFO->getMomentum()[ 2 ] << "	, " << linkedPFO->getEnergy() << "	, " << linkedPFO->getCharge() << "	)" << std::endl;
							}
						}
					}
					if ( getNeutralTLV && fabs( daughter->getCharge() ) < 0.1 )
					{
						VisibleFourMomentum += TLorentzVector( daughter->getMomentum()[ 0 ] , daughter->getMomentum()[ 1 ] , daughter->getMomentum()[ 2 ] , daughter->getEnergy() );
						streamlog_out(DEBUG0) << "		Neutral:" << std::endl;
						streamlog_out(DEBUG0) << "			True:(	" << daughter->getPDG() << "	, " << daughter->getMass() << "	, " << daughter->getMomentum()[ 0 ] << "	, " << daughter->getMomentum()[ 1 ] << "	, " << daughter->getMomentum()[ 2 ] << "	, " << daughter->getEnergy() << "	, " << daughter->getCharge() << "	)" << std::endl;
					}
				}
			}
			else
			{
				VisibleFourMomentum += getVisibleFourMomentum( pLCEvent , SLDLepton , daughter , getChargedTLV , getNeutralTLV );
			}
			if ( abs( daughter->getPDG() ) == 12 && abs( daughter->getPDG() ) == 14 && abs( daughter->getPDG() ) == 16 )
			{
				m_nNeutrino += 1;
			}

		}
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Lepton for semi-leptonic decay not found" << std::endl;
        }
	return VisibleFourMomentum;
}

TLorentzVector SLDCorrection::getTrueNeutrinoFourMomentum( MCParticle *SLDLepton )
{
	TLorentzVector InisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	try
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		int nNeutrinos = 0;
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( MotherHadron->getDaughters() ).size() ; ++i_daughter )
		{
			const EVENT::MCParticle *daughter = MotherHadron->getDaughters()[ i_daughter ];
			if ( daughter->getGeneratorStatus() == 1 && ( abs( daughter->getPDG() ) == abs( SLDLepton->getPDG() ) + 1 ) )
			{
				streamlog_out(DEBUG0) << "		Neutrino:" << std::endl;
				streamlog_out(DEBUG0) << "			True:(	" << daughter->getPDG() << "	, " << daughter->getMass() << "	, " << daughter->getMomentum()[ 0 ] << "	, " << daughter->getMomentum()[ 1 ] << "	, " << daughter->getMomentum()[ 2 ] << "	, " << daughter->getEnergy() << "	, " << daughter->getCharge() << "	)" << std::endl;
				InisibleFourMomentum = TLorentzVector( daughter->getMomentum()[ 0 ] , daughter->getMomentum()[ 1 ] , daughter->getMomentum()[ 2 ] , daughter->getEnergy() );
			}
		}
		++nNeutrinos;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	True Neutrino for semi-leptonic decay not found in MCParticles" << std::endl;
        }
	return InisibleFourMomentum;
}

std::vector<MCParticle*> SLDCorrection::getMCLeptonVec( EVENT::LCEvent *pLCEvent , int parentHadronPDG )
{
	LCCollection *MCParticleCollection{};
	std::vector<MCParticle*> mcLeptonVecSLD;
	std::string hadron = "";
	if ( parentHadronPDG == 4 ) hadron = "C-Hadron";
	if ( parentHadronPDG == 5 ) hadron = "B-Hadron";
	try
        {
		MCParticleCollection = pLCEvent->getCollection( m_mcParticleCollection );
		int nMCP = MCParticleCollection->getNumberOfElements();
		MCParticle *testLepton = NULL;
		for ( int i_mcp = 0 ; i_mcp < nMCP ; ++i_mcp )
		{
			testLepton = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_mcp ) );
			int nLeptonsInVertexSLD = 0;
			bool foundSLD = false;
			if ( ( abs( testLepton->getPDG() ) == 11 || abs( testLepton->getPDG() ) == 13 ) && ( testLepton->getGeneratorStatus() == 1 ) )
			{
				for ( long unsigned int i_parent = 0 ; i_parent < ( testLepton->getParents() ).size() ; ++i_parent )
				{
					const EVENT::MCParticle *testMotherHadron = testLepton->getParents()[ i_parent ];
					if ( floor( abs( testMotherHadron->getPDG() ) / 100 ) == parentHadronPDG || ( floor( abs( testMotherHadron->getPDG() ) / 1000 ) == parentHadronPDG ) )
					{
						for ( long unsigned int i_daughter = 0 ; i_daughter < ( testMotherHadron->getDaughters() ).size() ; ++i_daughter )
						{
							if ( ( abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 11 || abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 13 ) && ( testMotherHadron->getDaughters()[ i_daughter ] )->getGeneratorStatus() == 1 )
							{
								++nLeptonsInVertexSLD;
							}
							if ( abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == abs( testLepton->getPDG() ) + 1 )
							{
								mcLeptonVecSLD.push_back( testLepton );
								foundSLD = true;
							}
						}
					}
				}
				if ( foundSLD ) streamlog_out(DEBUG0) << "	One Semi-Leptonic Decay of " << hadron << " is found ; Number of leptons associated to SLD vertex : " << nLeptonsInVertexSLD << std::endl;
			}
		}
		streamlog_out(DEBUG0) << "	" << mcLeptonVecSLD.size() << " Semi-Leptonic Decay(s) of " << hadron << " found" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input " << MCParticleCollection << " collection not found in event " << pLCEvent->getEventNumber() << std::endl;
        }
	return mcLeptonVecSLD;
}

std::vector<MCParticle*> SLDCorrection::getMCVisibles( MCParticle *mcParentHadron , bool getChargedTLV , bool getNeutralTLV )
{
	std::vector<MCParticle*> mcVisibles;
	std::vector<MCParticle*> tempVisibles;
	MCParticle *daughter = NULL;
	for ( long unsigned int i_daughter = 0 ; i_daughter < ( mcParentHadron->getDaughters() ).size() ; ++i_daughter )
	{
		daughter = mcParentHadron->getDaughters()[ i_daughter ];
		if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && abs( daughter->getPDG() ) != 11 && abs( daughter->getPDG() ) != 13 && abs( daughter->getPDG() ) != 15 )
		{
			if ( getChargedTLV && fabs( daughter->getCharge() ) < 0.1 )
			{
				if ( daughter->getGeneratorStatus() == 1 )
				{
					mcVisibles.push_back( daughter );
					streamlog_out(DEBUG0) << "		Daughter[ "  << i_daughter << " ]" << std::endl;
					streamlog_out(DEBUG0) << "				" << daughter->getPDG() << "(" << daughter->getGeneratorStatus() << ")	, " << daughter->getMass() << "	, " << daughter->getMomentum()[ 0 ] << "	, " << daughter->getMomentum()[ 1 ] << "	, " << daughter->getMomentum()[ 2 ] << "	, " << daughter->getEnergy() << "	, " << daughter->getCharge() << "	)" << std::endl;
				}
				else
				{
					tempVisibles = getMCVisibles( daughter , getChargedTLV , getNeutralTLV );
				}
			}
			if ( getNeutralTLV && fabs( daughter->getCharge() ) > 0.1 )
			{
				if ( daughter->getGeneratorStatus() == 1 )
				{
					mcVisibles.push_back( daughter );
					streamlog_out(DEBUG0) << "		Daughter[ "  << i_daughter << " ]" << std::endl;
					streamlog_out(DEBUG0) << "				" << daughter->getPDG() << "(" << daughter->getGeneratorStatus() << ")	, " << daughter->getMass() << "	, " << daughter->getMomentum()[ 0 ] << "	, " << daughter->getMomentum()[ 1 ] << "	, " << daughter->getMomentum()[ 2 ] << "	, " << daughter->getEnergy() << "	, " << daughter->getCharge() << "	)" << std::endl;
				}
				else
				{
					tempVisibles = getMCVisibles( daughter , getChargedTLV , getNeutralTLV );
				}
			}
		}
	}
	for ( long unsigned int i_temp = 0 ; i_temp < tempVisibles.size() ; ++i_temp )
	{
		mcVisibles.push_back( tempVisibles.at( i_temp ) );
	}
	return mcVisibles;
}

TLorentzVector SLDCorrection::getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector FourMomentumLepton , TLorentzVector VisibleFourMomentumCharged , TLorentzVector VisibleFourMomentumNeutral , double ParentHadronMass , int solutionSign )
{
	int sign = ( solutionSign != 0 ? solutionSign / abs( solutionSign ) : 1 );
	const char *solSign = ( sign >= 0 ? "+" : "-" );
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "		--------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "		Calculate Neutrino 4-Momentum for " << solSign << " solution" << std::endl;
	streamlog_out(DEBUG1) << "		--------------------------------------------" << std::endl;

	flightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG1) << "		flightDirection:			( " << flightDirection.Px() << "	, " << flightDirection.Py() << "	, " << flightDirection.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG1) << "		Parent Hadron Mass =	 " << ParentHadronMass << std::endl;

	TLorentzVector visible_tlv	= VisibleFourMomentumCharged + VisibleFourMomentumNeutral + FourMomentumLepton;
	streamlog_out(DEBUG1) << "		Visible 4-Momentum:			( " << visible_tlv.Px() << "	, " << visible_tlv.Py() << "	, " << visible_tlv.Pz() << "	, " << visible_tlv.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Lepton			( " << FourMomentumLepton.Px() << "	, " << FourMomentumLepton.Py() << "	, " << FourMomentumLepton.Pz() << "	, " << FourMomentumLepton.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Charged		( " << VisibleFourMomentumCharged.Px() << "	, " << VisibleFourMomentumCharged.Py() << "	, " << VisibleFourMomentumCharged.Pz() << "	, " << VisibleFourMomentumCharged.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Neutral		( " << VisibleFourMomentumNeutral.Px() << "	, " << VisibleFourMomentumNeutral.Py() << "	, " << VisibleFourMomentumNeutral.Pz() << "	, " << VisibleFourMomentumNeutral.E() << " )" << std::endl;

	double visible_mass		= visible_tlv.M();
	streamlog_out(DEBUG1) << "		Visible Inv Mass:	" << visible_mass << std::endl;

	double visible_E		= visible_tlv.E();
	streamlog_out(DEBUG1) << "		Visible Energy:									" << visible_E << std::endl;

	TVector3 visible_p		= TVector3( visible_tlv.Px() , visible_tlv.Py() , visible_tlv.Pz() );
	streamlog_out(DEBUG1) << "		Visible Momentum:			( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG1) << "		(Visible Momentum).(FlightDirection):							" << visible_p.Dot( flightDirection ) << std::endl;

	TVector3 visible_p_par		= visible_p.Dot( flightDirection ) * flightDirection;
	streamlog_out(DEBUG1) << "		Visible Momentum (par):		( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	TVector3 visible_p_nor		= visible_p - visible_p_par;
	streamlog_out(DEBUG1) << "		Visible Momentum (nor):		( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;

	double visible_E_prime		= ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass );
	streamlog_out(DEBUG1) << "		Visible Energy (prime):								" << visible_E_prime << std::endl;

	TVector3 visible_p_par_prime	= sign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) * flightDirection;
	if ( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) < visible_p_nor.Mag2() )
	{
		visible_p_par_prime	= std::numeric_limits<double>::min() * flightDirection;
	}
	streamlog_out(DEBUG1) << "		Visible Momentum (par-prime):		( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;

	double parent_hadron_E		= ( ( visible_tlv.E() * ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) ) - visible_p.Dot( flightDirection ) * sign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) ) * ParentHadronMass / ( pow( visible_mass , 2 ) + visible_p_nor.Mag2() );
	TVector3 parent_hadron_p	= sqrt( pow( parent_hadron_E , 2 ) - pow( ParentHadronMass , 2 ) ) * flightDirection;
	streamlog_out(DEBUG2) << "		Parent Hadron Momentum:		( " << parent_hadron_p.Px() << "	, " << parent_hadron_p.Py() << "	, " << parent_hadron_p.Pz() << "	, " << parent_hadron_E << " )" << std::endl;

	double Neutrino_E		= parent_hadron_E - visible_E;
	streamlog_out(DEBUG1) << "		Neutrino Energy:									" << Neutrino_E << std::endl;

	TVector3 Neutrino_p_nor		= -1 * visible_p_nor;
	streamlog_out(DEBUG1) << "		Neutrino Momentum (nor):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p_par		= sqrt( pow( Neutrino_E , 2 ) - Neutrino_p_nor.Mag2() ) * flightDirection;
	streamlog_out(DEBUG1) << "		Neutrino Momentum (par):		( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p		= Neutrino_p_nor + Neutrino_p_par;
	streamlog_out(DEBUG1) << "		Neutrino Momentum:			( " << Neutrino_p_par.Px() << "	, " << Neutrino_p_par.Py() << "	, " << Neutrino_p_par.Pz() << "	)" << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino Mass:	" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino 4-Momentum:		( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
	return Neutrino_tlv;
}

std::vector< float > SLDCorrection::getNeutrinoCovMat( TVector3 flightDirection , TLorentzVector fourMomentumLepton , TLorentzVector visibleFourMomentumCharged , TLorentzVector visibleFourMomentumNeutral , double parentHadronMass , std::vector< float > flightDirectionCovMat , std::vector< float > LeptonCovMat , std::vector< float > visibleChargedCovMat , std::vector< float > visibleNeutralCovMat , double sigmaParentHadronMass )
{
	std::vector< float > NeutrinoCovMat( 10 , 0.0 );
	std::vector< float > ParentHadronCovMat( 10 , 0.0 );

	TLorentzVector visibleFourMomentum = fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral;
	double Mvis			= visibleFourMomentum.M();
	double Evis			= visibleFourMomentum.E();
	double EvisPrime		= ( pow( parentHadronMass , 2 ) + pow( Mvis , 2 ) ) / ( 2 * parentHadronMass );
	TVector3 visible_p		= visibleFourMomentum.Vect();
	TVector3 visible_p_par		= visible_p.Dot( flightDirection ) * flightDirection;
	TVector3 visible_p_nor		= visible_p - visible_p_par;
	TVector3 visible_p_par_prime	= sqrt( pow( ( pow( parentHadronMass , 2 ) - pow( Mvis , 2 ) ) / ( 2 * parentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) * flightDirection;
	double parentHadronEnergy	= ( ( Evis * EvisPrime ) - visible_p_par.Dot( visible_p_par_prime ) ) * parentHadronMass / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );
	std::vector< float > visibleCovMat( 10 , 0.0 );
	for ( int i_element = 0 ; i_element < 10 ; ++i_element )
	{
		visibleCovMat[ i_element ] = LeptonCovMat[ i_element ] + visibleChargedCovMat[ i_element ] + visibleNeutralCovMat[ i_element ];
	}

	double sigmaMvis		= getSigmaVisibleMass( visibleFourMomentum , visibleCovMat );
	double sigmaEvis		= std::sqrt( visibleCovMat[ 9 ] );
	double sigmaEvis_prime		= std::sqrt( 0.5 * pow( 1 - pow( Mvis / parentHadronMass , 2 ) , 2 ) * pow( sigmaParentHadronMass , 2 ) + pow( Mvis / parentHadronMass , 2 ) * pow( sigmaMvis , 2 ) );
	double sigmaPvisPar		= getSigmaPvisPar( flightDirection , visibleFourMomentum , flightDirectionCovMat , visibleCovMat );
	double sigmaPvisPar_prime	= getSigmaVisiblePpar_prime( flightDirection , visibleFourMomentum , parentHadronMass , flightDirectionCovMat , visibleCovMat , sigmaParentHadronMass , sigmaMvis );
	double sigmaPvisNor		= getSigmaPvisNor( visible_p , visible_p_par.Mag() , visibleCovMat , sigmaPvisPar );

	double dEparHad_dEvis			= parentHadronMass * EvisPrime / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );
	double dEparHad_dEvisPrime		= parentHadronMass * Evis / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );
	double dEparHad_dPvisPar		= -parentHadronMass * visible_p_par_prime.Mag() / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );
	double dEparHad_dPvisParPrime		= -parentHadronMass * visible_p_par.Mag() / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );
	double dEparHad_dPvisNor		= -2 * visible_p_nor.Mag() * parentHadronMass * ( Evis * EvisPrime - visible_p_par.Dot( visible_p_par_prime ) ) / pow( pow( Mvis , 2 ) + visible_p_nor.Mag2() , 2 );
	double dEparHad_dMvis			= -2 * Mvis * parentHadronMass * ( Evis * EvisPrime - visible_p_par.Dot( visible_p_par_prime ) ) / pow( pow( Mvis , 2 ) + visible_p_nor.Mag2() , 2 );
	double dEparHad_dMparentHadron		= ( Evis * EvisPrime - visible_p_par.Dot( visible_p_par_prime ) ) / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );

	double sigmaEparHad			= std::sqrt( 	pow( dEparHad_dEvis , 2 ) * pow( sigmaEvis , 2 ) + pow( dEparHad_dEvisPrime , 2 ) * pow( sigmaEvis_prime , 2 ) +
 								pow( dEparHad_dPvisPar , 2 ) * pow( sigmaPvisPar , 2 ) + pow( dEparHad_dPvisParPrime , 2 ) * pow( sigmaPvisPar_prime , 2 ) +
								pow( dEparHad_dMvis , 2 ) * pow( sigmaMvis , 2 ) + pow( dEparHad_dMparentHadron , 2 ) * pow( sigmaParentHadronMass , 2 ) +
								pow( dEparHad_dPvisNor , 2 ) * pow( sigmaPvisNor , 2 ) 	);
	ParentHadronCovMat		= getParentHadronCovMat( flightDirection , parentHadronEnergy , parentHadronMass , flightDirectionCovMat , sigmaEparHad );

	for ( int i_element = 0 ; i_element < 10 ; ++i_element )
	{
		NeutrinoCovMat[ i_element ] = ParentHadronCovMat[ i_element ] + visibleCovMat[ i_element ];
	}

	return NeutrinoCovMat;
}

double SLDCorrection::getSigmaVisibleMass( TLorentzVector visibleFourMomentum , std::vector< float > visibleCovMat )
{
	double sigma_Mvis	= 0.0;
	double Mvis		= visibleFourMomentum.M();
	double Evis		= visibleFourMomentum.E();
	double Pxvis		= visibleFourMomentum.Px();
	double Pyvis		= visibleFourMomentum.Py();
	double Pzvis		= visibleFourMomentum.Pz();
	double dM_dPx		= -Pxvis / Mvis;
	double dM_dPy		= -Pyvis / Mvis;
	double dM_dPz		= -Pzvis / Mvis;
	double dM_dE		= Evis / Mvis;
	double sigmaPx2		= visibleCovMat[ 0 ];
	double sigmaPxPy	= visibleCovMat[ 1 ];
	double sigmaPy2		= visibleCovMat[ 2 ];
	double sigmaPxPz	= visibleCovMat[ 3 ];
	double sigmaPyPz	= visibleCovMat[ 4 ];
	double sigmaPz2		= visibleCovMat[ 5 ];
	double sigmaPxE		= visibleCovMat[ 6 ];
	double sigmaPyE		= visibleCovMat[ 7 ];
	double sigmaPzE		= visibleCovMat[ 8 ];
	double sigmaE2		= visibleCovMat[ 9 ];
	sigma_Mvis		= std::sqrt( 	pow( dM_dPx , 2 ) * sigmaPx2 + pow( dM_dPy , 2 ) * sigmaPy2 + pow( dM_dPz , 2 ) * sigmaPz2 + pow( dM_dE , 2 ) * sigmaE2 +
	 					dM_dPx * dM_dPy * sigmaPxPy + dM_dPx * dM_dPz * sigmaPxPz + dM_dPx * dM_dE * sigmaPxE +
						dM_dPy * dM_dPx * sigmaPxPy + dM_dPy * dM_dPz * sigmaPyPz + dM_dPy * dM_dE * sigmaPyE +
						dM_dPz * dM_dPx * sigmaPxPz + dM_dPz * dM_dPy * sigmaPyPz + dM_dPz * dM_dE * sigmaPzE +
						dM_dE * dM_dPx * sigmaPxE + dM_dE * dM_dPy * sigmaPyE + dM_dE * dM_dPz * sigmaPzE );
	return sigma_Mvis;
}

double SLDCorrection::getSigmaVisiblePpar_prime( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , std::vector< float > flightDirectionCovMat , std::vector< float > visibleCovMat , double sigmaParentHadronMass , double sigma_Mvis )
{
	flightDirection.SetMag( 1.0 );
	double visible_mass		= visibleFourMomentum.M();
	TVector3 visible_p		= visibleFourMomentum.Vect();
	TVector3 visible_p_par		= visible_p.Dot( flightDirection ) * flightDirection;
	double visible_p_nor		= ( visible_p - visible_p_par ).Mag();
	double visible_p_par_prime	= sqrt( pow( ( pow( parentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * parentHadronMass ) , 2 ) - pow( visible_p_nor , 2 ) );
	if ( pow( ( pow( parentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * parentHadronMass ) , 2 ) < pow( visible_p_nor , 2 ) )
	{
		visible_p_par_prime	= std::numeric_limits<double>::min();
	}
	double sigmaPvisPar		= getSigmaPvisPar( flightDirection , visibleFourMomentum , flightDirectionCovMat , visibleCovMat );
	double sigmaPvisNor		= getSigmaPvisNor( visible_p , visible_p_par.Mag() , visibleCovMat , sigmaPvisPar );

	double dPvisParPrime_dMvis	= -visible_mass * ( 1 - pow( visible_mass / parentHadronMass , 2 ) ) / ( 2 * visible_p_par_prime );
	double dPvisParPrime_dMparent	= parentHadronMass * ( 1 - pow( visible_mass / parentHadronMass , 4 ) ) / ( 4 * visible_p_par_prime );
	double dPvisParPrime_dPvisNor	= -visible_p_nor / visible_p_par_prime;
	double sigmaPvisPar_prime	= std::sqrt( pow( dPvisParPrime_dMvis , 2 ) * pow( sigma_Mvis , 2 ) + pow( dPvisParPrime_dMparent , 2 ) * pow( sigmaParentHadronMass , 2 ) + pow( dPvisParPrime_dPvisNor , 2 ) * pow( sigmaPvisNor , 2 ) );
	return sigmaPvisPar_prime;
}

double SLDCorrection::getSigmaPvisPar( TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > flightDirectionCovMat , std::vector< float > visibleCovMat )
{
	flightDirection.SetMag( 1.0 );
	TVector3 visible_p		= visibleFourMomentum.Vect();

	double sigmaX2			= flightDirectionCovMat[ 0 ];
	double sigmaXY			= flightDirectionCovMat[ 1 ];
	double sigmaY2			= flightDirectionCovMat[ 2 ];
	double sigmaXZ			= flightDirectionCovMat[ 3 ];
	double sigmaYZ			= flightDirectionCovMat[ 4 ];
	double sigmaZ2			= flightDirectionCovMat[ 5 ];

	double sigmaPx2			= visibleCovMat[ 0 ];
	double sigmaPxPy		= visibleCovMat[ 1 ];
	double sigmaPy2			= visibleCovMat[ 2 ];
	double sigmaPxPz		= visibleCovMat[ 3 ];
	double sigmaPyPz		= visibleCovMat[ 4 ];
	double sigmaPz2			= visibleCovMat[ 5 ];

	double vertexX			= flightDirection.X();
	double vertexY			= flightDirection.Y();
	double vertexZ			= flightDirection.Z();
	double Px			= visibleFourMomentum.Px();
	double Py			= visibleFourMomentum.Py();
	double Pz			= visibleFourMomentum.Pz();

	double dPvisPar_dPx		= vertexX;
	double dPvisPar_dPy		= vertexY;
	double dPvisPar_dPz		= vertexZ;
	double dPvisPar_dX		= Px;
	double dPvisPar_dY		= Py;
	double dPvisPar_dZ		= Pz;

	double sigmaPvisPar	= std::sqrt( 	pow( dPvisPar_dPx , 2 ) * sigmaPx2 + pow( dPvisPar_dPy , 2 ) * sigmaPy2 + pow( dPvisPar_dPz , 2 ) * sigmaPz2 +
 						pow( dPvisPar_dX , 2 ) * sigmaX2 + pow( dPvisPar_dY , 2 ) * sigmaY2 + pow( dPvisPar_dZ , 2 ) * sigmaZ2 +
						2 * dPvisPar_dPx * dPvisPar_dPy * sigmaPxPy + 2 * dPvisPar_dPx * dPvisPar_dPz * sigmaPxPz + 2 * dPvisPar_dPy * dPvisPar_dPz * sigmaPyPz +
						2 * dPvisPar_dX * dPvisPar_dY * sigmaXY + 2 * dPvisPar_dX * dPvisPar_dZ * sigmaXZ + 2 * dPvisPar_dY * dPvisPar_dZ * sigmaYZ
					);
	return sigmaPvisPar;
}

double SLDCorrection::getSigmaPvisNor( TVector3 visibleMomentum , double visibleMomentumPar , std::vector< float > visibleCovMat , double sigmaPvisPar )
{
	double sigmaPx2			= visibleCovMat[ 0 ];
	double sigmaPxPy		= visibleCovMat[ 1 ];
	double sigmaPy2			= visibleCovMat[ 2 ];
	double sigmaPxPz		= visibleCovMat[ 3 ];
	double sigmaPyPz		= visibleCovMat[ 4 ];
	double sigmaPz2			= visibleCovMat[ 5 ];

	double visPx			= visibleMomentum.Px();
	double visPy			= visibleMomentum.Py();
	double visPz			= visibleMomentum.Pz();

	double PvisNor			= std::sqrt( visibleMomentum.Mag2() - pow( visibleMomentumPar , 2 ) );
	double dPvisNor_dPx		= visPx / PvisNor;
	double dPvisNor_dPy		= visPy / PvisNor;
	double dPvisNor_dPz		= visPz / PvisNor;
	double dPvisNor_dPvisPar	= -visibleMomentumPar / PvisNor;
	double sigmaPvisNor		= std::sqrt(	pow( dPvisNor_dPx , 2 ) * sigmaPx2 + pow( dPvisNor_dPy , 2 ) * sigmaPy2 + pow( dPvisNor_dPz , 2 ) * sigmaPz2 + pow( dPvisNor_dPvisPar , 2 ) * pow( sigmaPvisPar , 2 ) +
							2 * dPvisNor_dPx * dPvisNor_dPy * sigmaPxPy + 2 * dPvisNor_dPx * dPvisNor_dPz * sigmaPxPz + 2 * dPvisNor_dPy * dPvisNor_dPz * sigmaPyPz );
	return sigmaPvisNor;
}

std::vector<float> SLDCorrection::getParentHadronCovMat( TVector3 flightDirection , double parentHadronEnergy , double parentHadronMass , std::vector<float> flightDirectionCovMat , double parentHadronSigmaE )
{
	const int rows			= 4; // n rows jacobian
	const int columns		= 4; // n columns jacobian
	const int kspace_time_dim	= 4;

	TMatrixD covMatrixMomenta(kspace_time_dim,kspace_time_dim);
	std::vector<float> covP;

	float vertexX		=	flightDirection.X();
	float vertexY		=	flightDirection.Y();
	float vertexZ		=	flightDirection.Z();
	float vertexR		=	std::sqrt( pow( vertexX , 2 ) + pow( vertexY , 2 ) + pow( vertexZ , 2 ) );
	float vertexX2		=	pow( vertexX , 2 );
	float vertexY2		=	pow( vertexY , 2 );
	float vertexZ2		=	pow( vertexZ , 2 );
	float vertexR2		=	pow( vertexR , 2 );
	float vertexR3		=	pow( vertexR , 3 );
	float SigmaX2		=	flightDirectionCovMat[ 0 ];
	float SigmaXY		=	flightDirectionCovMat[ 1 ];
	float SigmaY2		=	flightDirectionCovMat[ 2 ];
	float SigmaXZ		=	flightDirectionCovMat[ 3 ];
	float SigmaYZ		=	flightDirectionCovMat[ 4 ];
	float SigmaZ2		=	flightDirectionCovMat[ 5 ];
	float SigmaE2		=	pow( parentHadronSigmaE , 2 );

	double parentHadronPmag	=	std::sqrt( pow( parentHadronEnergy , 2 ) - pow( parentHadronMass , 2 ) );

//	Define array with jacobian matrix elements by rows
	double jacobian_by_rows[rows*columns] =
	{
		parentHadronPmag * ( vertexR2 - vertexX2 ) / vertexR3		,	-parentHadronPmag * vertexX * vertexY / vertexR3		,	-parentHadronPmag * vertexX * vertexZ / vertexR3		,	0			,
		-parentHadronPmag * vertexY * vertexX / vertexR3		,	parentHadronPmag * ( vertexR2 - vertexY2 ) / vertexR3		,	-parentHadronPmag * vertexY * vertexZ / vertexR3		,	0			,
		-parentHadronPmag * vertexZ * vertexX / vertexR3		,	-parentHadronPmag * vertexZ * vertexY / vertexR3		,	parentHadronPmag * ( vertexR2 - vertexZ2 ) / vertexR3		,	0			,
		parentHadronEnergy * vertexX / ( parentHadronPmag * vertexR )	,	parentHadronEnergy * vertexY / ( parentHadronPmag * vertexR )	,	parentHadronEnergy * vertexZ / ( parentHadronPmag * vertexR )	,	1.0
	};

//	construct the Jacobian using previous array ("F" if filling by columns, "C" if filling by rows, $ROOTSYS/math/matrix/src/TMatrixT.cxx)
	TMatrixD jacobian(rows,columns, jacobian_by_rows, "C");
	streamlog_out(DEBUG0) << "	Jacobian array converted to Jacobian matrix" << std::endl;

//	cluster covariance matrix by rows
	double parentHadron_CovMat_by_rows[rows*rows] =
			{
				SigmaX2		,	SigmaXY		,	SigmaXZ		,	0	,
				SigmaXY		,	SigmaY2		,	SigmaYZ		,	0	,
				SigmaXZ		,	SigmaYZ		,	SigmaZ2		,	0	,
				0		,	0		,	0		,	SigmaE2
			};

	TMatrixD covMatrix_parentHadron(rows,rows, parentHadron_CovMat_by_rows, "C");

	covMatrixMomenta.Mult( TMatrixD( jacobian ,
					TMatrixD::kTransposeMult ,
					covMatrix_parentHadron) ,
					jacobian
					);

	covP.push_back( covMatrixMomenta(0,0) ); // x-x
	covP.push_back( covMatrixMomenta(1,0) ); // y-x
	covP.push_back( covMatrixMomenta(1,1) ); // y-y
	covP.push_back( covMatrixMomenta(2,0) ); // z-x
	covP.push_back( covMatrixMomenta(2,1) ); // z-y
	covP.push_back( covMatrixMomenta(2,2) ); // z-z
	covP.push_back( covMatrixMomenta(3,0) ); // e-x
	covP.push_back( covMatrixMomenta(3,1) ); // e-y
	covP.push_back( covMatrixMomenta(3,2) ); // e-z
	covP.push_back( covMatrixMomenta(3,3) ); // e-e
	streamlog_out(DEBUG0) << "	FourMomentumCovarianceMatrix for parent Hadron is Filled succesfully" << std::endl;

	return covP;

}

ReconstructedParticle* SLDCorrection::getLinkedPFO( EVENT::LCEvent *pLCEvent , MCParticle *visibleMCP , bool getChargedTLV , bool getNeutralTLV )
{
	streamlog_out(DEBUG0) << "	Look for PFO linked to MCParticle (PDG: " << visibleMCP->getPDG() << " , GenStat: " << visibleMCP->getGeneratorStatus() << ")" << std::endl;
	ReconstructedParticle* linkedPFO{};
	bool foundlinkedPFO = false;
	try
	{
		LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( m_RecoMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( m_MCTruthRecoLinkCollection ) );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator for PFO <-> MCParticle" << std::endl;
		return NULL;
	}
	LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( m_MCTruthRecoLinkCollection ) );
	const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( visibleMCP );
	const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( visibleMCP );
	streamlog_out(DEBUG0) << "	Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
	double maxweightPFOtoMCP = 0.;
	double maxweightMCPtoPFO = 0.;
	int iPFOtoMCPmax = -1;
	int iMCPtoPFOmax = -1;
	for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
	{
		double pfo_weight = 0.0;
		if ( getChargedTLV )
		{
			pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
		}
		else if ( getNeutralTLV )
		{
			pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
		}
		streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << int( PFOweightvec.at( i_pfo ) ) % 10000 / 1000.0 << " , Cluster: " << int( PFOweightvec.at( i_pfo ) ) / 10000 / 1000.0 << ")" << std::endl;
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
		if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightMCPtoPFO = pfo_weight;
			iMCPtoPFOmax = i_pfo;
			streamlog_out(DEBUG0) << "	PFO at index: " << i_pfo << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
		}
	}
	if ( iMCPtoPFOmax != -1 )
	{
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( iMCPtoPFOmax );
		LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( m_RecoMCTruthLinkCollection ) );
		const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( testPFO );
		const EVENT::FloatVec&  MCPweightvec = RecoMCParticleNav.getRelatedToWeights( testPFO );
		for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
		{
			double mcp_weight = 0.0;
			if ( getChargedTLV )
			{
				mcp_weight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
			}
			else if ( getNeutralTLV )
			{
				mcp_weight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
			}
			MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
			if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightPFOtoMCP = mcp_weight;
				iPFOtoMCPmax = i_mcp;
				streamlog_out(DEBUG0) << "	MCParticle at index: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
			}
		}
		if ( iPFOtoMCPmax != -1 )
		{
			if ( MCPvec.at( iPFOtoMCPmax ) == visibleMCP )
			{
				streamlog_out(DEBUG0) << "	Linked PFO to MCParticle found successfully " << std::endl;
				linkedPFO = testPFO;
				foundlinkedPFO = true;
			}
		}
	}

	if( foundlinkedPFO )
	{
		streamlog_out(DEBUG2) << "	Found linked RecoLepton (px,py,pz,E): ( " << linkedPFO->getMomentum()[ 0 ] << " 	, " << linkedPFO->getMomentum()[ 1 ] << " 	, " << linkedPFO->getMomentum()[ 2 ] << " 	, " << linkedPFO->getEnergy() << " 	)" << std::endl;
		return linkedPFO;
	}
	else
	{
		streamlog_out(DEBUG2) << "	Couldn't Find a PFO linked to MCParticle" << std::endl;
		return NULL;
	}
}

ReconstructedParticle* SLDCorrection::getLinkedChargedPFO( EVENT::LCEvent *pLCEvent , MCParticle *chargedMCP )
{
	streamlog_out(DEBUG0) << "	Look for charged PFO linked to MCParticle (" << chargedMCP->getPDG() << " , " << chargedMCP->getGeneratorStatus() << ")" << std::endl;
	LCCollection *inputPfoCollection{};
	ReconstructedParticle* linkedChargedPFO{};
	int i_linkedPFO = -1;
	Track *linkedTrack{};
	bool foundLinkedTrack = false;
	bool foundlinkedChargedPFO = false;
	int PFOTypes[6]{11,13,211,22,310,3122};
	try
	{
		LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator for Track <-> MCParticle" << std::endl;
		return NULL;
	}
	LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
	const EVENT::LCObjectVec& trkvec = MCParticleTrackNav.getRelatedToObjects( chargedMCP );
	const EVENT::FloatVec&  trkweightvec = MCParticleTrackNav.getRelatedToWeights( chargedMCP );
	streamlog_out(DEBUG0) << "	Charged visible MCParticle is linked to " << trkvec.size() << " tracks" << std::endl;
	double maxweightTRKtoMCP = 0.;
	double maxweightMCPtoTRK = 0.;
	int iTRKtoMCPmax = -1;
	int iMCPtoTRKmax = -1;
	for ( unsigned int i_trk = 0; i_trk < trkvec.size(); i_trk++ )
	{
		double track_weight = trkweightvec.at( i_trk );
		if ( track_weight > maxweightMCPtoTRK )//&& track_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightMCPtoTRK = track_weight;
			iMCPtoTRKmax = i_trk;
			streamlog_out(DEBUG0) << "	Track at index: " << i_trk << " has highest link weight to lepton = " << track_weight << std::endl;
		}
	}
	if ( iMCPtoTRKmax != -1 )
	{
		Track *testTrack = (Track *) trkvec.at( iMCPtoTRKmax );
		LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
		const EVENT::LCObjectVec& mcpvec = TrackMCParticleNav.getRelatedToObjects( testTrack );
		const EVENT::FloatVec&  mcpweightvec = TrackMCParticleNav.getRelatedToWeights( testTrack );
		for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
		{
			double mcp_weight = mcpweightvec.at(i_mcp);
			MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
			if ( mcp_weight > maxweightTRKtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightTRKtoMCP = mcp_weight;
				iTRKtoMCPmax = i_mcp;
				streamlog_out(DEBUG0) << "	MCParticle at index: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and MCParticle to Track link weight is " << mcp_weight << std::endl;
			}
		}
		if ( iTRKtoMCPmax != -1 )
		{
			if ( mcpvec.at( iTRKtoMCPmax ) == chargedMCP )
			{
				streamlog_out(DEBUG0) << "	Linked track to lepton found successfully " << std::endl;
				linkedTrack = testTrack;
				foundLinkedTrack = true;
			}
		}
	}
	if ( foundLinkedTrack )
	{
		inputPfoCollection = pLCEvent->getCollection( m_inputPfoCollection );
		int n_PFO = inputPfoCollection->getNumberOfElements();
		streamlog_out(DEBUG0) << "	Looking for linkedTrack in " << n_PFO << " PFOs" << std::endl;
		for (int i_pfo = 0; i_pfo < n_PFO ; ++i_pfo)
		{
			ReconstructedParticle* inputPFO = dynamic_cast<ReconstructedParticle*>( inputPfoCollection->getElementAt( i_pfo ) );
			const EVENT::TrackVec& inputPFOtrkvec = inputPFO->getTracks();
			for ( long unsigned int i_trk = 0 ; i_trk < inputPFOtrkvec.size() ; ++i_trk )
			{
				if ( linkedTrack == inputPFOtrkvec.at( i_trk ) )
				{
					streamlog_out(DEBUG0) << "	Found linkedPFOTrack" << std::endl;
					foundlinkedChargedPFO = true;
					i_linkedPFO = i_pfo;
				}
			}
		}
	}
	if( foundlinkedChargedPFO )
	{
		linkedChargedPFO = dynamic_cast<ReconstructedParticle*>( inputPfoCollection->getElementAt( i_linkedPFO ) );
		streamlog_out(DEBUG0) << "	Found linked RecoLepton (px,py,pz,E): ( " << linkedChargedPFO->getMomentum()[ 0 ] << " 	, " << linkedChargedPFO->getMomentum()[ 1 ] << " 	, " << linkedChargedPFO->getMomentum()[ 2 ] << " 	, " << linkedChargedPFO->getEnergy() << " 	)" << std::endl;
		bool knownPFO = false;
		for ( int l = 0 ; l < 6 ; ++l )
		{
			if ( abs( linkedChargedPFO->getType() ) == PFOTypes[ l ] )
			{
				knownPFO = true;
				if ( abs( chargedMCP->getPDG() ) == 11 )
				{
					h_recoPFOLinkedToElectron_Type->Fill( l );
				}
				else if ( abs( chargedMCP->getPDG() ) == 13 )
				{
					h_recoPFOLinkedToMuon_Type->Fill( l );
				}
			}
		}
		if ( !knownPFO )
		{
			if ( abs( chargedMCP->getPDG() ) == 11 )
			{
				h_recoPFOLinkedToElectron_Type->Fill( 7 );
			}
			else if ( abs( chargedMCP->getPDG() ) == 13 )
			{
				h_recoPFOLinkedToMuon_Type->Fill( 7 );
			}
		}
		return linkedChargedPFO;
	}
	else
	{
		streamlog_out(DEBUG0) << "	Couldn't find linked RecoLepton" << std::endl;
		if ( abs( chargedMCP->getPDG() ) == 11 )
		{
			h_recoPFOLinkedToElectron_Type->Fill( 8 );
		}
		else if ( abs( chargedMCP->getPDG() ) == 13 )
		{
			h_recoPFOLinkedToMuon_Type->Fill( 8 );
		}
		return NULL;
	}
}

ReconstructedParticleImpl* SLDCorrection::getLinkedTrack4MomOfPFO( EVENT::LCEvent *pLCEvent , MCParticle *chargedMCP )
{
	streamlog_out(DEBUG0) << "	Look for charged PFO linked to MCParticle (" << chargedMCP->getPDG() << " , " << chargedMCP->getGeneratorStatus() << ")" << std::endl;
	ReconstructedParticleImpl* linkedChargedPFO = new ReconstructedParticleImpl;
	TLorentzVector chargedPFOFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	Track *linkedTrack{};
	bool foundLinkedTrack = false;
	try
	{
		LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator for Track <-> MCParticle" << std::endl;
		return NULL;
	}
	LCRelationNavigator MCParticleTrackNav( pLCEvent->getCollection( m_MCTruthTrackLinkCollection ) );
	const EVENT::LCObjectVec& trkvec = MCParticleTrackNav.getRelatedToObjects( chargedMCP );
	const EVENT::FloatVec&  trkweightvec = MCParticleTrackNav.getRelatedToWeights( chargedMCP );
	streamlog_out(DEBUG0) << "	Charged visible MCParticle is linked to " << trkvec.size() << " tracks" << std::endl;
	double maxweightTRKtoMCP = 0.;
	double maxweightMCPtoTRK = 0.;
	int iTRKtoMCPmax = -1;
	int iMCPtoTRKmax = -1;
	for ( unsigned int i_trk = 0; i_trk < trkvec.size(); i_trk++ )
	{
		double track_weight = trkweightvec.at( i_trk );
		if ( track_weight > maxweightMCPtoTRK )//&& track_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightMCPtoTRK = track_weight;
			iMCPtoTRKmax = i_trk;
			streamlog_out(DEBUG0) << "	Track at index: " << i_trk << " has highest link weight to lepton = " << track_weight << std::endl;
		}
	}
	if ( iMCPtoTRKmax != -1 )
	{
		Track *testTrack = (Track *) trkvec.at( iMCPtoTRKmax );
		LCRelationNavigator TrackMCParticleNav( pLCEvent->getCollection( m_TrackMCTruthLinkCollection ) );
		const EVENT::LCObjectVec& mcpvec = TrackMCParticleNav.getRelatedToObjects( testTrack );
		const EVENT::FloatVec&  mcpweightvec = TrackMCParticleNav.getRelatedToWeights( testTrack );
		for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
		{
			double mcp_weight = mcpweightvec.at(i_mcp);
			MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
			if ( mcp_weight > maxweightTRKtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightTRKtoMCP = mcp_weight;
				iTRKtoMCPmax = i_mcp;
				streamlog_out(DEBUG0) << "	MCParticle at index: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and MCParticle to Track link weight is " << mcp_weight << std::endl;
			}
		}
		if ( iTRKtoMCPmax != -1 )
		{
			if ( mcpvec.at( iTRKtoMCPmax ) == chargedMCP )
			{
				streamlog_out(DEBUG0) << "	Linked track to lepton found successfully " << std::endl;
				linkedTrack = testTrack;
				foundLinkedTrack = true;
			}
		}
	}
	if ( foundLinkedTrack )
	{
		chargedPFOFourMomentum = getTrackFourMomentum( linkedTrack , 0.13957018 );
		streamlog_out(DEBUG0) << "	Track parametrs vonverted to 4-momentum of charged visible " << std::endl;
		double Momentum[3]{ chargedPFOFourMomentum.Px() , chargedPFOFourMomentum.Py() , chargedPFOFourMomentum.Pz() };
		double Energy = chargedPFOFourMomentum.E();
		double Mass = chargedPFOFourMomentum.M();
		streamlog_out(DEBUG0) << "	4-momentum of charged visible is formed" << std::endl;
		linkedChargedPFO->setMomentum( Momentum );
		linkedChargedPFO->setEnergy( Energy );
		linkedChargedPFO->setMass( Mass );
		streamlog_out(DEBUG0) << "	Linked charged visible is formed " << std::endl;
		return linkedChargedPFO;
	}
	else
	{
		return NULL;
	}

}

TLorentzVector SLDCorrection::getTrackFourMomentum( EVENT::Track* inputTrk , double trackMass )
{
	streamlog_out(DEBUG1) << "	------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "	Calculating PFO 4-momentum from track parameters" << std::endl;
	streamlog_out(DEBUG1) << "	------------------------------------------------" << std::endl;
	double Phi = inputTrk->getPhi();
	double Omega = inputTrk->getOmega();
	double tanLambda = inputTrk->getTanLambda();
	streamlog_out(DEBUG0) << "	Track parameters obtained" << std::endl;
	double pT = eB / fabs( Omega );
	double px = pT * TMath::Cos( Phi );
	double py = pT * TMath::Sin( Phi );
	double pz = pT * tanLambda;
	double E = sqrt( pow( trackMass , 2 ) + px * px + py * py + pz * pz);
	streamlog_out(DEBUG0) << "	Track parameters is converted to (p,E)" << std::endl;
	TLorentzVector trackFourMomentum( px , py , pz , E );
	return trackFourMomentum;
}

void SLDCorrection::plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat )
{
	double NuPxResidual = FourMomentumNuClose.Px() - trueFourMomentumNeutrino.Px(); m_NuPxResidual.push_back( NuPxResidual );
	double NuPyResidual = FourMomentumNuClose.Py() - trueFourMomentumNeutrino.Py(); m_NuPyResidual.push_back( NuPyResidual );
	double NuPzResidual = FourMomentumNuClose.Pz() - trueFourMomentumNeutrino.Pz(); m_NuPzResidual.push_back( NuPzResidual );
	double NuEResidual = FourMomentumNuClose.E() - trueFourMomentumNeutrino.E(); m_NuEResidual.push_back( NuEResidual );
	double NuPxNormalizedResidual = NuPxResidual / sqrt( NeutrinoCovMat[ 0 ] ); m_NuPxNormalizedResidual.push_back( NuPxNormalizedResidual );
	double NuPyNormalizedResidual = NuPyResidual / sqrt( NeutrinoCovMat[ 2 ] ); m_NuPyNormalizedResidual.push_back( NuPyNormalizedResidual );
	double NuPzNormalizedResidual = NuPzResidual / sqrt( NeutrinoCovMat[ 5 ] ); m_NuPzNormalizedResidual.push_back( NuPzNormalizedResidual );
	double NuENormalizedResidual = NuEResidual / sqrt( NeutrinoCovMat[ 9 ] ); m_NuENormalizedResidual.push_back( NuENormalizedResidual );
	h_NuPxResidual->Fill( NuPxResidual ); ++n_NuPxResidual;
	h_NuPyResidual->Fill( NuPyResidual ); ++n_NuPyResidual;
	h_NuPzResidual->Fill( NuPzResidual ); ++n_NuPzResidual;
	h_NuEResidual->Fill( NuEResidual ); ++n_NuEResidual;
	h_NuPxNormalizedResidual->Fill( NuPxNormalizedResidual ); ++n_NuPxNormalizedResidual;
	h_NuPyNormalizedResidual->Fill( NuPyNormalizedResidual ); ++n_NuPyNormalizedResidual;
	h_NuPzNormalizedResidual->Fill( NuPzNormalizedResidual ); ++n_NuPzNormalizedResidual;
	h_NuENormalizedResidual->Fill( NuENormalizedResidual ); ++n_NuENormalizedResidual;
	h_recoNuPx_mcNuPx->Fill( trueFourMomentumNeutrino.Px() , FourMomentumNuClose.Px() );
	h_recoNuPy_mcNuPy->Fill( trueFourMomentumNeutrino.Py() , FourMomentumNuClose.Py() );
	h_recoNuPz_mcNuPz->Fill( trueFourMomentumNeutrino.Pz() , FourMomentumNuClose.Pz() );
	h_recoNuE_mcNuE->Fill( trueFourMomentumNeutrino.E() , FourMomentumNuClose.E() );
}

void SLDCorrection::InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle )
{
	histogram->Scale( 1.0 / scale );
	histogram->SetLineColor( color );
	histogram->SetLineWidth( lineWidth );
	histogram->SetMarkerSize( markerSize );
	histogram->SetMarkerStyle( markerStyle );
	histogram->SetMarkerColor( color );
	float fit_range = 4.0;
	float fit_min = -2.0;
	float fit_max = 2.0;
	doProperGaussianFit( histogram , fit_min , fit_max , fit_range );
	histogram->GetFunction("gaus")->SetLineColor( color );
	float y_max = 1.2 * histogram->GetMaximum();
	histogram->GetYaxis()->SetRangeUser(0.0, y_max);
	histogram->Write();
}

void SLDCorrection::doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange )
{
	float Chi2 = 0.0;
	float NDF = 0.0;
	for ( int i_fit = 0 ; i_fit < 3 ; ++i_fit )
	{
		histogram->Fit( "gaus" , "" , "" , fitMin , fitMax );
		TF1 *fitFunction = (TF1 *)histogram->GetFunction("gaus");
		double fitMean = fitFunction->GetParameter( 1 );
		double fitSigma = fitFunction->GetParameter( 2 );
		fitMin = fitMean - fitRange * fitSigma;
		fitMax = fitMean + fitRange * fitSigma;
		Chi2 = fitFunction->GetChisquare();
		NDF = fitFunction->GetNDF();
	}
	streamlog_out(DEBUG2) << "	FIT : CHI2(" << Chi2 << ") / NDF(" << NDF << ") = " << Chi2 / NDF << " 	, fitrange = " << fitRange << std::endl;
	streamlog_out(DEBUG2) << "" << std::endl;
	if ( Chi2 != 0.0 && NDF != 0.0 && Chi2 / NDF > 2.0 && fitRange >= 0.5 )
	{
		doProperGaussianFit( histogram , fitMin , fitMax , fitRange - 0.1 );
	}
}

void SLDCorrection::check( EVENT::LCEvent *pLCEvent )
{
	LCCollection *inputPfoCollection{};
	LCCollection *SLDNuCollection{};
	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		SLDNuCollection = pLCEvent->getCollection(m_SLDNuCollection);
		int n_inputPFOs = inputPfoCollection->getNumberOfElements();
		int n_outputPFOs = SLDNuCollection->getNumberOfElements();
		streamlog_out(DEBUG4) << "	CHECK : processed events: " << m_nEvtSum << " (Number of inputPFOS: " << n_inputPFOs << " , Number of outputPFOs: " << n_outputPFOs <<")" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "	Input/Output collection not found in event " << m_nEvt << std::endl;
        }
}

void SLDCorrection::end()
{

	if ( m_fillRootTree )
	{
		m_pTFile->cd();
		m_pTTree->Write();
		InitializeHistogram( h_NuPxResidual , n_NuPxResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPyResidual , n_NuPyResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuPzResidual , n_NuPzResidual , 4 , 1 , 1.0 , 1 );
		InitializeHistogram( h_NuEResidual , n_NuEResidual , 4 , 1 , 1.0 , 1 );
		if( !m_cheatLepton4momentum )
		{
/*
			InitializeHistogram( h_NuPxNormalizedResidual , n_NuPxNormalizedResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuPyNormalizedResidual , n_NuPyNormalizedResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuPzNormalizedResidual , n_NuPzNormalizedResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuENormalizedResidual , n_NuENormalizedResidual , 4 , 1 , 1.0 , 1 );
*/
		}
		h_recoNuPx_mcNuPx->Write();
		h_recoNuPy_mcNuPy->Write();
		h_recoNuPz_mcNuPz->Write();
		h_recoNuE_mcNuE->Write();
		h_parentPx_daughtersPx->Write();
		h_parentPy_daughtersPy->Write();
		h_parentPz_daughtersPz->Write();
		h_parentE_daughtersE->Write();
		h_recoPFOLinkedToElectron_Type->Write();
		h_recoPFOLinkedToMuon_Type->Write();
		h_SLDecayOrder->Write();
		h_foundVertex->Write();
		h_secondaryVertex->Write();
		h_parentHadronCharge->Write();
		h_MCPTracks->Write();
		h_MCPTracks_Eweighted->Write();
		h_MCPTracks_Ptweighted->Write();

		m_pTFile->Close();
	}
	delete m_pTFile;

	streamlog_out(MESSAGE) << " " << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	processed events: 	" << m_nEvtSum << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << " " << std::endl;

}
