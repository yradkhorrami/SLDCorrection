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

	registerProcessorParameter(	"cheatLepton4momentum",
					"Cheat FourMomentum of lepton in semi-leptonic decays",
					m_cheatLepton4momentum,
					bool(true)
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
	printParameters();
	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");
	m_pTTree = new TTree("SLDCorrection", "SLDCorrection");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("NuPxResidual", &m_NuPxResidual);
	m_pTTree->Branch("NuPyResidual", &m_NuPyResidual);
	m_pTTree->Branch("NuPzResidual", &m_NuPzResidual);
	m_pTTree->Branch("NuEResidual", &m_NuEResidual);
	m_pTTree->Branch("NuPxNormalizedResidual", &m_NuPxNormalizedResidual);
	m_pTTree->Branch("NuPyNormalizedResidual", &m_NuPyNormalizedResidual);
	m_pTTree->Branch("NuPzNormalizedResidual", &m_NuPzNormalizedResidual);
	m_pTTree->Branch("NuENormalizedResidual", &m_NuENormalizedResidual);
	h_NuPxResidual = new TH1F( "" , "; _{}p_{x,#nu}^{REC} - p_{x,#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPxResidual = 0;
	h_NuPyResidual = new TH1F( "" , "; _{}p_{y,#nu}^{REC} - p_{y,#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPyResidual = 0;
	h_NuPzResidual = new TH1F( "" , "; _{}p_{z,#nu}^{REC} - p_{z,#nu}^{MC}  [GeV]; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPzResidual = 0;
	h_NuEResidual = new TH1F( "" , "; _{}E_{#nu}^{REC} - E_{#nu}^{MC} [GeV]; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuEResidual = 0;
	h_NuPxNormalizedResidual = new TH1F( "" , "; ( _{}p_{x,#nu}^{REC} - p_{x,#nu}^{MC} ) / #sigma_{p_{x,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPxNormalizedResidual = 0;
	h_NuPyNormalizedResidual = new TH1F( "" , "; ( _{}p_{y,#nu}^{REC} - p_{y,#nu}^{MC} ) / #sigma_{p_{y,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPyNormalizedResidual = 0;
	h_NuPzNormalizedResidual = new TH1F( "" , "; ( _{}p_{z,#nu}^{REC} - p_{z,#nu}^{MC} ) / #sigma_{p_{z,#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuPzNormalizedResidual = 0;
	h_NuENormalizedResidual = new TH1F( "" , "; ( _{}E_{#nu}^{REC} - E_{#nu}^{MC} ) / #sigma_{E_{#nu}}; Normalized Entries / 0.1" , 200 , -10.0 , 10.0 ); n_NuENormalizedResidual = 0;
	h_recoNuPx_mcNuPx = new TH2F( "" , "; _{}p_{x,#nu}^{MC} [GeV] ;  _{}p_{x,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
	h_recoNuPy_mcNuPy = new TH2F( "" , "; _{}p_{y,#nu}^{MC} [GeV] ;  _{}p_{y,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
	h_recoNuPz_mcNuPz = new TH2F( "" , "; _{}p_{z,#nu}^{MC} [GeV] ;  _{}p_{z,#nu}^{REC} [GeV]" , 120 , -60.0 , 60.0 , 120 , -60.0 , 60.0 );
	h_recoNuE_mcNuE = new TH2F( "" , "; _{}E_{#nu}^{MC} [GeV] ;  _{}E_{#nu}^{REC} [GeV]" , 100 , 0.0 , 100.0 , 100 , 0.0 , 100.0 );
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
}

void SLDCorrection::Clear()
{
	m_NuPxResidual.clear();
	m_NuPyResidual.clear();
	m_NuPzResidual.clear();
	m_NuEResidual.clear();
	m_NuPxNormalizedResidual.clear();
	m_NuPyNormalizedResidual.clear();
	m_NuPzNormalizedResidual.clear();
	m_NuENormalizedResidual.clear();
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
	LCCollection *RecoJetCollection{};
	++m_nEvtSum;
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

        try
        {
		std::vector<MCParticle*> mcLeptonVecBSLD;
		std::vector<MCParticle*> mcLeptonVecCSLD;
		mcLeptonVecBSLD = getMCLeptonVec( pLCEvent , 5 );
		mcLeptonVecCSLD = getMCLeptonVec( pLCEvent , 4 );
		streamlog_out(DEBUG4) << "	" << mcLeptonVecBSLD.size() << " semi-leptonic decay of B-Hadron found in MCTruth " << std::endl;
		for ( long unsigned int i_BSLD = 0 ; i_BSLD < mcLeptonVecBSLD.size() ; ++i_BSLD )
		{
			MCParticle *mcLepton = mcLeptonVecBSLD.at( i_BSLD );
			streamlog_out(DEBUG3) << "		" << i_BSLD+1 << ") PDG code of Lepton: 	" << mcLepton->getPDG() << std::endl;
			streamlog_out(DEBUG3) << "		" << i_BSLD+1 << ") Momentum (px,py,pz): 	( " << mcLepton->getMomentum()[ 0 ] << "	, " << mcLepton->getMomentum()[ 1 ] << "	 , " << mcLepton->getMomentum()[ 2 ] << "	)" << std::endl;
			streamlog_out(DEBUG3) << "		" << i_BSLD+1 << ") Vertex (x,y,z): 		( " << mcLepton->getVertex()[ 0 ] << "	, " << mcLepton->getVertex()[ 1 ] << "	 , " << mcLepton->getVertex()[ 2 ] << "	)" << std::endl;
		}
		streamlog_out(DEBUG4) << "	" << mcLeptonVecCSLD.size() << " semi-leptonic decay of C-Hadron found in MCTruth " << std::endl;
		for ( long unsigned int i_CSLD = 0 ; i_CSLD < mcLeptonVecCSLD.size() ; ++i_CSLD )
		{
			MCParticle *mcLepton = mcLeptonVecCSLD.at( i_CSLD );
			streamlog_out(DEBUG3) << "		" << i_CSLD+1 << ") PDG code of Lepton: 	" << mcLepton->getPDG() << std::endl;
			streamlog_out(DEBUG3) << "		" << i_CSLD+1 << ") Momentum (px,py,pz): 	( " << mcLepton->getMomentum()[ 0 ] << "	, " << mcLepton->getMomentum()[ 1 ] << "	 , " << mcLepton->getMomentum()[ 2 ] << "	)" << std::endl;
			streamlog_out(DEBUG3) << "		" << i_CSLD+1 << ") Vertex (x,y,z): 		( " << mcLepton->getVertex()[ 0 ] << "	, " << mcLepton->getVertex()[ 1 ] << "	 , " << mcLepton->getVertex()[ 2 ] << "	)" << std::endl;
		}
		std::vector<MCParticle*> mcSLDLeptons;
		for ( long unsigned int i_BSLD = 0 ; i_BSLD < mcLeptonVecBSLD.size() ; ++i_BSLD )
		{
			mcSLDLeptons.push_back( mcLeptonVecBSLD.at( i_BSLD ) );
		}
		for ( long unsigned int i_CSLD = 0 ; i_CSLD < mcLeptonVecCSLD.size() ; ++i_CSLD )
		{
			mcSLDLeptons.push_back( mcLeptonVecCSLD.at( i_CSLD ) );
		}
		std::vector<TVector3> trueFlightDirections;
		std::vector<TLorentzVector> trueFourMomentumLeptons;
		std::vector<TLorentzVector> trueVisibleFourMomentumChargeds;
		std::vector<TLorentzVector> trueVisibleFourMomentumNeutrals;
		std::vector<TVector3> recoFlightDirections;
		std::vector<TLorentzVector> recoFourMomentumLeptons;
		std::vector<TLorentzVector> recoVisibleFourMomentumChargeds;
		std::vector<TLorentzVector> recoVisibleFourMomentumNeutrals;
		std::vector< std::vector< float > > recoFlightDirectionsCovMat;
		std::vector< std::vector< float > > recoLeptonsCovMat;
		std::vector< std::vector< float > > recoVisibleChargedsCovMat;
		std::vector< std::vector< float > > recoVisibleNeutralsCovMat;
		std::vector< bool > foundLinkedRecoLeptons;
		std::vector<TLorentzVector> FourMomentumNeutrinosPos;
		std::vector<TLorentzVector> FourMomentumNeutrinosNeg;
		std::vector<TLorentzVector> FourMomentumNeutrinosClose;
		std::vector<double> trueParentHadronMasses;
		std::vector<TLorentzVector> trueFourMomentumNeutrinos;
		std::vector< std::vector< float > > NeutrinosCovMat;
		for ( long unsigned int i_SLD = 0 ; i_SLD < mcSLDLeptons.size() ; ++i_SLD )
		{
			ReconstructedParticle* linkedRecoLepton = getLinkedRecoLepton( pLCEvent , mcSLDLeptons.at( i_SLD ) );
			std::vector< float > recoLeptonCovMat( 10 , 0.0 );
			TLorentzVector recoFourMomentumLepton( 0.0 , 0.0 , 0.0 , 0.0 );
			bool foundLinkedRecoLepton = false;
			if ( linkedRecoLepton != NULL )
			{
				recoFourMomentumLepton = TLorentzVector( linkedRecoLepton->getMomentum()[ 0 ] , linkedRecoLepton->getMomentum()[ 1 ] , linkedRecoLepton->getMomentum()[ 2 ] , linkedRecoLepton->getEnergy() );
				recoLeptonCovMat = linkedRecoLepton->getCovMatrix();
				foundLinkedRecoLepton = true;
			}
			else
			{
				recoFourMomentumLepton = TLorentzVector( ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 0 ] , ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 1 ] , ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 2 ] , ( mcSLDLeptons.at( i_SLD ) )->getEnergy() );
				foundLinkedRecoLepton = false;
			}
			recoFourMomentumLeptons.push_back( recoFourMomentumLepton );
			recoLeptonsCovMat.push_back( recoLeptonCovMat );
			foundLinkedRecoLeptons.push_back( foundLinkedRecoLepton );

			TVector3 trueFlightDirection = getTrueFlightDirection( mcSLDLeptons.at( i_SLD ) );
			trueFlightDirections.push_back( trueFlightDirection );
			TLorentzVector trueFourMomentumLepton = TLorentzVector( ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 0 ] , ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 1 ] , ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 2 ] , ( mcSLDLeptons.at( i_SLD ) )->getEnergy() );
			trueFourMomentumLeptons.push_back( trueFourMomentumLepton );
			TLorentzVector trueVisibleFourMomentumCharged = getTrueVisibleFourMomentum( mcSLDLeptons.at( i_SLD ) , true , false );
			trueVisibleFourMomentumChargeds.push_back( trueVisibleFourMomentumCharged );
			TLorentzVector trueVisibleFourMomentumNeutral = getTrueVisibleFourMomentum( mcSLDLeptons.at( i_SLD ) , false , true );
			trueVisibleFourMomentumNeutrals.push_back( trueVisibleFourMomentumNeutral );
			TLorentzVector trueFourMomentumNeutrino = getTrueFourMomentumNeutrino( mcSLDLeptons.at( i_SLD ) );
			trueFourMomentumNeutrinos.push_back( trueFourMomentumNeutrino );
			TLorentzVector trueVisibleFourMomentum = trueFourMomentumLepton + trueVisibleFourMomentumCharged + trueVisibleFourMomentumNeutral;
			streamlog_out(DEBUG3) << "		              ( PDG	, Mass		, Px		, Py		, Pz		, E		)" << std::endl;
			streamlog_out(DEBUG3) << "		Lepton 4-Mom"  <<  std::endl;
			streamlog_out(DEBUG3) << "				" << ( mcSLDLeptons.at( i_SLD ) )->getPDG() << "	  " << trueFourMomentumLepton.Mag() << "	, " << trueFourMomentumLepton.Px() << "	, " << trueFourMomentumLepton.Py() << "	, " << trueFourMomentumLepton.Pz() << "	, " << trueFourMomentumLepton.E() << "	)" << std::endl;
			streamlog_out(DEBUG3) << "		Visible 4-Mom"  <<  std::endl;
			streamlog_out(DEBUG3) << "					  " << trueVisibleFourMomentum.Mag() << "	, " << trueVisibleFourMomentum.Px() << "	, " << trueVisibleFourMomentum.Py() << "	, " << trueVisibleFourMomentum.Pz() << "	, " << trueVisibleFourMomentum.E() << "	)" << std::endl;
			streamlog_out(DEBUG3) << "		Invisible 4-Mom"  <<  std::endl;
			streamlog_out(DEBUG3) << "					  " << trueFourMomentumNeutrino.Mag() << "	, " << trueFourMomentumNeutrino.Px() << "	, " << trueFourMomentumNeutrino.Py() << "	, " << trueFourMomentumNeutrino.Pz() << "	, " << trueFourMomentumNeutrino.E() << "	)" << std::endl;
			TLorentzVector trueFourMomentum = trueVisibleFourMomentum + trueFourMomentumNeutrino;
			streamlog_out(DEBUG3) << "		Total 4-Mom"  <<  std::endl;
			streamlog_out(DEBUG3) << "					  " << trueFourMomentum.Mag() << "	, " << trueFourMomentum.Px() << "	, " << trueFourMomentum.Py() << "	, " << trueFourMomentum.Pz() << "	, " << trueFourMomentum.E() << "	)" << std::endl;
			trueParentHadronMasses.push_back( getTrueParentHadronMass( mcSLDLeptons.at( i_SLD ) ) );
		}
		int nSLD = mcSLDLeptons.size();
		for ( int i_SLD = 0 ; i_SLD < nSLD ; ++i_SLD )
		{
			double sigmaParentHadronMass = 0.0;
			std::vector< float > flightDirectionCovMat( 6 , 0.0 );
			std::vector< float > LeptonCovMat( 10 , 0.0 );
			std::vector< float > visibleChargedCovMat( 10 , 0.0 );
			std::vector< float > visibleNeutralCovMat( 10 , 0.0 );
			std::vector< float > NeutrinoCovMat( 10 , 0.0 );

			TLorentzVector fourMomentumLepton;// = trueFourMomentumLeptons.at( i_SLD );
			if ( m_cheatLepton4momentum || !foundLinkedRecoLeptons.at( i_SLD ) )
			{
				fourMomentumLepton = trueFourMomentumLeptons.at( i_SLD );
			}
			else
			{
				fourMomentumLepton = recoFourMomentumLeptons.at( i_SLD );
				LeptonCovMat = recoLeptonsCovMat.at( i_SLD );
			}

			TVector3 flightDirection = trueFlightDirections.at( i_SLD );
			TLorentzVector visibleFourMomentumCharged = trueVisibleFourMomentumChargeds.at( i_SLD );
			TLorentzVector visibleFourMomentumNeutral = trueVisibleFourMomentumNeutrals.at( i_SLD );
			double parentHadronMass = trueParentHadronMasses.at( i_SLD );
			TLorentzVector FourMomentumNuPos;
			TLorentzVector FourMomentumNuNeg;
			TLorentzVector FourMomentumNuClose;
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "	Calculate Neutrino 4-Momentum ( + solution )" << std::endl;
			FourMomentumNuPos = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , +1 );
			FourMomentumNeutrinosPos.push_back( FourMomentumNuPos );
			streamlog_out(DEBUG4) << "	----------------------------------------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "	Calculate Neutrino 4-Momentum ( - solution )" << std::endl;
			FourMomentumNuNeg = getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , -1 );
			FourMomentumNeutrinosNeg.push_back( FourMomentumNuNeg );
			FourMomentumNuClose = ( abs( FourMomentumNuPos.E() - ( trueFourMomentumNeutrinos.at( i_SLD ) ).E() ) < abs( FourMomentumNuNeg.E() - ( trueFourMomentumNeutrinos.at( i_SLD ) ).E() ) ? FourMomentumNuPos : FourMomentumNuNeg );
			FourMomentumNeutrinosClose.push_back( FourMomentumNuClose );
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "	Closest Neutrino 4-Momentum:		( " << FourMomentumNuClose.Px() << "	, " << FourMomentumNuClose.Py() << "	, " << FourMomentumNuClose.Pz() << "	, " << FourMomentumNuClose.E() << " )" << std::endl;
			streamlog_out(DEBUG4) << "	True Neutrino 4-Momentum:		( " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Px() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Py() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Pz() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).E() << " )" << std::endl;
			if ( abs( FourMomentumNuClose.E() - ( trueFourMomentumNeutrinos.at( i_SLD ) ).E() ) > 1.0 )
			{
				streamlog_out(DEBUG2) << "	!!! Big Difference between true and reco neutrino Energy : " << FourMomentumNuClose.E() - ( trueFourMomentumNeutrinos.at( i_SLD ) ).E() << "  GeV" << std::endl;
			}
			streamlog_out(DEBUG4) << "	----------------------------------------------------------------------------------------------------" << std::endl;
			NeutrinoCovMat = getNeutrinoCovMat( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , flightDirectionCovMat , LeptonCovMat , visibleChargedCovMat , visibleNeutralCovMat , sigmaParentHadronMass );
			streamlog_out(DEBUG4) << "	Reconstructed Neutrino CovMat:" << std::endl;
			streamlog_out(DEBUG4) << "						" << NeutrinoCovMat[ 0 ] << std::endl;
			streamlog_out(DEBUG4) << "						" << NeutrinoCovMat[ 1 ] << "	, " << NeutrinoCovMat[ 2 ] << std::endl;
			streamlog_out(DEBUG4) << "						" << NeutrinoCovMat[ 3 ] << "	, " << NeutrinoCovMat[ 4 ]  << "	, " << NeutrinoCovMat[ 5 ] << std::endl;
			streamlog_out(DEBUG4) << "						" << NeutrinoCovMat[ 6 ] << "	, " << NeutrinoCovMat[ 7 ]  << "	, " << NeutrinoCovMat[ 8 ] << "	, " << NeutrinoCovMat[ 9 ]  << std::endl;
			streamlog_out(DEBUG4) << "	----------------------------------------------------------------------------------------------------" << std::endl;
			NeutrinosCovMat.push_back( NeutrinoCovMat );
			plotHistograms( trueFourMomentumNeutrinos.at( i_SLD ) , FourMomentumNuClose , NeutrinoCovMat );
		}

		RecoJetCollection = pLCEvent->getCollection( m_inputJetCollection );
/*
		int nJets = RecoJetCollection->getNumberOfElements();
		int nPFOs = 0;

		for ( int i_jet = 0 ; i_jet < nJets ; ++i_jet )
		{
			ReconstructedParticle* Jet = dynamic_cast<ReconstructedParticle*>( RecoJetCollection->getElementAt( i_jet ) );
			ReconstructedParticleVec jetPFOs = Jet->getParticles();
       		 	nPFOs = jetPFOs.size();
       		 	streamlog_out(DEBUG4) << "	Number of PFOs in the jet[ " << i_jet << " ]: " << nPFOs << std::endl;
			for ( int i_pfo = 0 ; i_pfo < nPFOs ; ++i_pfo )
			{
				ReconstructedParticle *testPFO = jetPFOs[ i_pfo ];
				if ( abs( testPFO->getType() ) == 11 || abs( testPFO->getType() ) == 13 )
				{
					const EVENT::TrackVec& PFOtrkvec = testPFO->getTracks();
					Track *inputTrk = (Track*)PFOtrkvec.at(0);
					Vertex *pfoVertex = testPFO->getStartVertex();
//					inputTrk->getRadiusOfInnermostHit();
					streamlog_out(DEBUG4) << "	Found a lepton in the jet[ " << i_jet << " ]: " << " with TYPE: " << testPFO->getType() << std::endl;
					streamlog_out(DEBUG3) << "		TYPE: 			" << testPFO->getType() << std::endl;
					if ( PFOtrkvec.size() != 1 )
					{
						streamlog_out(DEBUG6) << "		*************** Number of tracks:	" << PFOtrkvec.size() << " ***************" << std::endl;
					}
					streamlog_out(DEBUG3) << "		Number of tracks:	" << PFOtrkvec.size() << std::endl;
					streamlog_out(DEBUG3) << "		RadiusOfInnermostHit:	" << inputTrk->getRadiusOfInnermostHit() << std::endl;
					streamlog_out(DEBUG3) << "		Momentum (px,py,pz): 	( " << testPFO->getMomentum()[ 0 ] << "	, " << testPFO->getMomentum()[ 1 ] << "	 , " << testPFO->getMomentum()[ 2 ] << "	)" << std::endl;
//					streamlog_out(DEBUG3) << "		Vertex (x,y,z): 	( " << ( testPFO->getStartVertex() )->getPosition()[ 0 ] << "	, " << ( testPFO->getStartVertex() )->getPosition()[ 1 ] << "	 , " << ( testPFO->getStartVertex() )->getPosition()[ 2 ] << "	)" << std::endl;
				}
			}
		}
*/
		m_pTTree->Fill();
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
        }

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

TVector3 SLDCorrection::getTrueFlightDirection( MCParticle *SLDLepton )
{
	TVector3 flightDirection( 0.0 , 0.0 , 0.0 );
	try
	{
		int nParentHadrons = ( SLDLepton->getParents() ).size();
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		flightDirection = TVector3( MotherHadron->getMomentum()[ 0 ] , MotherHadron->getMomentum()[ 1 ] , MotherHadron->getMomentum()[ 2 ] );
		flightDirection.SetMag( 1.0 );
		streamlog_out(DEBUG0) << "	" << nParentHadrons << " True Parent Hadron(s) for semi-leptonic decay found in MCParticles" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Lepton for semi-leptonic decay not found" << std::endl;
        }
	return flightDirection;
}

double SLDCorrection::getTrueParentHadronMass( MCParticle *SLDLepton )
{
	double TrueParentHadronMass = 0.0;
	try
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		TrueParentHadronMass = MotherHadron->getMass();
		streamlog_out(DEBUG3) << "		Parent Hadron" << std::endl;
		streamlog_out(DEBUG3) << "				" << MotherHadron->getPDG() << "	, " << MotherHadron->getMass() << "	, " << MotherHadron->getMomentum()[ 0 ] << "	, " << MotherHadron->getMomentum()[ 1 ] << "	, " << MotherHadron->getMomentum()[ 2 ] << "	, " << MotherHadron->getEnergy() << "	)" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Lepton for semi-leptonic decay not found" << std::endl;
        }
	return TrueParentHadronMass;
}

TLorentzVector SLDCorrection::getTrueVisibleFourMomentum( MCParticle *SLDLepton , bool getChargedTLV , bool getNeutralTLV )
{
	TLorentzVector VisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	try
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		streamlog_out(DEBUG0) << "		|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||" << std::endl;
		streamlog_out(DEBUG0) << "		Parent Hadron ( PDG	, Mass		, Px		, Py		, Pz		, E		)" << std::endl;
		streamlog_out(DEBUG0) << "				" << MotherHadron->getPDG() << "	, " << MotherHadron->getMass() << "	, " << MotherHadron->getMomentum()[ 0 ] << "	, " << MotherHadron->getMomentum()[ 1 ] << "	, " << MotherHadron->getMomentum()[ 2 ] << "	, " << MotherHadron->getEnergy() << "	)" << std::endl;
		streamlog_out(DEBUG0) << "		              ( PDG	, Mass		, Px		, Py		, Pz		, E		, Charge	)" << std::endl;
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( MotherHadron->getDaughters() ).size() ; ++i_daughter )
		{
			const EVENT::MCParticle *daughter = MotherHadron->getDaughters()[ i_daughter ];
			if ( abs( daughter->getPDG() ) != 12 && abs( daughter->getPDG() ) != 14 && abs( daughter->getPDG() ) != 16 && abs( daughter->getPDG() ) != 11 && abs( daughter->getPDG() ) != 13 && abs( daughter->getPDG() ) != 15 )
			{
				if ( getChargedTLV && ( daughter->getCharge() != 0 ) )
				{
					VisibleFourMomentum += TLorentzVector( daughter->getMomentum()[ 0 ] , daughter->getMomentum()[ 1 ] , daughter->getMomentum()[ 2 ] , daughter->getEnergy() );
					streamlog_out(DEBUG0) << "		Daughter[ "  << i_daughter << " ]" << std::endl;
					streamlog_out(DEBUG0) << "				" << daughter->getPDG() << "	, " << daughter->getMass() << "	, " << daughter->getMomentum()[ 0 ] << "	, " << daughter->getMomentum()[ 1 ] << "	, " << daughter->getMomentum()[ 2 ] << "	, " << daughter->getEnergy() << "	, " << daughter->getCharge() << "	)" << std::endl;
				}
				if ( getNeutralTLV && ( daughter->getCharge() == 0 ) )
				{
					VisibleFourMomentum += TLorentzVector( daughter->getMomentum()[ 0 ] , daughter->getMomentum()[ 1 ] , daughter->getMomentum()[ 2 ] , daughter->getEnergy() );
					streamlog_out(DEBUG0) << "		Daughter[ "  << i_daughter << " ]" << std::endl;
					streamlog_out(DEBUG0) << "				" << daughter->getPDG() << "	, " << daughter->getMass() << "	, " << daughter->getMomentum()[ 0 ] << "	, " << daughter->getMomentum()[ 1 ] << "	, " << daughter->getMomentum()[ 2 ] << "	, " << daughter->getEnergy() << "	, " << daughter->getCharge() << "	)" << std::endl;
				}
			}

		}
		streamlog_out(DEBUG0) << "		Visible 4-Mom ( PDG	, Mass		, Px		, Py		, Pz		, E	)" << std::endl;
		streamlog_out(DEBUG0) << "				" << MotherHadron->getPDG() << "	, " << VisibleFourMomentum.Mag() << "	, " << VisibleFourMomentum.Px() << "	, " << VisibleFourMomentum.Py() << "	, " << VisibleFourMomentum.Pz() << "	, " << VisibleFourMomentum.E() << "	)" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Lepton for semi-leptonic decay not found" << std::endl;
        }
	return VisibleFourMomentum;
}

TLorentzVector SLDCorrection::getTrueFourMomentumNeutrino( MCParticle *SLDLepton )
{
	TLorentzVector InisibleFourMomentum( 0.0 , 0.0 , 0.0 , 0.0 );
	try
	{
		const EVENT::MCParticle *MotherHadron = SLDLepton->getParents()[ 0 ];
		int nNeutrinos = 0;
		for ( long unsigned int i_daughter = 0 ; i_daughter < ( MotherHadron->getDaughters() ).size() ; ++i_daughter )
		{
			if ( ( MotherHadron->getDaughters()[ i_daughter ] )->getGeneratorStatus() == 1 && ( abs( ( MotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 12 || abs( ( MotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 14 || abs( ( MotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 16 ) )
			{
				InisibleFourMomentum = TLorentzVector( ( MotherHadron->getDaughters()[ i_daughter ] )->getMomentum()[ 0 ] , ( MotherHadron->getDaughters()[ i_daughter ] )->getMomentum()[ 1 ] , ( MotherHadron->getDaughters()[ i_daughter ] )->getMomentum()[ 2 ] , ( MotherHadron->getDaughters()[ i_daughter ] )->getEnergy() );
			}
		}
		++nNeutrinos;
		streamlog_out(DEBUG0) << "	" << nNeutrinos << " True Neutrino(s) for semi-leptonic decay found in MCParticles" << std::endl;
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	True Neutrino for semi-leptonic decay not found in MCParticles" << std::endl;
        }
	return InisibleFourMomentum;
}

TLorentzVector SLDCorrection::getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector FourMomentumLepton , TLorentzVector VisibleFourMomentumCharged , TLorentzVector VisibleFourMomentumNeutral , double ParentHadronMass , int solutionSign )
{
	int sign = ( solutionSign != 0 ? solutionSign / abs( solutionSign ) : 1 );
	const char *solSign = ( sign >= 0 ? "+" : "-" );
	streamlog_out(DEBUG1) << "		Calculate Neutrino 4-Momentum for " << solSign << " solution" << std::endl;

	flightDirection.SetMag( 1.0 );
	streamlog_out(DEBUG1) << "		flightDirection:		( " << flightDirection.Px() << "	, " << flightDirection.Py() << "	, " << flightDirection.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG1) << "		Parent Hadron Mass =	 	" << ParentHadronMass << std::endl;

	TLorentzVector visible_tlv	= VisibleFourMomentumCharged + VisibleFourMomentumNeutral + FourMomentumLepton;
	streamlog_out(DEBUG1) << "		Visible 4-Momentum:		( " << visible_tlv.Px() << "	, " << visible_tlv.Py() << "	, " << visible_tlv.Pz() << "	, " << visible_tlv.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Lepton		( " << FourMomentumLepton.Px() << "	, " << FourMomentumLepton.Py() << "	, " << FourMomentumLepton.Pz() << "	, " << FourMomentumLepton.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Charged	( " << VisibleFourMomentumCharged.Px() << "	, " << VisibleFourMomentumCharged.Py() << "	, " << VisibleFourMomentumCharged.Pz() << "	, " << VisibleFourMomentumCharged.E() << " )" << std::endl;
	streamlog_out(DEBUG1) << "				Neutral	( " << VisibleFourMomentumNeutral.Px() << "	, " << VisibleFourMomentumNeutral.Py() << "	, " << VisibleFourMomentumNeutral.Pz() << "	, " << VisibleFourMomentumNeutral.E() << " )" << std::endl;

	double visible_mass		= visible_tlv.Mag();
	streamlog_out(DEBUG1) << "		Visible Invariant Mass:	" << visible_mass << std::endl;

	double visible_E		= visible_tlv.E();
	streamlog_out(DEBUG1) << "		Visible Energy:		" << visible_E << std::endl;

	TVector3 visible_p		= visible_tlv.Vect();
	streamlog_out(DEBUG1) << "		Visible Momentum:		( " << visible_p.Px() << "	, " << visible_p.Py() << "	, " << visible_p.Pz() << "	)" << std::endl;
	streamlog_out(DEBUG1) << "		(Visible Momentum).(FlightDirection):		" << visible_p.Dot( flightDirection ) << std::endl;

	TVector3 visible_p_par		= visible_p.Dot( flightDirection ) * flightDirection;
	streamlog_out(DEBUG1) << "		Visible Momentum (par):	( " << visible_p_par.Px() << "	, " << visible_p_par.Py() << "	, " << visible_p_par.Pz() << "	)" << std::endl;

	TVector3 visible_p_nor		= visible_p - visible_p_par;
	streamlog_out(DEBUG1) << "		Visible Momentum (nor):	( " << visible_p_nor.Px() << "	, " << visible_p_nor.Py() << "	, " << visible_p_nor.Pz() << "	)" << std::endl;

	double visible_E_prime		= ( pow( ParentHadronMass , 2 ) + pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass );
	streamlog_out(DEBUG1) << "		Visible Energy (prime):	" << visible_E_prime << std::endl;

	TVector3 visible_p_par_prime	= sign * sqrt( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) * flightDirection;
	if ( pow( ( pow( ParentHadronMass , 2 ) - pow( visible_mass , 2 ) ) / ( 2 * ParentHadronMass ) , 2 ) < visible_p_nor.Mag2() )
	{
		visible_p_par_prime	= std::numeric_limits<double>::min() * flightDirection;
	}
	streamlog_out(DEBUG1) << "		Visible Momentum (par-prime):	( " << visible_p_par_prime.Px() << "	, " << visible_p_par_prime.Py() << "	, " << visible_p_par_prime.Pz() << "	)" << std::endl;

	double parent_hadron_E		= ( ( visible_E * visible_E_prime ) - visible_p_par.Dot( visible_p_par_prime ) ) * ParentHadronMass / ( pow( visible_mass , 2 ) + visible_p_nor.Mag2() );
	streamlog_out(DEBUG1) << "		Parent Hadron Energy:		" << parent_hadron_E << std::endl;

	TVector3 parent_hadron_p	= sqrt( pow( parent_hadron_E , 2 ) - pow( ParentHadronMass , 2 ) ) * flightDirection;
	streamlog_out(DEBUG1) << "		Parent Hadron Momentum:	( " << parent_hadron_p.Px() << "	, " << parent_hadron_p.Py() << "	, " << parent_hadron_p.Pz() << "	)" << std::endl;

	double Neutrino_E		= parent_hadron_E - visible_E;
	streamlog_out(DEBUG1) << "		Neutrino Energy:		" << Neutrino_E << std::endl;

	TVector3 Neutrino_p_nor		= -1 * visible_p_nor;
	streamlog_out(DEBUG1) << "		Neutrino Momentum (nor):	( " << Neutrino_p_nor.Px() << "	, " << Neutrino_p_nor.Py() << "	, " << Neutrino_p_nor.Pz() << "	)" << std::endl;

	TVector3 Neutrino_p		= parent_hadron_p - visible_p;
	streamlog_out(DEBUG1) << "		Neutrino Momentum:		( " << Neutrino_p.Px() << "	, " << Neutrino_p.Py() << "	, " << Neutrino_p.Pz() << "	)" << std::endl;

	TLorentzVector Neutrino_tlv( Neutrino_p , Neutrino_E );
	streamlog_out(DEBUG2) << "" << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino Mass:		" << Neutrino_tlv.Mag() << std::endl;
	streamlog_out(DEBUG2) << "	Reconstructed Neutrino 4-Momentum:	( " << Neutrino_tlv.Px() << "	, " << Neutrino_tlv.Py() << "	, " << Neutrino_tlv.Pz() << "	, " << Neutrino_tlv.E() << " )" << std::endl;
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

ReconstructedParticle* SLDCorrection::getLinkedRecoLepton( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton )
{
	LCCollection *inputPfoCollection{};
	ReconstructedParticle* linkedRecoLepton{};
	int i_linkedPFO = -1;
	Track *linkedTrack{};
	bool foundLinkedTrack = false;
	bool foundLinkedRecoLepton = false;
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
	const EVENT::LCObjectVec& trkvec = MCParticleTrackNav.getRelatedToObjects( SLDLepton );
	const EVENT::FloatVec&  trkweightvec = MCParticleTrackNav.getRelatedToWeights( SLDLepton );
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
			if ( mcpvec.at( iTRKtoMCPmax ) == SLDLepton )
			{
				streamlog_out(DEBUG0) << "	Linked track to lepton found successfully " << std::endl;
				linkedTrack = testTrack;
				foundLinkedTrack = true;
			}
		}
	}
	if ( foundLinkedTrack )
	{
		streamlog_out(DEBUG0) << "	linkedTrack: " << linkedTrack << std::endl;
		inputPfoCollection = pLCEvent->getCollection( m_inputPfoCollection );
		int n_PFO = inputPfoCollection->getNumberOfElements();
		streamlog_out(DEBUG0) << "	Looking for linkedTrack in " << n_PFO << " PFOs" << std::endl;
		for (int i_pfo = 0; i_pfo < n_PFO ; ++i_pfo)
		{
			ReconstructedParticle* inputPFO = dynamic_cast<ReconstructedParticle*>( inputPfoCollection->getElementAt( i_pfo ) );
			const EVENT::TrackVec& inputPFOtrkvec = inputPFO->getTracks();
			for ( long unsigned int i_trk = 0 ; i_trk < inputPFOtrkvec.size() ; ++i_trk )
			{
				streamlog_out(DEBUG0) << "	pfoTrack: " << inputPFOtrkvec.at( i_trk ) << std::endl;
				if ( linkedTrack == inputPFOtrkvec.at( i_trk ) )
				{
					streamlog_out(DEBUG0) << "	Found linkedPFOTrack" << std::endl;
					foundLinkedRecoLepton = true;
					i_linkedPFO = i_pfo;
//					continue;
//					linkedRecoLepton = inputPfoCollection->getElementAt( i_pfo );
				}
			}
		}
	}
	if( foundLinkedRecoLepton )
	{
		linkedRecoLepton = dynamic_cast<ReconstructedParticle*>( inputPfoCollection->getElementAt( i_linkedPFO ) );
		streamlog_out(DEBUG0) << "	Found linked RecoLepton (px,py,pz,E): ( " << linkedRecoLepton->getMomentum()[ 0 ] << " 	, " << linkedRecoLepton->getMomentum()[ 1 ] << " 	, " << linkedRecoLepton->getMomentum()[ 2 ] << " 	, " << linkedRecoLepton->getEnergy() << " 	)" << std::endl;
		bool knownPFO = false;
		for ( int l = 0 ; l < 6 ; ++l )
		{
			if ( abs( linkedRecoLepton->getType() ) == PFOTypes[ l ] )
			{
				knownPFO = true;
				if ( abs( SLDLepton->getPDG() ) == 11 )
				{
					h_recoPFOLinkedToElectron_Type->Fill( l );
				}
				else if ( abs( SLDLepton->getPDG() ) == 13 )
				{
					h_recoPFOLinkedToMuon_Type->Fill( l );
				}
			}
		}
		if ( !knownPFO )
		{
			if ( abs( SLDLepton->getPDG() ) == 11 )
			{
				h_recoPFOLinkedToElectron_Type->Fill( 7 );
			}
			else if ( abs( SLDLepton->getPDG() ) == 13 )
			{
				h_recoPFOLinkedToMuon_Type->Fill( 7 );
			}
		}
		return linkedRecoLepton;
	}
	else
	{
		streamlog_out(DEBUG0) << "	Couldn't find linked RecoLepton" << std::endl;
		if ( abs( SLDLepton->getPDG() ) == 11 )
		{
			h_recoPFOLinkedToElectron_Type->Fill( 8 );
		}
		else if ( abs( SLDLepton->getPDG() ) == 13 )
		{
			h_recoPFOLinkedToMuon_Type->Fill( 8 );
		}
		return NULL;
	}
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
	float fit_range = 2.0;
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
	streamlog_out(DEBUG4) << "	FIT : CHI2(" << Chi2 << ") / NDF(" << NDF << ") = " << Chi2 / NDF << " 	, fitrange = " << fitRange << std::endl;
	streamlog_out(DEBUG4) << "" << std::endl;
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
		if( !m_cheatLepton4momentum )
		{
			InitializeHistogram( h_NuPxResidual , n_NuPxResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuPyResidual , n_NuPyResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuPzResidual , n_NuPzResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuEResidual , n_NuEResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuPxNormalizedResidual , n_NuPxNormalizedResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuPyNormalizedResidual , n_NuPyNormalizedResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuPzNormalizedResidual , n_NuPzNormalizedResidual , 4 , 1 , 1.0 , 1 );
			InitializeHistogram( h_NuENormalizedResidual , n_NuENormalizedResidual , 4 , 1 , 1.0 , 1 );
		}
		h_recoNuPx_mcNuPx->Write();
		h_recoNuPy_mcNuPy->Write();
		h_recoNuPz_mcNuPz->Write();
		h_recoNuE_mcNuE->Write();
		h_recoPFOLinkedToElectron_Type->Write();
		h_recoPFOLinkedToMuon_Type->Write();
		m_pTFile->Close();
		delete m_pTFile;
	}

	streamlog_out(MESSAGE) << " " << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	processed events: 	" << m_nEvtSum << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << " " << std::endl;

}
