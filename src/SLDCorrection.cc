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
	m_nEvtSum(0)
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
	m_pTTree = new TTree("PFOswithRFT", "PFOswithRFT");
	m_pTTree->SetDirectory(m_pTFile);
}

void SLDCorrection::Clear()
{

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
		std::vector<TLorentzVector> FourMomentumNeutrinosPos;
		std::vector<TLorentzVector> FourMomentumNeutrinosNeg;
		std::vector<TLorentzVector> FourMomentumNeutrinosClose;
		std::vector<double> trueParentHadronMasses;
		std::vector<TLorentzVector> trueFourMomentumNeutrinos;
		for ( long unsigned int i_SLD = 0 ; i_SLD < mcSLDLeptons.size() ; ++i_SLD )
		{
			trueFourMomentumLeptons.push_back( TLorentzVector( ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 0 ] , ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 1 ] , ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 2 ] , ( mcSLDLeptons.at( i_SLD ) )->getEnergy() ) );
			TLorentzVector trueFourMomentumLepton = TLorentzVector( ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 0 ] , ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 1 ] , ( mcSLDLeptons.at( i_SLD ) )->getMomentum()[ 2 ] , ( mcSLDLeptons.at( i_SLD ) )->getEnergy() );
			trueFlightDirections.push_back( getTrueFlightDirection( mcSLDLeptons.at( i_SLD ) ) );
			trueVisibleFourMomentumChargeds.push_back( getTrueVisibleFourMomentum( mcSLDLeptons.at( i_SLD ) , true , false ) );
			trueVisibleFourMomentumNeutrals.push_back( getTrueVisibleFourMomentum( mcSLDLeptons.at( i_SLD ) , false , true ) );
			trueFourMomentumNeutrinos.push_back( getTrueFourMomentumNeutrino( mcSLDLeptons.at( i_SLD ) ) );
			TLorentzVector trueVisibleFourMomentum = getTrueVisibleFourMomentum( mcSLDLeptons.at( i_SLD ) , true , false ) + getTrueVisibleFourMomentum( mcSLDLeptons.at( i_SLD ) , false , true );
			streamlog_out(DEBUG3) << "		              ( PDG	, Mass		, Px		, Py		, Pz		, E		)" << std::endl;
			streamlog_out(DEBUG3) << "		Lepton 4-Mom"  <<  std::endl;// 	( 	  Mass		, Px		, Py		, Pz		, E	)" << std::endl;
			streamlog_out(DEBUG3) << "				" << ( mcSLDLeptons.at( i_SLD ) )->getPDG() << "	  " << trueFourMomentumLepton.Mag() << "	, " << trueFourMomentumLepton.Px() << "	, " << trueFourMomentumLepton.Py() << "	, " << trueFourMomentumLepton.Pz() << "	, " << trueFourMomentumLepton.E() << "	)" << std::endl;
			streamlog_out(DEBUG3) << "		Visible 4-Mom"  <<  std::endl;// 	( 	  Mass		, Px		, Py		, Pz		, E	)" << std::endl;
			streamlog_out(DEBUG3) << "					  " << trueVisibleFourMomentum.Mag() << "	, " << trueVisibleFourMomentum.Px() << "	, " << trueVisibleFourMomentum.Py() << "	, " << trueVisibleFourMomentum.Pz() << "	, " << trueVisibleFourMomentum.E() << "	)" << std::endl;
			TLorentzVector trueInvisibleFourMomentum = getTrueFourMomentumNeutrino( mcSLDLeptons.at( i_SLD ) );
			streamlog_out(DEBUG3) << "		Inisible 4-Mom"  <<  std::endl;//	(	  Mass		, Px		, Py		, Pz		, E	)" << std::endl;
			streamlog_out(DEBUG3) << "					  " << trueInvisibleFourMomentum.Mag() << "	, " << trueInvisibleFourMomentum.Px() << "	, " << trueInvisibleFourMomentum.Py() << "	, " << trueInvisibleFourMomentum.Pz() << "	, " << trueInvisibleFourMomentum.E() << "	)" << std::endl;
			TLorentzVector trueFourMomentum = trueVisibleFourMomentum + trueInvisibleFourMomentum + trueFourMomentumLepton;
			streamlog_out(DEBUG3) << "		Total 4-Mom"  <<  std::endl;// 	( 	  Mass		, Px		, Py		, Pz		, E	)" << std::endl;
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

			TVector3 flightDirection = trueFlightDirections.at( i_SLD );
			TLorentzVector fourMomentumLepton = trueFourMomentumLeptons.at( i_SLD );
			TLorentzVector visibleFourMomentumCharged = trueVisibleFourMomentumChargeds.at( i_SLD );
			TLorentzVector visibleFourMomentumNeutral = trueVisibleFourMomentumNeutrals.at( i_SLD );
			double parentHadronMass = trueParentHadronMasses.at( i_SLD );
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "	Calculate Neutrino 4-Momentum ( + solution )" << std::endl;
			FourMomentumNeutrinosPos.push_back( getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , +1 ) );
			streamlog_out(DEBUG4) << "	True Neutrino 4-Momentum:		( " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Px() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Py() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Pz() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).E() << " )" << std::endl;
			streamlog_out(DEBUG4) << "	----------------------------------------------------------------------------------------------------" << std::endl;
			streamlog_out(DEBUG4) << "" << std::endl;
			streamlog_out(DEBUG4) << "	Calculate Neutrino 4-Momentum ( - solution )" << std::endl;
			FourMomentumNeutrinosNeg.push_back( getNeutrinoFourMomentum( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , -1 ) );
			streamlog_out(DEBUG4) << "	True Neutrino 4-Momentum:		( " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Px() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Py() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).Pz() << "	, " << ( trueFourMomentumNeutrinos.at( i_SLD ) ).E() << " )" << std::endl;
			streamlog_out(DEBUG4) << "	----------------------------------------------------------------------------------------------------" << std::endl;
			NeutrinoCovMat = getNeutrinoCovMat( flightDirection , fourMomentumLepton , visibleFourMomentumCharged , visibleFourMomentumNeutral , parentHadronMass , flightDirectionCovMat , LeptonCovMat , visibleChargedCovMat , visibleNeutralCovMat , sigmaParentHadronMass );
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

	TLorentzVector visibleFourMomentum = fourMomentumLepton + visibleFourMomentumCharged + visibleFourMomentumNeutral;
	double Mvis			= visibleFourMomentum.M();
	double Evis			= visibleFourMomentum.E();
	double EvisPrime		= ( pow( parentHadronMass , 2 ) + pow( Mvis , 2 ) ) / ( 2 * parentHadronMass );
	TVector3 visible_p		= visibleFourMomentum.Vect();
	TVector3 visible_p_par		= visible_p.Dot( flightDirection ) * flightDirection;
	TVector3 visible_p_nor		= visible_p - visible_p_par;
	TVector3 visible_p_par_prime	= sqrt( pow( ( pow( parentHadronMass , 2 ) - pow( Mvis , 2 ) ) / ( 2 * parentHadronMass ) , 2 ) - visible_p_nor.Mag2() ) * flightDirection;
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

	double dEnu_dEvis		= parentHadronMass * EvisPrime / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() ) - 1.0;
	double dEnu_dEvisPrime		= parentHadronMass * Evis / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );
	double dEnu_dPvisPar		= -parentHadronMass * visible_p_par_prime.Mag() / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );
	double dEnu_dPvisParPrime	= -parentHadronMass * visible_p_par.Mag() / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );
	double dEnu_dPvisNor		= -2 * visible_p_nor.Mag() * parentHadronMass * ( Evis * EvisPrime - visible_p_par.Dot( visible_p_par_prime ) ) / pow( pow( Mvis , 2 ) + visible_p_nor.Mag2() , 2 );
	double dEnu_dMvis		= -2 * Mvis * parentHadronMass * ( Evis * EvisPrime - visible_p_par.Dot( visible_p_par_prime ) ) / pow( pow( Mvis , 2 ) + visible_p_nor.Mag2() , 2 );
	double dEnu_dMparentHadron	= ( Evis * EvisPrime - visible_p_par.Dot( visible_p_par_prime ) ) / ( pow( Mvis , 2 ) + visible_p_nor.Mag2() );

	double sigmaEnu			= std::sqrt( 	pow( dEnu_dEvis , 2 ) * pow( sigmaEvis , 2 ) + pow( dEnu_dEvisPrime , 2 ) * pow( sigmaEvis_prime , 2 ) +
 							pow( dEnu_dPvisPar , 2 ) * pow( sigmaPvisPar , 2 ) + pow( dEnu_dPvisParPrime , 2 ) * pow( sigmaPvisPar_prime , 2 ) +
							pow( dEnu_dMvis , 2 ) * pow( sigmaMvis , 2 ) + pow( dEnu_dMparentHadron , 2 ) * pow( sigmaParentHadronMass , 2 ) +
							pow( dEnu_dPvisNor , 2 ) * pow( sigmaPvisNor , 2 ) 	);

	NeutrinoCovMat[ 9 ] = sigmaEnu;
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
	}

	streamlog_out(MESSAGE) << " " << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	processed events: 	" << m_nEvtSum << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << " " << std::endl;

}
