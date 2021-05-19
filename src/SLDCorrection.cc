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
		RecoJetCollection = pLCEvent->getCollection( m_inputJetCollection );
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
					streamlog_out(DEBUG3) << "		Number of tracks:	" << PFOtrkvec.size() << std::endl;
					streamlog_out(DEBUG3) << "		RadiusOfInnermostHit:	" << inputTrk->getRadiusOfInnermostHit() << std::endl;
					streamlog_out(DEBUG3) << "		Momentum (px,py,pz): 	( " << testPFO->getMomentum()[ 0 ] << "	, " << testPFO->getMomentum()[ 1 ] << "	 , " << testPFO->getMomentum()[ 2 ] << "	)" << std::endl;
//					streamlog_out(DEBUG3) << "		Vertex (x,y,z): 	( " << ( testPFO->getStartVertex() )->getPosition()[ 0 ] << "	, " << ( testPFO->getStartVertex() )->getPosition()[ 1 ] << "	 , " << ( testPFO->getStartVertex() )->getPosition()[ 2 ] << "	)" << std::endl;
				}
			}
		}


	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
        }

}

std::vector<MCParticle*> SLDCorrection::getMCLeptonVec( EVENT::LCEvent *pLCEvent , int ParentHadronPDG )
{
	LCCollection *MCParticleCollection{};
	std::vector<MCParticle*> mcLeptonVecSLD;
	std::string hadron = "";
	if ( ParentHadronPDG == 4 ) hadron = "C-Hadron";
	if ( ParentHadronPDG == 5 ) hadron = "B-Hadron";
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
					if ( floor( abs( testMotherHadron->getPDG() ) / 100 ) == ParentHadronPDG || ( floor( abs( testMotherHadron->getPDG() ) / 1000 ) == ParentHadronPDG ) )
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
