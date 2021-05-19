#include "SLDCorrection.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include "EVENT/LCCollection.h"
#include "IMPL/LCCollectionVec.h"
#include "EVENT/MCParticle.h"
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
		std::vector<const EVENT::MCParticle*> mcLeptonVecBSLD;
		std::vector<const EVENT::MCParticle*> mcLeptonVecCSLD;
		getMCLeptonVec( pLCEvent , mcLeptonVecBSLD , mcLeptonVecCSLD );
		RecoJetCollection = pLCEvent->getCollection( m_inputJetCollection );

	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
        }

}

void SLDCorrection::getMCLeptonVec( EVENT::LCEvent *pLCEvent , std::vector<const EVENT::MCParticle*> mcLeptonVecBSLD , std::vector<const EVENT::MCParticle*> mcLeptonVecCSLD )
{
	LCCollection *MCParticleCollection{};
	try
        {
		MCParticleCollection = pLCEvent->getCollection( m_mcParticleCollection );
		int nMCP = MCParticleCollection->getNumberOfElements();
		for ( int i_mcp = 0 ; i_mcp < nMCP ; ++i_mcp )
		{
			const EVENT::MCParticle *testLepton = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_mcp ) );
			int nLeptonsInVertexBSLD = 0;
			int nLeptonsInVertexCSLD = 0;
			bool foundBSLD = false;
			bool foundCSLD = false;
			if ( ( abs( testLepton->getPDG() ) == 11 || abs( testLepton->getPDG() ) == 13 ) && ( testLepton->getGeneratorStatus() == 1 ) )
			{
				for ( long unsigned int i_parent = 0 ; i_parent < ( testLepton->getParents() ).size() ; ++i_parent )
				{
					const EVENT::MCParticle *testMotherHadron = testLepton->getParents()[ i_parent ];
					if ( floor( abs( testMotherHadron->getPDG() ) / 100 ) == 5 || ( floor( abs( testMotherHadron->getPDG() ) / 1000 ) == 5 ) )
					{
						for ( long unsigned int i_daughter = 0 ; i_daughter < ( testMotherHadron->getDaughters() ).size() ; ++i_daughter )
						{
							if ( ( abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 11 || abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 13 ) && ( testMotherHadron->getDaughters()[ i_daughter ] )->getGeneratorStatus() == 1 )
							{
								++nLeptonsInVertexBSLD;
							}
							if ( abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == abs( testLepton->getPDG() ) + 1 )
							{
								mcLeptonVecBSLD.push_back( testLepton );
								foundBSLD = true;
							}
						}
					}
					if ( floor( abs( testMotherHadron->getPDG() ) / 100 ) == 4 || ( floor( abs( testMotherHadron->getPDG() ) / 1000 ) == 4 ) )
					{
						for ( long unsigned int i_daughter = 0 ; i_daughter < ( testMotherHadron->getDaughters() ).size() ; ++i_daughter )
						{
							if ( ( abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 11 || abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == 13 ) && ( testMotherHadron->getDaughters()[ i_daughter ] )->getGeneratorStatus() == 1 )
							{
								++nLeptonsInVertexCSLD;
							}
							if ( abs( ( testMotherHadron->getDaughters()[ i_daughter ] )->getPDG() ) == abs( testLepton->getPDG() ) + 1 )
							{
								mcLeptonVecCSLD.push_back( testLepton );
								foundCSLD = true;
							}
						}
					}
				}
				if ( foundBSLD ) streamlog_out(DEBUG0) << "	One Semi-Leptonic Decay of B-Hadron is found ; Number of leptons associated to SLD vertex : " << nLeptonsInVertexBSLD << std::endl;
				if ( foundCSLD ) streamlog_out(DEBUG0) << "	One Semi-Leptonic Decay of C-Hadron is found ; Number of leptons associated to SLD vertex : " << nLeptonsInVertexCSLD << std::endl;
			}
		}
	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input " << MCParticleCollection << " collection not found in event " << pLCEvent->getEventNumber() << std::endl;
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
	}

	streamlog_out(MESSAGE) << " " << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	processed events: 	" << m_nEvtSum << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << " " << std::endl;

}
