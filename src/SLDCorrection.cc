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
	++m_nEvtSum;
	this->Clear();
	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

        try
        {

	}
	catch(DataNotAvailableException &e)
        {
        	streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
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
