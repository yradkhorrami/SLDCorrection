#ifndef SLDCorrection_h
#define SLDCorrection_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include <IMPL/ReconstructedParticleImpl.h>
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixD.h"
class TFile;
class TDirectory;
class TH1F;
class TH1I;
class TH2I;
class TH2F;
class TTree;
class TF1;

using namespace lcio ;
using namespace marlin ;

class SLDCorrection : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new SLDCorrection;
		}
		SLDCorrection();
		virtual ~SLDCorrection() = default;
		SLDCorrection(const SLDCorrection&) = delete;
		SLDCorrection& operator=(const SLDCorrection&) = delete;
		virtual void init();
		virtual void Clear();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
//		virtual void getMCLeptonVec( EVENT::LCEvent *pLCEvent , EVENT::MCParticleVec &mcLeptonVecBSLD , EVENT::MCParticleVec &mcLeptonVecCSLD );
		std::vector<MCParticle*> getMCLeptonVec( EVENT::LCEvent *pLCEvent , int parentHadronPDG );
		TVector3 getTrueFlightDirection( MCParticle *SLDLepton );
		TLorentzVector getTrueVisibleFourMomentum( MCParticle *SLDLepton , bool getChargedTLV , bool getNeutralTLV );
		TLorentzVector getTrueFourMomentumNeutrino( MCParticle *SLDLepton );
		TLorentzVector getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector fourMomentumLepton , TLorentzVector visibleFourMomentumCharged , TLorentzVector visibleFourMomentumNeutral , double parentHadronMass , int solutionSign );
		std::vector< float > getNeutrinoCovMat( TVector3 flightDirection , TLorentzVector fourMomentumLepton , TLorentzVector visibleFourMomentumCharged , TLorentzVector visibleFourMomentumNeutral , double parentHadronMass , std::vector< float > flightDirectionCovMat , std::vector< float > leptonCovMat , std::vector< float > visibleChargedCovMat , std::vector< float > visibleNeutralCovMat , double sigmaParentHadronMass );
		double getSigmaVisibleMass( TLorentzVector visibleFourMomentum , std::vector< float > visibleCovMat );
		double getSigmaVisiblePpar_prime( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , std::vector< float > flightDirectionCovMat , std::vector< float > visibleCovMat , double sigmaParentHadronMass , double sigma_Mvis );
		double getSigmaPvisPar( TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > flightDirectionCovMat , std::vector< float > visibleCovMat );
		double getSigmaPvisNor( TVector3 visibleMomentum , double visibleMomentumPar , std::vector< float > visibleCovMat , double sigmaPvisPar );
		std::vector<float> getParentHadronCovMat( TVector3 flightDirection , double parentHadronEnergy , double parentHadronMass , std::vector<float> flightDirectionCovMat , double parentHadronSigmaE );
		double getTrueParentHadronMass( MCParticle *SLDLepton );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:

		typedef std::vector<int>		IntVector;
		typedef std::vector<double>		DoubleVector;
		typedef std::vector<float>		FloatVector;

		std::string				m_mcParticleCollection{};
		std::string				m_inputPfoCollection{};
		std::string				m_inputJetCollection{};
		std::string				m_SLDNuCollection{};
		std::string				m_rootFile{};

		bool					m_cheatSLDLeptons = true;
		bool					m_cheatFlightDirection = true;
		bool					m_fillRootTree = true;

		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		TFile					*m_pTFile{};
		TTree					*m_pTTree{};

};

#endif
