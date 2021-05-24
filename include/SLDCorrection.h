#ifndef SLDCorrection_h
#define SLDCorrection_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include "UTIL/LCRelationNavigator.h"
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
		ReconstructedParticle* getLinkedRecoLepton( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton );
		virtual void plotHistograms( TLorentzVector trueFourMomentumNeutrino , TLorentzVector FourMomentumNuClose , std::vector<float> NeutrinoCovMat );
		virtual void InitializeHistogram( TH1F *histogram , int scale , int color , int lineWidth , int markerSize , int markerStyle );
		virtual void doProperGaussianFit( TH1F *histogram , float fitMin , float fitMax , float fitRange );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();

	private:

		typedef std::vector<int>		IntVector;
		typedef std::vector<double>		DoubleVector;
		typedef std::vector<float>		FloatVector;

		std::string				m_mcParticleCollection{};
		std::string				m_inputPfoCollection{};
		std::string				m_inputJetCollection{};
		std::string				m_TrackMCTruthLinkCollection{};
		std::string				m_MCTruthTrackLinkCollection{};
		std::string				m_ClusterMCTruthLinkCollection{};
		std::string				m_MCTruthClusterLinkCollection{};
		std::string				m_SLDNuCollection{};
		std::string				m_rootFile{};

		bool					m_cheatSLDLeptons = true;
		bool					m_cheatFlightDirection = true;
		bool					m_cheatLepton4momentum = true;
		bool					m_fillRootTree = true;

		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		int					n_NuPxResidual;
		int					n_NuPyResidual;
		int					n_NuPzResidual;
		int					n_NuEResidual;
		int					n_NuPxNormalizedResidual;
		int					n_NuPyNormalizedResidual;
		int					n_NuPzNormalizedResidual;
		int					n_NuENormalizedResidual;
		DoubleVector				m_NuPxResidual{};
		DoubleVector				m_NuPyResidual{};
		DoubleVector				m_NuPzResidual{};
		DoubleVector				m_NuEResidual{};
		DoubleVector				m_NuPxNormalizedResidual{};
		DoubleVector				m_NuPyNormalizedResidual{};
		DoubleVector				m_NuPzNormalizedResidual{};
		DoubleVector				m_NuENormalizedResidual{};
		TH1F					*h_NuPxResidual{};
		TH1F					*h_NuPyResidual{};
		TH1F					*h_NuPzResidual{};
		TH1F					*h_NuEResidual{};
		TH1F					*h_NuPxNormalizedResidual{};
		TH1F					*h_NuPyNormalizedResidual{};
		TH1F					*h_NuPzNormalizedResidual{};
		TH1F					*h_NuENormalizedResidual{};
		TH2F					*h_recoNuPx_mcNuPx{};
		TH2F					*h_recoNuPy_mcNuPy{};
		TH2F					*h_recoNuPz_mcNuPz{};
		TH2F					*h_recoNuE_mcNuE{};
		TH1I					*h_recoPFOLinkedToElectron_Type{};
		TH1I					*h_recoPFOLinkedToMuon_Type{};
		TFile					*m_pTFile{};
		TTree					*m_pTTree{};

};

#endif
