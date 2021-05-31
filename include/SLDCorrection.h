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
		bool hasPrimarySLDecay( MCParticle *parentHadron );
		bool hasSecondaySLDecay( MCParticle *parentHadron );
		virtual void doSLDCorrection( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton );
		TVector3 getFlightDirection( MCParticle *SLDLepton );
		double getParentHadronMass( MCParticle *SLDLepton );
		TLorentzVector getLeptonFourMomentum( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton );
		TLorentzVector getVisibleFourMomentum( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton , MCParticle *parentHadron , bool getChargedTLV , bool getNeutralTLV );
		TLorentzVector getTrueNeutrinoFourMomentum( MCParticle *SLDLepton );

		TLorentzVector getParentHadron4mom( MCParticle *SLDLepton );
		TLorentzVector getDaughters4mom( MCParticle *SLDLepton );

		TLorentzVector visibleFourMomentumCharged( MCParticle *SLDLepton , bool getChargedTLV , bool getNeutralTLV );


//		virtual void getMCLeptonVec( EVENT::LCEvent *pLCEvent , EVENT::MCParticleVec &mcLeptonVecBSLD , EVENT::MCParticleVec &mcLeptonVecCSLD );
		std::vector<MCParticle*> getMCLeptonVec( EVENT::LCEvent *pLCEvent , int parentHadronPDG );
		std::vector<MCParticle*> getMCVisibles( MCParticle *mcParentHadron , bool getChargedTLV , bool getNeutralTLV );
		TLorentzVector getNeutrinoFourMomentum( TVector3 flightDirection , TLorentzVector fourMomentumLepton , TLorentzVector visibleFourMomentumCharged , TLorentzVector visibleFourMomentumNeutral , double parentHadronMass , int solutionSign );
		std::vector< float > getNeutrinoCovMat( TVector3 flightDirection , TLorentzVector fourMomentumLepton , TLorentzVector visibleFourMomentumCharged , TLorentzVector visibleFourMomentumNeutral , double parentHadronMass , std::vector< float > flightDirectionCovMat , std::vector< float > leptonCovMat , std::vector< float > visibleChargedCovMat , std::vector< float > visibleNeutralCovMat , double sigmaParentHadronMass );
		double getSigmaVisibleMass( TLorentzVector visibleFourMomentum , std::vector< float > visibleCovMat );
		double getSigmaVisiblePpar_prime( TVector3 flightDirection , TLorentzVector visibleFourMomentum , double parentHadronMass , std::vector< float > flightDirectionCovMat , std::vector< float > visibleCovMat , double sigmaParentHadronMass , double sigma_Mvis );
		double getSigmaPvisPar( TVector3 flightDirection , TLorentzVector visibleFourMomentum , std::vector< float > flightDirectionCovMat , std::vector< float > visibleCovMat );
		double getSigmaPvisNor( TVector3 visibleMomentum , double visibleMomentumPar , std::vector< float > visibleCovMat , double sigmaPvisPar );
		std::vector<float> getParentHadronCovMat( TVector3 flightDirection , double parentHadronEnergy , double parentHadronMass , std::vector<float> flightDirectionCovMat , double parentHadronSigmaE );
		ReconstructedParticle* getLinkedChargedPFO( EVENT::LCEvent *pLCEvent , MCParticle *SLDLepton );
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
		bool					m_cheatCharged4momentum = true;
		bool					m_cheatNeutral4momentum = true;
		bool					m_fillRootTree = true;

		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		int					m_nTauSLDecay;
		int					m_nNeutrino;
		int					m_nChargedPFOwoTrack;
		IntVector				m_nSLD_chargedMCPwoTrack{};
		IntVector				m_GenStatParentHadron{};
		int					n_NuPxResidual;
		int					n_NuPyResidual;
		int					n_NuPzResidual;
		int					n_NuEResidual;
		int					n_NuPxNormalizedResidual;
		int					n_NuPyNormalizedResidual;
		int					n_NuPzNormalizedResidual;
		int					n_NuENormalizedResidual;
		DoubleVector				m_BSLDecayX{};
		DoubleVector				m_BSLDecayY{};
		DoubleVector				m_BSLDecayZ{};
		DoubleVector				m_BSLDecayR{};
		DoubleVector				m_CSLDecayX{};
		DoubleVector				m_CSLDecayY{};
		DoubleVector				m_CSLDecayZ{};
		DoubleVector				m_CSLDecayR{};
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
		TH2F					*h_parentPx_daughtersPx{};
		TH2F					*h_parentPy_daughtersPy{};
		TH2F					*h_parentPz_daughtersPz{};
		TH2F					*h_parentE_daughtersE{};
		TH1I					*h_recoPFOLinkedToElectron_Type{};
		TH1I					*h_recoPFOLinkedToMuon_Type{};
		TH1I					*h_SLDecayOrder{};
		TFile					*m_pTFile{};
		TTree					*m_pTTree{};

};

#endif
