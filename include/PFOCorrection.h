#ifndef MyTrackParameters_h
#define MyTrackParameters_h 1
#include <marlin/Processor.h>
#include <marlin/Global.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
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

class TFile;
class TH1F;
class TTree;

using namespace lcio ;
using namespace marlin ;

class PFOCorrection : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new PFOCorrection;
		}
		PFOCorrection();
		virtual ~PFOCorrection() = default;
		PFOCorrection(const PFOCorrection&) = delete;
		PFOCorrection& operator=(const PFOCorrection&) = delete;
		virtual void init();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();
		std::vector<float> CheckPfoPID(EVENT::LCEvent *sLCEvent, EVENT::Track* inputTrk);
		int FindTrackIndex(EVENT::LCEvent *sLCEvent, EVENT::Track* inputTrk);
		void Clear();

	private:

		typedef std::vector<int>		IntVector;
		typedef std::vector<double>		DoubleVector;
		typedef std::vector<float>		FloatVector;

		std::string				m_mcParticleCollection{};
		std::string				m_inputPfoCollection{};
		std::string				m_MarlinTrkTracks{};
		std::string				m_MarlinTrkTracksKAON{};
		std::string				m_MarlinTrkTracksPROTON{};
		std::string				m_MCParticleRecoRelCol{};
		std::string				m_RecoMCParticleRelCol{};
		std::string				m_MCParticleTrackRelCol{};
		std::string				m_TrackMCParticleRelCol{};
		std::string				m_outputPfoCollection{};
		std::string				m_ModPFORecoLink{};
		std::string				m_RecoModPFOLink{};
		std::string				m_rootFile{};

		bool					m_updatePFO = true;
		bool					m_fillRootTree = true;


		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		float					m_MCTruthRecoLinkWeight;
		float					m_RecoMCTruthLinkWeight;
		float					m_MCTruthMarlinTrkTracksLinkWeight;
		float					m_MarlinTrkTracksMCTruthLinkWeight;
		int					m_nPfosTotal;
		int					m_nPfosTracks;
		int					m_nPfosPhotons;
		int					m_nPfosNeutralHadrons;
		float					m_pfoMomentumTotal;
		float					m_pfoEnergyTotal;
		float					m_pfoEnergyTracks;
		float					m_pfoEnergyPhotons;
		float					m_pfoEnergyNeutralHadrons;
		float					m_pfoPxTotal;
		float					m_pfoPyTotal;
		float					m_pfoPzTotal;
		float					m_pfoMomentumTotalOld;
		float					m_pfoEnergyTotalOld;
		float					m_pfoEnergyTracksOld;
		float					m_pfoPxTotalOld;
		float					m_pfoPyTotalOld;
		float					m_pfoPzTotalOld;
		float					m_Bfield;
		double					c;
		double					mm2m;
		double					eV2GeV;
		double					eB;
		double					proton_mass;
		double					proton_mass_sq;
		double					kaon_mass;
		double					kaon_mass_sq;
		double					pion_mass;
		double					pion_mass_sq;
		FloatVector				m_maxWeight{};
		FloatVector				m_pfoPDGCode{};
		FloatVector				m_pfoCosTheta{};
		FloatVector				m_pfoPx{};
		FloatVector				m_pfoPy{};
		FloatVector				m_pfoPz{};
		FloatVector				m_pfoPt{};
		FloatVector				m_pfoMomentum{};
		FloatVector				m_pfoEnergy{};
		FloatVector				m_pfoPxOld{};
		FloatVector				m_pfoPyOld{};
		FloatVector				m_pfoPzOld{};
		FloatVector				m_pfoPtOld{};
		FloatVector				m_pfoMomentumOld{};
		FloatVector				m_pfoEnergyOld{};
		IntVector				m_nTracks{};
		TFile					*m_pTFile{};
	        TTree					*m_pTTree{};

};

#endif
