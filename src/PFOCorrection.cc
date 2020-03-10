#include "PFOCorrection.h"
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
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;

PFOCorrection aPFOCorrection;

PFOCorrection::PFOCorrection() :

	Processor("PFOCorrection"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	m_MCTruthRecoLinkWeight(0.f),
	m_RecoMCTruthLinkWeight(0.f),
	m_MCTruthMarlinTrkTracksLinkWeight(0.f),
	m_MarlinTrkTracksMCTruthLinkWeight(0.f),
	m_nPfosTotal(0),
	m_nPfosTracks(0),
	m_nPfosPhotons(0),
	m_nPfosNeutralHadrons(0),
	m_pfoMomentumTotal(0.f),
	m_pfoEnergyTotal(0.f),
	m_pfoEnergyTracks(0.f),
	m_pfoEnergyPhotons(0.f),
	m_pfoEnergyNeutralHadrons(0.f),
	m_pfoPxTotal(0.f),
	m_pfoPyTotal(0.f),
	m_pfoPzTotal(0.f),
	m_pfoMomentumTotalOld(0.f),
	m_pfoEnergyTotalOld(0.f),
	m_pfoEnergyTracksOld(0.f),
	m_pfoPxTotalOld(0.f),
	m_pfoPyTotalOld(0.f),
	m_pfoPzTotalOld(0.f),
	m_Bfield(0.f),
	c(0.),
	mm2m(0.),
	eV2GeV(0.),
	eB(0.),
	proton_mass(0.),
	proton_mass_sq(0.),
	kaon_mass(0.),
	kaon_mass_sq(0.),
	pion_mass(0.),
	pion_mass_sq(0.)
{
	_description = "PFOCorrection creates new PandoraPFOs collection using tracks refitted with true mass of protons and kaons";

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

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollection" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracks,
					std::string("MarlinTrkTracks")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollectionKaon" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracksKAON,
					std::string("MarlinTrkTracksKaon")
				);

	registerInputCollection(	LCIO::TRACK,
					"MarlinTrkTracksCollectionProton" ,
					"Name of the MarlinTrkTracks collection"  ,
					m_MarlinTrkTracksPROTON,
					std::string("MarlinTrkTracksProton")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthRecoLink" ,
					"Name of the MCParticle-ReconstructedParticle Relations collection"  ,
					m_MCParticleRecoRelCol ,
					std::string("MCTruthRecoLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"RecoMCTruthLink" ,
					"Name of the ReconstructedParticle-MCParticle Relations collection"  ,
					m_RecoMCParticleRelCol ,
					std::string("RecoMCTruthLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthMarlinTrkTracksLink" ,
					"Name of the MCParticle-Track Relations collection for input tracks"  ,
					m_MCParticleTrackRelCol ,
					std::string("MCTruthMarlinTrkTracksLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MarlinTrkTracksMCTruthLink" ,
					"Name of the Track-MCParticle Relations collection for input tracks"  ,
					m_TrackMCParticleRelCol ,
					std::string("MarlinTrkTracksMCTruthLink")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"CorrectedPfoCollection",
					"Name of output pfo collection",
					m_outputPfoCollection,
					std::string("CorrectedPfoCollection")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"ModPFORecoLink" ,
					"Name of the PfoCollectionModified-ReconstructedParticle Relations collection"  ,
					m_ModPFORecoLink ,
					std::string("ModPFORecoLink")
				);

	registerOutputCollection(	LCIO::LCRELATION,
					"RecoModPFOLink" ,
					"Name of the ReconstructedParticle-PfoCollectionModified Relations collection"  ,
					m_RecoModPFOLink ,
					std::string("RecoModPFOLink")
				);

	registerProcessorParameter(	"MCTruthRecoLinkWeight" ,
					"Minimum acceptable weight for MCParticle-ReconstructedParticle Link"  ,
					m_MCTruthRecoLinkWeight ,
					float(0.9f)
				);

	registerProcessorParameter(	"RecoMCTruthLinkWeight" ,
					"Minimum acceptable weight for ReconstructedParticle-MCParticle Link"  ,
					m_RecoMCTruthLinkWeight ,
					float(0.9f)
				);

	registerProcessorParameter(	"MCTruthMarlinTrkTracksLinkWeight" ,
					"Minimum acceptable weight for Track-MCParticle Link"  ,
					m_MCTruthMarlinTrkTracksLinkWeight ,
					float(0.9f)
				);

	registerProcessorParameter(	"MarlinTrkTracksMCTruthLinkWeight" ,
					"Minimum acceptable weight for MCParticle-Track Link"  ,
					m_MarlinTrkTracksMCTruthLinkWeight ,
					float(0.9f)
				);

	registerProcessorParameter(	"updatePFO",
					"Update PFO.",
					m_updatePFO,
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
					std::string("PFOCorrection.root"));

}

void PFOCorrection::init()
{
	streamlog_out(DEBUG) << "   init called  " << std::endl;
	m_Bfield = MarlinUtil::getBzAtOrigin();
	printParameters();
	streamlog_out(DEBUG) << " BField =  "<< m_Bfield << " Tesla" << std::endl ;
	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	this->Clear();
	c = 2.99792458e8;
	mm2m = 1e-3;
	eV2GeV = 1e-9;
	eB = m_Bfield * c * mm2m * eV2GeV;
	proton_mass = 0.938272088;
	proton_mass_sq = proton_mass * proton_mass;
	kaon_mass = 0.493677;
	kaon_mass_sq = kaon_mass * kaon_mass;
	pion_mass = 0.13957018;
	pion_mass_sq = pion_mass * pion_mass;

	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");

	m_pTTree = new TTree("PfoCorrectionTree", "PfoCorrectionTree");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("maxWeight", &m_maxWeight);
	m_pTTree->Branch("pfoPDGCode", &m_pfoPDGCode);
	m_pTTree->Branch("pfoCosTheta", &m_pfoCosTheta);
	m_pTTree->Branch("pfoPx", &m_pfoPx);
        m_pTTree->Branch("pfoPy", &m_pfoPy);
        m_pTTree->Branch("pfoPz", &m_pfoPz);
	m_pTTree->Branch("pfoPt", &m_pfoPt);
	m_pTTree->Branch("pfoMomentum", &m_pfoMomentum);
	m_pTTree->Branch("pfoEnergy", &m_pfoEnergy);
	m_pTTree->Branch("pfoPxOld", &m_pfoPxOld);
        m_pTTree->Branch("pfoPyOld", &m_pfoPyOld);
        m_pTTree->Branch("pfoPzOld", &m_pfoPzOld);
	m_pTTree->Branch("pfoPtOld", &m_pfoPtOld);
	m_pTTree->Branch("pfoMomentumOld", &m_pfoMomentumOld);
	m_pTTree->Branch("pfoEnergyOld", &m_pfoEnergyOld);
	m_pTTree->Branch("pfoPxTotal", &m_pfoPxTotal, "pfoPxTotal/F");
	m_pTTree->Branch("pfoPyTotal", &m_pfoPyTotal, "pfoPyTotal/F");
	m_pTTree->Branch("pfoPzTotal", &m_pfoPzTotal, "pfoPzTotal/F");
	m_pTTree->Branch("pfoMomentumTotal", &m_pfoMomentumTotal, "pfoMomentumTotal/F");
	m_pTTree->Branch("pfoPxTotalOld", &m_pfoPxTotalOld, "pfoPxTotalOld/F");
	m_pTTree->Branch("pfoPyTotalOld", &m_pfoPyTotalOld, "pfoPyTotalOld/F");
	m_pTTree->Branch("pfoPzTotalOld", &m_pfoPzTotalOld, "pfoPzTotalOld/F");
	m_pTTree->Branch("pfoMomentumTotalOld", &m_pfoMomentumTotalOld, "pfoMomentumTotalOld/F");
	m_pTTree->Branch("nPfosTotal", &m_nPfosTotal, "nPfosTotal/I");
	m_pTTree->Branch("nPfosTracks", &m_nPfosTracks, "nPfosTracks/I");
        m_pTTree->Branch("nPfosPhotons", &m_nPfosPhotons, "nPfosPhotons/I");
	m_pTTree->Branch("nPfosNeutralHadrons", &m_nPfosNeutralHadrons, "nPfosNeutralHadrons/I");
	m_pTTree->Branch("nTracks", &m_nTracks);
        m_pTTree->Branch("pfoEnergyTotal", &m_pfoEnergyTotal, "pfoEnergyTotal/F");
	m_pTTree->Branch("pfoEnergyTracks", &m_pfoEnergyTracks, "pfoEnergyTracks/F");
        m_pTTree->Branch("pfoEnergyPhotons", &m_pfoEnergyPhotons, "pfoEnergyPhotons/F");
	m_pTTree->Branch("pfoEnergyNeutralHadrons", &m_pfoEnergyNeutralHadrons, "pfoEnergyNeutralHadrons/F");
	m_pTTree->Branch("pfoEnergyTotalOld", &m_pfoEnergyTotalOld, "pfoEnergyTotalOld/F");
	m_pTTree->Branch("pfoEnergyTracksOld", &m_pfoEnergyTracksOld, "pfoEnergyTracksOld/F");

	this->Clear();
}

void PFOCorrection::processRunHeader()
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}

void PFOCorrection::processEvent( EVENT::LCEvent *pLCEvent )
{

	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	++m_nEvtSum;

	LCCollection *inputPfoCollection{};
	LCCollection *MarlinTrkTracks{};
	LCCollection *MarlinTrkTracksKAON{};
	LCCollection *MarlinTrkTracksPROTON{};
	int n_PFO = -1;
	int n_TRK = -1;
	int n_TRKp = -1;
	int n_TRKk = -1;
	this->Clear();

        try
        {
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		MarlinTrkTracks = pLCEvent->getCollection(m_MarlinTrkTracks);
		MarlinTrkTracksKAON = pLCEvent->getCollection(m_MarlinTrkTracksKAON);
		MarlinTrkTracksPROTON = pLCEvent->getCollection(m_MarlinTrkTracksPROTON);
        }
        catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input collection not found in event " << m_nEvt << std::endl;
        }
	n_PFO = inputPfoCollection->getNumberOfElements();
	n_TRK = MarlinTrkTracks->getNumberOfElements();
	n_TRKk = MarlinTrkTracksKAON->getNumberOfElements();
	n_TRKp = MarlinTrkTracksPROTON->getNumberOfElements();
	if ( n_PFO == -1 ) streamlog_out(DEBUG) << "Input PFO collection (" << m_inputPfoCollection << ") has no element (PFO) " << std::endl;
	if ( n_TRK == -1 ) streamlog_out(DEBUG) << "Input TRACK collection (" << m_MarlinTrkTracks << ") has no element (Track) " << std::endl;
	if ( n_TRK != n_TRKk ) streamlog_out(DEBUG) << "Input TRACK collection (" << m_MarlinTrkTracksKAON << ") has in-equal number of elements with Main Track Collection (" << m_MarlinTrkTracks << ")!" << std::endl;
	if ( n_TRK != n_TRKp ) streamlog_out(DEBUG) << "Input TRACK collection (" << MarlinTrkTracksPROTON << ") has in-equal number of elements with Main Track Collection (" << m_MarlinTrkTracks << ")!" << std::endl;

	streamlog_out(DEBUG) << "Total Number of PFOs: " << n_PFO << std::endl;
	streamlog_out(DEBUG) << "Total Number of Tracks: " << n_TRK << std::endl;
	streamlog_out(DEBUG) << "Total Number of KaonTracks: " << n_TRKk << std::endl;
	streamlog_out(DEBUG) << "Total Number of ProtonTracks: " << n_TRKp << std::endl;

	LCCollectionVec *m_col_outputPfo = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

	for (int i_pfo = 0; i_pfo < n_PFO ; ++i_pfo)
	{
		ReconstructedParticleImpl* inputPFO = dynamic_cast<ReconstructedParticleImpl*>(inputPfoCollection->getElementAt(i_pfo));
		ReconstructedParticleImpl* corPFO = new ReconstructedParticleImpl;
		const EVENT::TrackVec& PFOtrkvec = inputPFO->getTracks();
		int nTRKsofPFO = PFOtrkvec.size();
		Track *refittedTrack = NULL;
		double maxweightTRKtoMCP = 0.;
		int pfo_PID = 0;
		double pfoMass = 0;
		double pfoMass_sq = pfoMass * pfoMass;
		double newd0 = 0.;
		double newz0 = 0.;
		double newPhi = 0.;
		double newOmega = 0.;
		double newTanLambda = 0.;
		std::vector<float> newCovMatrix{};
		int signTanLambda = 0;
		int signBz = 0;
		int signOmega = 0;
		double newCharge = 0.;
		double newcosLambda = 0.;
		double newsinLambda = 0.;
		double newcosPhi = 0.;
		double newsinPhi = 0.;
		double pT = 0.;
		double p = 0.;
		double px = 0.;
		double py = 0.;
		double pz = 0.;
		double energy = 0.;
		double newMomentum[3]{0., 0., 0.};
		m_updatePFO = true;
		if (nTRKsofPFO == 0)
		{
			m_updatePFO = false;
		}
		else if (nTRKsofPFO == 1)
		{
			m_updatePFO = true;
			Track *inputTrk = (Track*)PFOtrkvec.at(0);
			std::vector<float> PFOMCPlink = this->CheckPfoPID(pLCEvent, inputTrk);
			pfo_PID = PFOMCPlink[0];
			maxweightTRKtoMCP = PFOMCPlink[1];
			int track_index = this->FindTrackIndex(pLCEvent, inputTrk);
			if ( fabs(pfo_PID) == 2212 && maxweightTRKtoMCP >= m_RecoMCTruthLinkWeight )
			{
				refittedTrack = dynamic_cast<EVENT::Track*>(MarlinTrkTracksPROTON->getElementAt(track_index));
				pfoMass = proton_mass;
				pfoMass_sq = proton_mass_sq;
			}
			else if ( fabs(pfo_PID) == 321 && maxweightTRKtoMCP >= m_RecoMCTruthLinkWeight )
			{
				refittedTrack = dynamic_cast<EVENT::Track*>(MarlinTrkTracksKAON->getElementAt(track_index));
				pfoMass = kaon_mass;
				pfoMass_sq = kaon_mass_sq;
			}
			else
			{
				refittedTrack = dynamic_cast<EVENT::Track*>(MarlinTrkTracks->getElementAt(track_index));
				pfoMass = pion_mass;
				pfoMass_sq = pion_mass_sq;
			}

			newd0 = refittedTrack->getD0();
			newz0 = refittedTrack->getZ0();
			newPhi = refittedTrack->getPhi();
			newOmega = refittedTrack->getOmega();
			newTanLambda = refittedTrack->getTanLambda();
			newCovMatrix = refittedTrack->getCovMatrix();
			signTanLambda = (newTanLambda > 0 ? 1 : -1);
			signBz = (m_Bfield > 0 ? 1 : -1);
			signOmega = (newOmega > 0 ? 1 : -1);
			newCharge = signBz * signOmega;
			newcosLambda = 1 / sqrt( 1 + pow(newTanLambda,2)) * signTanLambda;
			newsinLambda = newcosLambda * newTanLambda;
			newcosPhi = cos(newPhi);
			newsinPhi = sin(newPhi);
			pT = eB / fabs(newOmega);
			p = pT / newcosLambda;
			px = p * newcosLambda * newcosPhi;
			py = p * newcosLambda * newsinPhi;
			pz = p * newsinLambda;
			energy = sqrt( pfoMass_sq + px * px + py * py + pz * pz);
			newMomentum[0] = px;
			newMomentum[1] = py;
			newMomentum[2] = pz;

			streamlog_out(DEBUG) << "PFO has ONE track with PDG code: " << pfo_PID << " at Index: " << track_index << std::endl;

		}
		else
		{
			streamlog_out(DEBUG) << "PFO has mare than ONE track ; neglected at the moment " << std::endl;
			m_updatePFO = false;
		}


		if (m_updatePFO)
		{
			corPFO->setType(pfo_PID);
			corPFO->setMomentum(newMomentum);
			corPFO->setEnergy(energy);
			corPFO->setCovMatrix(inputPFO->getCovMatrix());
			corPFO->setMass(pfoMass);
			corPFO->setCharge(newCharge);
			corPFO->setReferencePoint(inputPFO->getReferencePoint());
			for (unsigned int j=0; j<inputPFO->getParticleIDs().size(); ++j)
			{
				ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>(inputPFO->getParticleIDs()[j]);
			        ParticleIDImpl* outPID = new ParticleIDImpl;
			        outPID->setType(inPID->getType());
			        outPID->setPDG(pfo_PID);
			        outPID->setLikelihood(inPID->getLikelihood());
			        outPID->setAlgorithmType(inPID->getAlgorithmType()) ;
			        for (unsigned int k=0; k<inPID->getParameters().size()  ; ++k) outPID->addParameter(inPID->getParameters()[k]) ;
			        corPFO->addParticleID(outPID);
			}
			corPFO->setParticleIDUsed(inputPFO->getParticleIDUsed());
			corPFO->setGoodnessOfPID(inputPFO->getGoodnessOfPID());
			for (unsigned int j=0; j<inputPFO->getParticles().size(); ++j)
			{
				corPFO->addParticle(inputPFO->getParticles()[j]);
			}
			for (unsigned int j=0; j<inputPFO->getClusters().size(); ++j)
			{
				corPFO->addCluster(inputPFO->getClusters()[j]);
			}
			corPFO->addTrack(refittedTrack);
			corPFO->setStartVertex(inputPFO->getStartVertex());
		}
		else
		{
			corPFO->setType(inputPFO->getType());
			corPFO->setMomentum(inputPFO->getMomentum());
			corPFO->setEnergy(inputPFO->getEnergy());
			corPFO->setCovMatrix(inputPFO->getCovMatrix());
			corPFO->setMass(inputPFO->getMass());
			corPFO->setCharge(inputPFO->getCharge());
			corPFO->setReferencePoint(inputPFO->getReferencePoint());
			for (unsigned int j=0; j<inputPFO->getParticleIDs().size(); ++j)
			{
				ParticleIDImpl* inPID = dynamic_cast<ParticleIDImpl*>(inputPFO->getParticleIDs()[j]);
			        ParticleIDImpl* outPID = new ParticleIDImpl;
			        outPID->setType(inPID->getType());
			        outPID->setPDG(inPID->getPDG());
			        outPID->setLikelihood(inPID->getLikelihood());
			        outPID->setAlgorithmType(inPID->getAlgorithmType()) ;
			        for (unsigned int k=0; k<inPID->getParameters().size()  ; ++k) outPID->addParameter(inPID->getParameters()[k]) ;
			        corPFO->addParticleID(outPID);
			}
			corPFO->setParticleIDUsed(inputPFO->getParticleIDUsed());
			corPFO->setGoodnessOfPID(inputPFO->getGoodnessOfPID());
			for (unsigned int j=0; j<inputPFO->getParticles().size(); ++j)
			{
				corPFO->addParticle(inputPFO->getParticles()[j]);
			}
			for (unsigned int j=0; j<inputPFO->getClusters().size(); ++j)
			{
				corPFO->addCluster(inputPFO->getClusters()[j]);
			}
			for (unsigned int j=0; j<inputPFO->getTracks().size(); ++j)
			{
				corPFO->addTrack(inputPFO->getTracks()[j]);
			}
			corPFO->setStartVertex(inputPFO->getStartVertex());
		}

		if (inputPFO->getTracks().size() == 0 && 22 == inputPFO->getType())
		{
			++m_nPfosPhotons;
			m_pfoEnergyPhotons += corPFO->getEnergy();
		}
		else if (inputPFO->getTracks().size() == 0 && 22 != inputPFO->getType())
		{
			m_nPfosNeutralHadrons++;
			m_pfoEnergyNeutralHadrons += corPFO->getEnergy();
		}
		else
		{
			++m_nPfosTracks;
	                m_pfoEnergyTracks += corPFO->getEnergy();
			m_pfoEnergyTracksOld += inputPFO->getEnergy();
		}
		m_nTracks.push_back(nTRKsofPFO);
		++m_nPfosTotal;
		TLorentzVector pfoFourMomentum(corPFO->getMomentum(),corPFO->getEnergy());
		TVector3 pfoMomentum(corPFO->getMomentum());
		TVector3 pfoMomentumOld(inputPFO->getMomentum());
		m_maxWeight.push_back(maxweightTRKtoMCP);
		m_pfoPDGCode.push_back(pfo_PID);
		m_pfoCosTheta.push_back(pfoFourMomentum.CosTheta());
		m_pfoPx.push_back(pfoMomentum[0]);
		m_pfoPy.push_back(pfoMomentum[1]);
		m_pfoPz.push_back(pfoMomentum[2]);
		m_pfoPt.push_back(pfoMomentum.Perp());
		m_pfoMomentum.push_back(pfoMomentum.Mag());
		m_pfoEnergy.push_back(corPFO->getEnergy());
		m_pfoPxTotal += pfoMomentum[0];
		m_pfoPyTotal += pfoMomentum[1];
		m_pfoPzTotal += pfoMomentum[2];
		m_pfoPxOld.push_back(pfoMomentumOld[0]);
		m_pfoPyOld.push_back(pfoMomentumOld[1]);
		m_pfoPzOld.push_back(pfoMomentumOld[2]);
		m_pfoPtOld.push_back(pfoMomentumOld.Perp());
		m_pfoMomentumTotal += pfoMomentum.Mag();
		m_pfoEnergyTotal += corPFO->getEnergy();
		m_pfoMomentumOld.push_back(pfoMomentumOld.Mag());
		m_pfoEnergyOld.push_back(inputPFO->getEnergy());
		m_pfoPxTotalOld += pfoMomentumOld[0];
		m_pfoPyTotalOld += pfoMomentumOld[1];
		m_pfoPzTotalOld += pfoMomentumOld[2];
		m_pfoMomentumTotalOld += pfoMomentumOld.Mag();
		m_pfoEnergyTotalOld += inputPFO->getEnergy();

		m_col_outputPfo->addElement( corPFO );

	}
	pLCEvent->addCollection( m_col_outputPfo , m_outputPfoCollection );

	if ((m_nEvtSum % 100) == 0)
		std::cout << " processed events: " << m_nEvtSum << std::endl;

	if (m_fillRootTree) m_pTTree->Fill();
	m_nEvt++;

}

std::vector<float> PFOCorrection::CheckPfoPID(EVENT::LCEvent *pLCEvent, EVENT::Track* inputTrk)
{
	std::vector<float> PFOMCPlink;
	try
	{
		LCRelationNavigator TrackMCParticleNav(pLCEvent->getCollection(m_TrackMCParticleRelCol));
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the LCRelationNavigator from Track to MCParticle" << std::endl;
		return {0.,0.};
	}
	LCRelationNavigator TrackMCParticleNav(pLCEvent->getCollection(m_TrackMCParticleRelCol));
	const EVENT::LCObjectVec& mcpvec = TrackMCParticleNav.getRelatedToObjects(inputTrk);
	const EVENT::FloatVec&  trkweightvec = TrackMCParticleNav.getRelatedToWeights(inputTrk);
	double maxweightTRKtoMCP = 0.;
	int iTRKtoMCPmax = -1;
	for ( unsigned int i_mcp = 0; i_mcp < mcpvec.size(); i_mcp++ )
	{
		double track_weight = trkweightvec.at(i_mcp);
		MCParticle *testMCP = (MCParticle *) mcpvec.at(i_mcp);
		if ( (fabs(testMCP->getPDG()) == 2212 || fabs(testMCP->getPDG()) == 321) && track_weight > maxweightTRKtoMCP)
		{
			maxweightTRKtoMCP = track_weight;
			iTRKtoMCPmax = i_mcp;
			streamlog_out(DEBUG) << "track at index: " << i_mcp << " has PDG: " << testMCP->getPDG() << " and weight to MCP = " << track_weight << std::endl;
		}
	}
	int pfo_PDG = 0;
	if (iTRKtoMCPmax != -1 )//&& maxweightTRKtoMCP >= m_RecoMCTruthLinkWeight )
	{
		MCParticle *linkedMCP = (MCParticle *) mcpvec.at(iTRKtoMCPmax);
		pfo_PDG = linkedMCP->getPDG();
	}
	else
	{
		pfo_PDG = 0;
	}
	PFOMCPlink.push_back(pfo_PDG);
	PFOMCPlink.push_back(maxweightTRKtoMCP);
	return PFOMCPlink;
}
int PFOCorrection::FindTrackIndex(EVENT::LCEvent *pLCEvent, EVENT::Track* inputTrk)
{
	LCCollection *MarlinTrkTracks{};
	try
	{
		MarlinTrkTracks = pLCEvent->getCollection(m_MarlinTrkTracks);
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(WARNING) << "Could not find the MarlinTrack Collection" << std::endl;
		return -1;
	}
	unsigned int nTRKs = MarlinTrkTracks->getNumberOfElements();
	int track_index = -1;
	for (unsigned int i_trk = 0; i_trk < nTRKs;  ++i_trk )
	{
		Track* pionTrack = dynamic_cast<EVENT::Track*>(MarlinTrkTracks->getElementAt(i_trk));
		if ( pionTrack == inputTrk)
		{
			track_index = i_trk;
		}
	}
	if (track_index == -1)
	{
		streamlog_out(DEBUG) << "Coudln't find track_index!!!  " << std::endl;
	}
	else
	{
		streamlog_out(DEBUG) << "Track index in " << m_MarlinTrkTracks << " collection is found; track_index = " << track_index << std::endl;
	}
	return track_index;
}
void PFOCorrection::Clear()
{
	m_pfoPx.clear();
	m_pfoPy.clear();
	m_pfoPz.clear();
	m_pfoPt.clear();
	m_pfoMomentum.clear();
	m_pfoEnergy.clear();
	m_pfoPxOld.clear();
	m_pfoPyOld.clear();
	m_pfoPzOld.clear();
	m_pfoPtOld.clear();
	m_pfoMomentumOld.clear();
	m_pfoEnergyOld.clear();
	m_maxWeight.clear();
	m_pfoPDGCode.clear();
	m_pfoCosTheta.clear();
	m_nTracks.clear();

	m_pfoEnergyTotal = 0.f;
	m_pfoMomentumTotal = 0.f;
	m_pfoEnergyTracks = 0.f;
	m_pfoEnergyPhotons = 0.f;
	m_pfoEnergyNeutralHadrons = 0.f;
	m_pfoEnergyTotalOld = 0.f;
	m_pfoMomentumTotalOld = 0.f;
	m_pfoEnergyTracksOld = 0.f;
	m_nPfosTotal = 0;
	m_nPfosTracks = 0;
	m_nPfosPhotons = 0;
	m_nPfosNeutralHadrons = 0;
	m_pfoPxTotal = 0.f;
	m_pfoPyTotal = 0.f;
	m_pfoPzTotal = 0.f;
	m_pfoPxTotalOld = 0.f;
	m_pfoPyTotalOld = 0.f;
	m_pfoPzTotalOld = 0.f;
}

void PFOCorrection::check(EVENT::LCEvent *pLCEvent)
{
	LCCollection *inputPfoCollection{};
	LCCollection *outputPfoCollection{};
	try
	{
		inputPfoCollection = pLCEvent->getCollection(m_inputPfoCollection);
		outputPfoCollection = pLCEvent->getCollection(m_outputPfoCollection);
	}
	catch(DataNotAvailableException &e)
        {
          streamlog_out(MESSAGE) << "Input/Output collection not found in event " << m_nEvt << std::endl;
        }
	int n_inputPFOs = inputPfoCollection->getNumberOfElements();
	int n_outputPFOs = outputPfoCollection->getNumberOfElements();
	streamlog_out(DEBUG) << " CHECK : processed events: " << m_nEvtSum << " (Number of inputPFOS: " << n_inputPFOs << " , Number of outputPFOs: " << n_outputPFOs <<")" << std::endl;
}
void PFOCorrection::end()
{

	m_pTFile->cd();
	m_pTTree->Write();
//
	m_pTFile->Close();
	delete m_pTFile;
	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}
/*
void PFOCorrection::Clear()
{
}
*/
