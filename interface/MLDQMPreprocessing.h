#ifndef MLDQMPreprocessing_h
#define MLDQMPreprocessing_h

#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"

#include "DQM/HcalCommon/interface/Utilities.h"
#include "HCALDQM/Preprocessing/interface/SubdetPMHist.h"

#include <string>
#include <map>
#include <vector>
#include <TTree.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>

enum HcalSubdetectorPM {
	HBM,
	HBP,
	HEM,
	HEP,
	HOM,
	HOP,
	HFM,
	HFP
};


// Thread histogram container
class MLDQMPreprocessingStreamData {
public:
	int _nevents;
	std::map<TString, SubdetPMHist*> _hists;

	inline void Init() {
		_nevents = 0;
		_hists["occupancy_20fC"] = new SubdetPMHist("occupancy_20fC", "occupancy_20fC");
		_hists["sumQ"]  = new SubdetPMHist("sumQ", "sumQ");
	}

	inline void Add(const MLDQMPreprocessingStreamData* data) {
		for (auto& it : _hists) {
			(it.second)->Add(data->_hists.at(it.first));
		}
		_nevents += data->_nevents;
	}

	inline void Finish() {
		// Normalize histograms to number of events
		if (_nevents > 0) {
			for (auto& it : _hists) {
				(it.second)->Scale(1. / _nevents);
			}
		}
	}

	inline void Reset() {
		for (auto& it : _hists) {
			(it.second)->Reset();
		}
		_nevents = 0;
	}

	inline void Write() {
		for (auto& it : _hists) {
			(it.second)->Write();
		}
	}
};

class MLDQMPreprocessing :  
	public edm::global::EDAnalyzer<edm::StreamCache<MLDQMPreprocessingStreamData>, edm::RunSummaryCache<MLDQMPreprocessingStreamData>> {
public:
	edm::InputTag		_tagHBHE;
	edm::EDGetTokenT<HBHEDigiCollection> _tokHBHE; // Old QIE8 HBHE digis

	edm::InputTag		_tagQIE11;
	edm::EDGetTokenT<QIE11DigiCollection> _tokQIE11; // New QIE11 HBHE digis

	edm::InputTag		_tagQIE10;
	edm::EDGetTokenT<QIE10DigiCollection> _tokQIE10; // New QIE10 HF digis

	edm::InputTag		_tagHO;
	edm::EDGetTokenT<HODigiCollection> _tokHO; // QIE8 HO digis

public:
	explicit MLDQMPreprocessing(const edm::ParameterSet& ps);// : pset(iConfig) {}
	~MLDQMPreprocessing();

 private:
	
	void beginJob();
	void beginRun(edm::Run const&, edm::EventSetup const&) {}
	std::shared_ptr<MLDQMPreprocessingStreamData> globalBeginRunSummary(edm::Run const&, edm::EventSetup const&) const;
	std::unique_ptr<MLDQMPreprocessingStreamData> beginStream(edm::StreamID) const;
	void streamBeginRun(edm::StreamID, edm::Run const&, edm::EventSetup const&) const;
	void analyze(edm::StreamID, const edm::Event&, const edm::EventSetup&) const;
	void streamEndRunSummary(edm::StreamID, edm::Run const&, edm::EventSetup const&, MLDQMPreprocessingStreamData*) const;
	void globalEndRunSummary(edm::Run const&, edm::EventSetup const&, MLDQMPreprocessingStreamData*) const;

	void endJob(){}
	void endRun(edm::Run const&, edm::EventSetup const&){}
	
	edm::Service<TFileService> _fs;

protected:
	edm::ESHandle<HcalDbService> _dbService;
	
 public:
	static const double adc2fC[];
};

#endif
