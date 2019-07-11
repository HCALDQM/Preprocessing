#include "HCALDQM/Preprocessing/interface/MLDQMPreprocessing.h"

#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/Framework/interface/ProductSelector.h"
#include "FWCore/Framework/interface/ProductSelectorRules.h"
#include "DataFormats/Provenance/interface/SelectedProducts.h"
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

#include <map>
#include "boost/foreach.hpp"
#include <TBranch.h>
#include <TLorentzVector.h>

// --------------------------------------------------------------------------------------------------------
// Updated using the example here:
//		https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/CalibTracker/SiStripCommon/plugins/ShallowTree.cc
//		https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/CalibTracker/SiStripCommon/interface/ShallowTree.h
// --------------------------------------------------------------------------------------------------------

const double MLDQMPreprocessing::adc2fC[] = {1.62, 4.86, 8.11, 11.35, 14.59, 17.84, 21.08, 24.32, 27.57, 30.81, 34.05, 37.30, 40.54, 43.78, 47.03, 50.27, 56.75, 63.24, 69.73, 76.21, 82.70, 89.19, 95.67, 102.2, 108.6, 115.1, 121.6, 128.1, 134.6, 141.1, 147.6, 154.0, 160.5, 167.0, 173.5, 180.0, 193.0, 205.9, 218.9, 231.9, 244.9, 257.8, 270.8, 283.8, 296.7, 309.7, 322.7, 335.7, 348.6, 361.6, 374.6, 387.6, 400.5, 413.5, 426.5, 439.4, 452.4, 478.4, 504.3, 530.3, 556.2, 582.1, 608.1, 634.0, 577.6, 603.0, 628.5, 654.0, 679.5, 705.0, 730.5, 756.0, 781.5, 806.9, 832.4, 857.9, 883.4, 908.9, 934.4, 959.9, 1010.9, 1061.8, 1112.8, 1163.8, 1214.8, 1265.7, 1316.7, 1367.7, 1418.7, 1469.6, 1520.6, 1571.6, 1622.6, 1673.5, 1724.5, 1775.5, 1826.5, 1877.5, 1928.4, 1979.4, 2081.4, 2183.3, 2285.3, 2387.2, 2489.2, 2591.1, 2693.1, 2795.0, 2897.0, 2998.9, 3100.9, 3202.8, 3304.8, 3406.8, 3508.7, 3610.7, 3712.6, 3814.6, 3916.5, 4018.5, 4120.4, 4324.3, 4528.2, 4732.1, 4936.1, 5140.0, 5343.9, 5547.8, 5331.9, 5542.5, 5753.1, 5963.7, 6174.3, 6384.9, 6595.5, 6806.1, 7016.7, 7227.3, 7437.9, 7648.4, 7859.0, 8069.6, 8280.2, 8490.8, 8912.0, 9333.2, 9754.3, 10175.5, 10596.7, 11017.9, 11439.1, 11860.3, 12281.4, 12702.6, 13123.8, 13545.0, 13966.2, 14387.3, 14808.5, 15229.7, 15650.9, 16072.1, 16493.2, 16914.4, 17756.8, 18599.1, 19441.5, 20283.9, 21126.2, 21968.6, 22811.0, 23653.3, 24495.7, 25338.0, 26180.4, 27022.8, 27865.1, 28707.5, 29549.9, 30392.2, 31234.6, 32076.9, 32919.3, 33761.7, 34604.0, 36288.8, 37973.5, 39658.2, 41342.9, 43027.6, 44712.4, 46397.1, 43321.6, 44990.1, 46658.6, 48327.1, 49995.7, 51664.2, 53332.7, 55001.2, 56669.7, 58338.2, 60006.7, 61675.2, 63343.7, 65012.3, 66680.8, 68349.3, 71686.3, 75023.3, 78360.3, 81697.4, 85034.4, 88371.4, 91708.4, 95045.4, 98382.5, 101719.5, 105056.5, 108393.5, 111730.6, 115067.6, 118404.6, 121741.6, 125078.6, 128415.7, 131752.7, 135089.7, 141763.8, 148437.8, 155111.8, 161785.9, 168459.9, 175134.0, 181808.0, 188482.1, 195156.1, 201830.1, 208504.2, 215178.2, 221852.3, 228526.3, 235200.4, 241874.4, 248548.4, 255222.5, 261896.5, 268570.6, 275244.6, 288592.7, 301940.8, 315288.9, 328637.0, 341985.1, 355333.1, 368681.2};

MLDQMPreprocessing::MLDQMPreprocessing(const edm::ParameterSet& ps) {
	_tagHBHE = ps.getUntrackedParameter<edm::InputTag>("tagHBHE", edm::InputTag("hcalDigis"));
	_tokHBHE = consumes<HBHEDigiCollection>(_tagHBHE);	

	_tagQIE11 = ps.getUntrackedParameter<edm::InputTag>("tagQIE11", edm::InputTag("hcalDigis"));
	_tokQIE11 = consumes<QIE11DigiCollection>(_tagQIE11);	

	_tagQIE10 = ps.getUntrackedParameter<edm::InputTag>("tagQIE10", edm::InputTag("hcalDigis"));
	_tokQIE10 = consumes<QIE10DigiCollection>(_tagQIE10);	

	_tagHO = ps.getUntrackedParameter<edm::InputTag>("tagHO", edm::InputTag("hcalDigis"));
	_tokHO = consumes<HODigiCollection>(_tagHO);	
}

MLDQMPreprocessing::~MLDQMPreprocessing() {
}

std::unique_ptr<MLDQMPreprocessingStreamData> MLDQMPreprocessing::beginStream(edm::StreamID sid) const {
	//return std::unique_ptr<MLDQMPreprocessing::ChannelHistogramMap>(new MLDQMPreprocessing::ChannelHistogramMap());
	MLDQMPreprocessingStreamData* data = new MLDQMPreprocessingStreamData();
	data->Init();
	return std::unique_ptr<MLDQMPreprocessingStreamData>(data);
}

std::shared_ptr<MLDQMPreprocessingStreamData> MLDQMPreprocessing::globalBeginRunSummary(edm::Run const& run, edm::EventSetup const& es) const {
	std::shared_ptr<MLDQMPreprocessingStreamData> returnData(new MLDQMPreprocessingStreamData());
	returnData->Init();
	return returnData;
}

void MLDQMPreprocessing::streamBeginRun(edm::StreamID sid, edm::Run const& run, edm::EventSetup const& es) const {
	// Book stream histograms
	streamCache(sid)->Init();
}

void MLDQMPreprocessing::beginJob() {
}

void MLDQMPreprocessing::analyze(edm::StreamID sid, const edm::Event& event, const edm::EventSetup& es) const {
	streamCache(sid)->_nevents += 1;

	// Load digis
	edm::Handle<HBHEDigiCollection> digis_hbhe;
	if (!event.getByToken(_tokHBHE, digis_hbhe)) {
		std::cout << "Collection HBHEDigiCollection isn't available" << _tagHBHE.label() << " " << _tagHBHE.instance() << std::endl;
		exit(1);
	}

	edm::Handle<QIE11DigiCollection> digis_qie11;
	if (!event.getByToken(_tokQIE11, digis_qie11)) {
		std::cout << "Collection QIE11DigiCollection isn't available" << _tagQIE11.label() << " " << _tagQIE11.instance() << std::endl;
		exit(1);
	}

	edm::Handle<QIE10DigiCollection> digis_qie10;
	if (!event.getByToken(_tokQIE10, digis_qie10)) {
		std::cout << "Collection QIE10DigiCollection isn't available" << _tagQIE10.label() << " " << _tagQIE10.instance() << std::endl;
		exit(1);
	}

	edm::Handle<HODigiCollection> digis_ho;
	if (!event.getByToken(_tokHO, digis_ho)) {
		std::cout << "Collection HODigiCollection isn't available" << _tagHO.label() << " " << _tagHO.instance() << std::endl;
		exit(1);
	}

	edm::ESHandle<HcalDbService> dbs;
	es.get<HcalDbRecord>().get(dbs);

	for (HBHEDigiCollection::const_iterator digi = digis_hbhe->begin(); digi != digis_hbhe->end(); ++digi) {
		HcalDetId const& did = digi->id();
		if (did.subdet() != HcalBarrel && did.subdet() != HcalEndcap) {
			continue;
		}
		CaloSamples digi_fC = hcaldqm::utilities::loadADC2fCDB<HBHEDataFrame>(dbs, did, *digi);
		double sumQ = hcaldqm::utilities::sumQDB<HBHEDataFrame>(dbs, digi_fC, did, *digi, 0, digi->size()-1);
		if (sumQ > 20) {
			streamCache(sid)->_hists["occupancy_20fC"]->Fill(did);
		}
		streamCache(sid)->_hists["sumQ"]->Fill(did, sumQ);
	}

	for (QIE11DigiCollection::const_iterator it = digis_qie11->begin(); it != digis_qie11->end(); ++it) {
		const QIE11DataFrame digi = static_cast<const QIE11DataFrame>(*it);

		HcalDetId const& did = digi.detid();
		if (did.subdet() != HcalBarrel && did.subdet() != HcalEndcap) {
			continue;
		}
		CaloSamples digi_fC = hcaldqm::utilities::loadADC2fCDB<QIE11DataFrame>(_dbService, did, digi);
		double sumQ = hcaldqm::utilities::sumQDB<QIE11DataFrame>(_dbService, digi_fC, did, digi, 0, digi.samples()-1);
		if (sumQ > 20) {
			streamCache(sid)->_hists["occupancy_20fC"]->Fill(did);
		}
		streamCache(sid)->_hists["sumQ"]->Fill(did, sumQ);
	}

	for (QIE10DigiCollection::const_iterator it = digis_qie10->begin(); it != digis_qie10->end(); ++it) {
		const QIE10DataFrame digi = static_cast<const QIE10DataFrame>(*it);

		HcalDetId const& did = digi.detid();
		if (did.subdet() != HcalForward) {
			continue;
		}
		CaloSamples digi_fC = hcaldqm::utilities::loadADC2fCDB<QIE10DataFrame>(_dbService, did, digi);
		double sumQ = hcaldqm::utilities::sumQDB<QIE10DataFrame>(_dbService, digi_fC, did, digi, 0, digi.samples()-1);
		if (sumQ > 20) {
			streamCache(sid)->_hists["occupancy_20fC"]->Fill(did);
		}
		streamCache(sid)->_hists["sumQ"]->Fill(did, sumQ);
	}

	for (HODigiCollection::const_iterator digi = digis_ho->begin(); digi != digis_ho->end(); ++digi) {
		HcalDetId const& did = digi->id();
		if (did.subdet() != HcalOuter) {
			continue;
		}
		CaloSamples digi_fC = hcaldqm::utilities::loadADC2fCDB<HODataFrame>(_dbService, did, *digi);
		double sumQ = hcaldqm::utilities::sumQDB<HODataFrame>(_dbService, digi_fC, did, *digi, 0, digi->size()-1);
		if (sumQ > 20) {
			streamCache(sid)->_hists["occupancy_20fC"]->Fill(did);
		}
		streamCache(sid)->_hists["sumQ"]->Fill(did, sumQ);
	}
}

void MLDQMPreprocessing::streamEndRunSummary(edm::StreamID sid, edm::Run const& run, edm::EventSetup const& es, MLDQMPreprocessingStreamData* globalData) const {
	//Add the values seen in this Stream to the total for this Run
	globalData->Add(streamCache(sid));
	streamCache(sid)->Reset();
}

void MLDQMPreprocessing::globalEndRunSummary(edm::Run const& run, edm::EventSetup const& es, MLDQMPreprocessingStreamData* globalData) const {
	globalData->Finish();
	_fs->cd();
	globalData->Write();
}
