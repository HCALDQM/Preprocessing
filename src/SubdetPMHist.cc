#ifndef SubdetPMHist_cc
#define SubdetPMHist_cc

#include "HCALDQM/Preprocessing/interface/SubdetPMHist.h"

const std::vector<SubdetPMHist::SubdetPM_t> SubdetPMHist::_subdetpm_list = {kHBM, kHBP, kHEM, kHEP, kHOM, kHOP, kHFM, kHFP};
const std::map<SubdetPMHist::SubdetPM_t, TString> SubdetPMHist::_subdetpm_names = {
	{kHBM, "HBM"},
	{kHBP, "HBP"},
	{kHEM, "HEM"},
	{kHEP, "HEP"},
	{kHOM, "HOM"},
	{kHOP, "HOP"},
	{kHFM, "HFM"},
	{kHFP, "HFP"},
};
const std::map<SubdetPMHist::SubdetPM_t, std::pair<int, int>> SubdetPMHist::_subdetpm_ieta = {
	{kHBM, std::pair(-16, -1)},
	{kHBP, std::pair(1, 16)},
	{kHEM, std::pair(-29, -16)},
	{kHEP, std::pair(16, 29)},
	{kHOM, std::pair(-16, -1)},
	{kHOP, std::pair(1, 16)},
	{kHFM, std::pair(-41, -29)},
	{kHFP, std::pair(29, 41)}
};

SubdetPMHist::SubdetPMHist(TString name, TString title) {
	for (auto& it : _subdetpm_list) {
		TString subdetpm_name = _subdetpm_names.at(it);
		TString hist_name = name + "_" + subdetpm_name;
		TString hist_title = title + "_" + subdetpm_name;
		Int_t ieta_min = _subdetpm_ieta.at(it).first;
		Int_t ieta_max = _subdetpm_ieta.at(it).second;
		_hists[it] = new TH2D(hist_name, hist_title, 
						ieta_max - ieta_min + 1, ieta_min - 0.5, ieta_max + 0.5, 
						72, 0.5, 72.5);
	}
}

SubdetPMHist::~SubdetPMHist() {
	_hists.clear();
}

SubdetPMHist::SubdetPM_t SubdetPMHist::GetSubdetPM(const HcalDetId hdid) {
	SubdetPM_t subdetpm = kNULL;
	if (hdid.subdet() == HcalBarrel) {
		if (hdid.ieta() > 0) {
			subdetpm = kHBP;
		} else {
			subdetpm = kHBM;
		}
	} else if (hdid.subdet() == HcalEndcap) {
		if (hdid.ieta() > 0) {
			subdetpm = kHEP;
		} else {
			subdetpm = kHEM;
		}
	} else if (hdid.subdet() == HcalOuter) {
		if (hdid.ieta() > 0) {
			subdetpm = kHOP;
		} else {
			subdetpm = kHOM;
		}
	} else if (hdid.subdet() == HcalForward) {
		if (hdid.ieta() > 0) {
			subdetpm = kHFP;
		} else {
			subdetpm = kHFM;
		}
	}
	return subdetpm;
}

void SubdetPMHist::Fill(const HcalDetId hdid, double w) {
	SubdetPM_t subdetpm = GetSubdetPM(hdid);
	if (subdetpm == kNULL) {
		std::cerr << "Attempt to fill SubdetPMHist with unknown HcalSubdetector: " << hdid.subdet() << std::endl;
		exit(1);
	}
	_hists[subdetpm]->Fill(hdid.ieta(), hdid.iphi(), w);
}

void SubdetPMHist::Add(const SubdetPMHist* hother) {
	for (auto& it : _subdetpm_list) {
		_hists[it]->Add(hother->at(it));
	}
}

void SubdetPMHist::Scale(const double scale) {
	for (auto& it : _subdetpm_list) {
		_hists[it]->Scale(scale);
	}
}

void SubdetPMHist::Reset() {
	for (auto& it : _subdetpm_list) {
		_hists[it]->Reset();
	}
}

void SubdetPMHist::Write() {
	for (auto& it : _subdetpm_list) {
		_hists[it]->Write();
	}
}


/*
std::pair<int, int> SubdetPMHist::GetBins(const HcalDetId hdid) {
	int bin_ieta = hdid.ieta() - _subdetpm_ieta[it].first + 1;
	int bin_iphi = hdid.iphi();
	return std::pair(bin_ieta, bin_iphi);
}
*/

#endif