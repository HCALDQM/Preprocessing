#ifndef SubdetPMHist_h
#define SubdetPMHist_h

#include <iostream>
#include <map>
#include <vector>

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include <TH2D.h>
#include <TString.h>

class SubdetPMHist {
public:
	enum SubdetPM_t {
		kNULL,
		kHBM,
		kHBP,
		kHEM,
		kHEP,
		kHOM,
		kHOP,
		kHFM,
		kHFP,
	};
public:
	SubdetPMHist(TString name, TString title);
	~SubdetPMHist();

	inline TH2D* operator [] (SubdetPM_t subdet) const {
		return _hists.at(subdet);
	}

	inline TH2D* at(SubdetPM_t subdet) const {
		return _hists.at(subdet);
	}

	void Fill(const HcalDetId hdid, Double_t w=1.);

	void Add(const SubdetPMHist* hother);

	void Scale(const double scale);

	void Reset();

	void Write();
	//std::pair<int, int> GetBins(const HcalDetId hdid);

private:
	SubdetPM_t GetSubdetPM(const HcalDetId hdid);


public:

private:


	std::map<SubdetPM_t, TH2D*> _hists;	

	static const std::vector<SubdetPM_t> _subdetpm_list;
	static const std::map<SubdetPM_t, TString> _subdetpm_names;
	static const std::map<SubdetPM_t, std::pair<int, int>> _subdetpm_ieta;

};

#endif