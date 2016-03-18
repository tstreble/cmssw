///
/// \class L1TMuonGlobalParamsHelper
///
/// Description: Wrapper for L1TMuonGlobalParams
///
/// Implementation:
///
/// \author: Thomas Reis
///

#ifndef L1TMuonGlobalParamsHelper_h
#define L1TMuonGlobalParamsHelper_h

#include <memory>
#include <iostream>
#include <vector>

#include "CondFormats/L1TObjects/interface/L1TMuonGlobalParams.h"

class L1TMuonGlobalParamsHelper : public L1TMuonGlobalParams {

public:
  enum {absIsoCheckMem=0,
        relIsoCheckMem=1,
        idxSelMemPhi=2,
        idxSelMemEta=3,
        fwdPosSingleMatchQual=4,
        fwdNegSingleMatchQual=5,
        ovlPosSingleMatchQual=6,
        ovlNegSingleMatchQual=7,
        bOPosMatchQual=8,
        bONegMatchQual=9,
        fOPosMatchQual=10,
        fONegMatchQual=11,
        bPhiExtrapolation=12,
        oPhiExtrapolation=13,
        fPhiExtrapolation=14,
        bEtaExtrapolation=15,
        oEtaExtrapolation=16,
        fEtaExtrapolation=17,
        sortRank=18,
        FWVERSION=19,
        INPUTS_TO_DISABLE=20,
        NUM_GMTPARAMNODES=21
  };

  // string parameters indices
  enum spIdx {fname=0};

  // unsigned parameters indices
  enum upIdx {ptFactor=0, qualFactor=1, FWVERSION_IDX=0, CALOINPUTS_TO_DISABLE=0, BMTFINPUTS_TO_DISABLE=1, OMTFINPUTS_TO_DISABLE=2, EMTFINPUTS_TO_DISABLE=3};

  // double parameters indices
  enum dpIdx {maxdr=0, maxdrEtaFine=1};

  // input enable indices
  enum linkNr {CALOLINK1=8, EMTFPLINK1=36, OMTFPLINK1=42, BMTFLINK1=48, OMTFNLINK1=60, EMTFNLINK1=66}; // link numbers start at 0

  L1TMuonGlobalParamsHelper() { pnodes_.resize(NUM_GMTPARAMNODES); }
  L1TMuonGlobalParamsHelper(const L1TMuonGlobalParams);
  ~L1TMuonGlobalParamsHelper() {}

  // FW version
  unsigned fwVersion() const { return pnodes_[FWVERSION].uparams_.size() > FWVERSION_IDX ? pnodes_[FWVERSION].uparams_[FWVERSION_IDX] : 0; }
  void setFwVersion(unsigned fwVersion);

  // Input disables
  std::bitset<72> inputsToDisable() const;
  void setInputsToDisable(const std::bitset<72> &inputsToDisable); 
  void setCaloInputsToDisable(const std::bitset<28> &disables);
  void setEmtfpInputsToDisable(const std::bitset<6> &disables) { setEOmtfInputsToDisable(0, EMTFINPUTS_TO_DISABLE, disables); };
  void setOmtfpInputsToDisable(const std::bitset<6> &disables) { setEOmtfInputsToDisable(0, OMTFINPUTS_TO_DISABLE, disables); };
  void setBmtfInputsToDisable(const std::bitset<12> &disables);
  void setOmtfnInputsToDisable(const std::bitset<6> &disables) { setEOmtfInputsToDisable(6, OMTFINPUTS_TO_DISABLE, disables); };
  void setEmtfnInputsToDisable(const std::bitset<6> &disables) { setEOmtfInputsToDisable(6, EMTFINPUTS_TO_DISABLE, disables); };

  // LUTs
  l1t::LUT* absIsoCheckMemLUT()        { return &pnodes_[absIsoCheckMem].LUT_; }
  l1t::LUT* relIsoCheckMemLUT()        { return &pnodes_[relIsoCheckMem].LUT_; }
  l1t::LUT* idxSelMemPhiLUT()          { return &pnodes_[idxSelMemPhi].LUT_; }
  l1t::LUT* idxSelMemEtaLUT()          { return &pnodes_[idxSelMemEta].LUT_; }
  l1t::LUT* fwdPosSingleMatchQualLUT() { return &pnodes_[fwdPosSingleMatchQual].LUT_; }
  l1t::LUT* fwdNegSingleMatchQualLUT() { return &pnodes_[fwdNegSingleMatchQual].LUT_; }
  l1t::LUT* ovlPosSingleMatchQualLUT() { return &pnodes_[ovlPosSingleMatchQual].LUT_; }
  l1t::LUT* ovlNegSingleMatchQualLUT() { return &pnodes_[ovlNegSingleMatchQual].LUT_; }
  l1t::LUT* bOPosMatchQualLUT()        { return &pnodes_[bOPosMatchQual].LUT_; }
  l1t::LUT* bONegMatchQualLUT()        { return &pnodes_[bONegMatchQual].LUT_; }
  l1t::LUT* fOPosMatchQualLUT()        { return &pnodes_[fOPosMatchQual].LUT_; }
  l1t::LUT* fONegMatchQualLUT()        { return &pnodes_[fONegMatchQual].LUT_; }
  l1t::LUT* bPhiExtrapolationLUT()     { return &pnodes_[bPhiExtrapolation].LUT_; }
  l1t::LUT* oPhiExtrapolationLUT()     { return &pnodes_[oPhiExtrapolation].LUT_; }
  l1t::LUT* fPhiExtrapolationLUT()     { return &pnodes_[fPhiExtrapolation].LUT_; }
  l1t::LUT* bEtaExtrapolationLUT()     { return &pnodes_[bEtaExtrapolation].LUT_; }
  l1t::LUT* oEtaExtrapolationLUT()     { return &pnodes_[oEtaExtrapolation].LUT_; }
  l1t::LUT* fEtaExtrapolationLUT()     { return &pnodes_[fEtaExtrapolation].LUT_; }
  l1t::LUT* sortRankLUT()              { return &pnodes_[sortRank].LUT_; }
  void setAbsIsoCheckMemLUT        (const l1t::LUT & lut) { pnodes_[absIsoCheckMem].type_ = "LUT"; pnodes_[absIsoCheckMem].LUT_ = lut; }
  void setRelIsoCheckMemLUT        (const l1t::LUT & lut) { pnodes_[relIsoCheckMem].type_ = "LUT"; pnodes_[relIsoCheckMem].LUT_ = lut; }
  void setIdxSelMemPhiLUT          (const l1t::LUT & lut) { pnodes_[idxSelMemPhi].type_ = "LUT"; pnodes_[idxSelMemPhi].LUT_ = lut; }
  void setIdxSelMemEtaLUT          (const l1t::LUT & lut) { pnodes_[idxSelMemEta].type_ = "LUT"; pnodes_[idxSelMemEta].LUT_ = lut; }
  void setFwdPosSingleMatchQualLUT (const l1t::LUT & lut) { pnodes_[fwdPosSingleMatchQual].type_ = "LUT"; pnodes_[fwdPosSingleMatchQual].LUT_ = lut; }
  void setFwdNegSingleMatchQualLUT (const l1t::LUT & lut) { pnodes_[fwdNegSingleMatchQual].type_ = "LUT"; pnodes_[fwdNegSingleMatchQual].LUT_ = lut; }
  void setOvlPosSingleMatchQualLUT (const l1t::LUT & lut) { pnodes_[ovlPosSingleMatchQual].type_ = "LUT"; pnodes_[ovlPosSingleMatchQual].LUT_ = lut; }
  void setOvlNegSingleMatchQualLUT (const l1t::LUT & lut) { pnodes_[ovlNegSingleMatchQual].type_ = "LUT"; pnodes_[ovlNegSingleMatchQual].LUT_ = lut; }
  void setBOPosMatchQualLUT        (const l1t::LUT & lut) { pnodes_[bOPosMatchQual].type_ = "LUT"; pnodes_[bOPosMatchQual].LUT_ = lut; }
  void setBONegMatchQualLUT        (const l1t::LUT & lut) { pnodes_[bONegMatchQual].type_ = "LUT"; pnodes_[bONegMatchQual].LUT_ = lut; }
  void setFOPosMatchQualLUT        (const l1t::LUT & lut) { pnodes_[fOPosMatchQual].type_ = "LUT"; pnodes_[fOPosMatchQual].LUT_ = lut; }
  void setFONegMatchQualLUT        (const l1t::LUT & lut) { pnodes_[fONegMatchQual].type_ = "LUT"; pnodes_[fONegMatchQual].LUT_ = lut; }
  void setBPhiExtrapolationLUT     (const l1t::LUT & lut) { pnodes_[bPhiExtrapolation].type_ = "LUT"; pnodes_[bPhiExtrapolation].LUT_ = lut; }
  void setOPhiExtrapolationLUT     (const l1t::LUT & lut) { pnodes_[oPhiExtrapolation].type_ = "LUT"; pnodes_[oPhiExtrapolation].LUT_ = lut; }
  void setFPhiExtrapolationLUT     (const l1t::LUT & lut) { pnodes_[fPhiExtrapolation].type_ = "LUT"; pnodes_[fPhiExtrapolation].LUT_ = lut; }
  void setBEtaExtrapolationLUT     (const l1t::LUT & lut) { pnodes_[bEtaExtrapolation].type_ = "LUT"; pnodes_[bEtaExtrapolation].LUT_ = lut; }
  void setOEtaExtrapolationLUT     (const l1t::LUT & lut) { pnodes_[oEtaExtrapolation].type_ = "LUT"; pnodes_[oEtaExtrapolation].LUT_ = lut; }
  void setFEtaExtrapolationLUT     (const l1t::LUT & lut) { pnodes_[fEtaExtrapolation].type_ = "LUT"; pnodes_[fEtaExtrapolation].LUT_ = lut; }
  void setSortRankLUT              (const l1t::LUT & lut) { pnodes_[sortRank].type_ = "LUT"; pnodes_[sortRank].LUT_ = lut; }

  // LUT paths
  std::string absIsoCheckMemLUTPath() const        { return pnodes_[absIsoCheckMem].sparams_.size() > spIdx::fname ? pnodes_[absIsoCheckMem].sparams_[spIdx::fname] : ""; }
  std::string relIsoCheckMemLUTPath() const        { return pnodes_[relIsoCheckMem].sparams_.size() > spIdx::fname ? pnodes_[relIsoCheckMem].sparams_[spIdx::fname] : ""; }
  std::string idxSelMemPhiLUTPath() const          { return pnodes_[idxSelMemPhi].sparams_.size() > spIdx::fname ? pnodes_[idxSelMemPhi].sparams_[spIdx::fname] : ""; }
  std::string idxSelMemEtaLUTPath() const          { return pnodes_[idxSelMemEta].sparams_.size() > spIdx::fname ? pnodes_[idxSelMemEta].sparams_[spIdx::fname] : ""; }
  std::string fwdPosSingleMatchQualLUTPath() const { return pnodes_[fwdPosSingleMatchQual].sparams_.size() > spIdx::fname ? pnodes_[fwdPosSingleMatchQual].sparams_[spIdx::fname] : ""; }
  std::string fwdNegSingleMatchQualLUTPath() const { return pnodes_[fwdNegSingleMatchQual].sparams_.size() > spIdx::fname ? pnodes_[fwdNegSingleMatchQual].sparams_[spIdx::fname] : ""; }
  std::string ovlPosSingleMatchQualLUTPath() const { return pnodes_[ovlPosSingleMatchQual].sparams_.size() > spIdx::fname ? pnodes_[ovlPosSingleMatchQual].sparams_[spIdx::fname] : ""; }
  std::string ovlNegSingleMatchQualLUTPath() const { return pnodes_[ovlNegSingleMatchQual].sparams_.size() > spIdx::fname ? pnodes_[ovlNegSingleMatchQual].sparams_[spIdx::fname] : ""; }
  std::string bOPosMatchQualLUTPath() const        { return pnodes_[bOPosMatchQual].sparams_.size() > spIdx::fname ? pnodes_[bOPosMatchQual].sparams_[spIdx::fname] : ""; }
  std::string bONegMatchQualLUTPath() const        { return pnodes_[bONegMatchQual].sparams_.size() > spIdx::fname ? pnodes_[bONegMatchQual].sparams_[spIdx::fname] : ""; }
  std::string fOPosMatchQualLUTPath() const        { return pnodes_[fOPosMatchQual].sparams_.size() > spIdx::fname ? pnodes_[fOPosMatchQual].sparams_[spIdx::fname] : ""; }
  std::string fONegMatchQualLUTPath() const        { return pnodes_[fONegMatchQual].sparams_.size() > spIdx::fname ? pnodes_[fONegMatchQual].sparams_[spIdx::fname] : ""; }
  std::string bPhiExtrapolationLUTPath() const     { return pnodes_[bPhiExtrapolation].sparams_.size() > spIdx::fname ? pnodes_[bPhiExtrapolation].sparams_[spIdx::fname] : ""; }
  std::string oPhiExtrapolationLUTPath() const     { return pnodes_[oPhiExtrapolation].sparams_.size() > spIdx::fname ? pnodes_[oPhiExtrapolation].sparams_[spIdx::fname] : ""; }
  std::string fPhiExtrapolationLUTPath() const     { return pnodes_[fPhiExtrapolation].sparams_.size() > spIdx::fname ? pnodes_[fPhiExtrapolation].sparams_[spIdx::fname] : ""; }
  std::string bEtaExtrapolationLUTPath() const     { return pnodes_[bEtaExtrapolation].sparams_.size() > spIdx::fname ? pnodes_[bEtaExtrapolation].sparams_[spIdx::fname] : ""; }
  std::string oEtaExtrapolationLUTPath() const     { return pnodes_[oEtaExtrapolation].sparams_.size() > spIdx::fname ? pnodes_[oEtaExtrapolation].sparams_[spIdx::fname] : ""; }
  std::string fEtaExtrapolationLUTPath() const     { return pnodes_[fEtaExtrapolation].sparams_.size() > spIdx::fname ? pnodes_[fEtaExtrapolation].sparams_[spIdx::fname] : ""; }
  std::string sortRankLUTPath() const              { return pnodes_[sortRank].sparams_.size() > spIdx::fname ? pnodes_[sortRank].sparams_[spIdx::fname] : ""; }
  void setAbsIsoCheckMemLUTPath        (std::string path) { pnodes_[absIsoCheckMem].sparams_.push_back(path); }
  void setRelIsoCheckMemLUTPath        (std::string path) { pnodes_[relIsoCheckMem].sparams_.push_back(path); }
  void setIdxSelMemPhiLUTPath          (std::string path) { pnodes_[idxSelMemPhi].sparams_.push_back(path); }
  void setIdxSelMemEtaLUTPath          (std::string path) { pnodes_[idxSelMemEta].sparams_.push_back(path); }
  void setFwdPosSingleMatchQualLUTPath (std::string path) { pnodes_[fwdPosSingleMatchQual].sparams_.push_back(path); }
  void setFwdNegSingleMatchQualLUTPath (std::string path) { pnodes_[fwdNegSingleMatchQual].sparams_.push_back(path); }
  void setOvlPosSingleMatchQualLUTPath (std::string path) { pnodes_[ovlPosSingleMatchQual].sparams_.push_back(path); }
  void setOvlNegSingleMatchQualLUTPath (std::string path) { pnodes_[ovlNegSingleMatchQual].sparams_.push_back(path); }
  void setBOPosMatchQualLUTPath        (std::string path) { pnodes_[bOPosMatchQual].sparams_.push_back(path); }
  void setBONegMatchQualLUTPath        (std::string path) { pnodes_[bONegMatchQual].sparams_.push_back(path); }
  void setFOPosMatchQualLUTPath        (std::string path) { pnodes_[fOPosMatchQual].sparams_.push_back(path); }
  void setFONegMatchQualLUTPath        (std::string path) { pnodes_[fONegMatchQual].sparams_.push_back(path); }
  void setBPhiExtrapolationLUTPath     (std::string path) { pnodes_[bPhiExtrapolation].sparams_.push_back(path); }
  void setOPhiExtrapolationLUTPath     (std::string path) { pnodes_[oPhiExtrapolation].sparams_.push_back(path); }
  void setFPhiExtrapolationLUTPath     (std::string path) { pnodes_[fPhiExtrapolation].sparams_.push_back(path); }
  void setBEtaExtrapolationLUTPath     (std::string path) { pnodes_[bEtaExtrapolation].sparams_.push_back(path); }
  void setOEtaExtrapolationLUTPath     (std::string path) { pnodes_[oEtaExtrapolation].sparams_.push_back(path); }
  void setFEtaExtrapolationLUTPath     (std::string path) { pnodes_[fEtaExtrapolation].sparams_.push_back(path); }
  void setSortRankLUTPath              (std::string path) { pnodes_[sortRank].sparams_.push_back(path); }

  // Cancel out LUT max dR
  double fwdPosSingleMatchQualLUTMaxDR() const { return pnodes_[fwdPosSingleMatchQual].dparams_.size() > dpIdx::maxdr ? pnodes_[fwdPosSingleMatchQual].dparams_[dpIdx::maxdr] : 0.; }
  double fwdNegSingleMatchQualLUTMaxDR() const { return pnodes_[fwdNegSingleMatchQual].dparams_.size() > dpIdx::maxdr ? pnodes_[fwdNegSingleMatchQual].dparams_[dpIdx::maxdr] : 0.; }
  double ovlPosSingleMatchQualLUTMaxDR() const { return pnodes_[ovlPosSingleMatchQual].dparams_.size() > dpIdx::maxdr ? pnodes_[ovlPosSingleMatchQual].dparams_[dpIdx::maxdr] : 0.; }
  double ovlNegSingleMatchQualLUTMaxDR() const { return pnodes_[ovlNegSingleMatchQual].dparams_.size() > dpIdx::maxdr ? pnodes_[ovlNegSingleMatchQual].dparams_[dpIdx::maxdr] : 0.; }
  double bOPosMatchQualLUTMaxDR() const        { return pnodes_[bOPosMatchQual].dparams_.size() > dpIdx::maxdr ? pnodes_[bOPosMatchQual].dparams_[dpIdx::maxdr] : 0.; }
  double bONegMatchQualLUTMaxDR() const        { return pnodes_[bONegMatchQual].dparams_.size() > dpIdx::maxdr ? pnodes_[bONegMatchQual].dparams_[dpIdx::maxdr] : 0.; }
  double bOPosMatchQualLUTMaxDREtaFine() const { return pnodes_[bOPosMatchQual].dparams_.size() > dpIdx::maxdrEtaFine ? pnodes_[bOPosMatchQual].dparams_[dpIdx::maxdrEtaFine] : 0.; }
  double bONegMatchQualLUTMaxDREtaFine() const { return pnodes_[bONegMatchQual].dparams_.size() > dpIdx::maxdrEtaFine ? pnodes_[bONegMatchQual].dparams_[dpIdx::maxdrEtaFine] : 0.; }
  double fOPosMatchQualLUTMaxDR() const        { return pnodes_[fOPosMatchQual].dparams_.size() > dpIdx::maxdr ? pnodes_[fOPosMatchQual].dparams_[dpIdx::maxdr] : 0.; }
  double fONegMatchQualLUTMaxDR() const        { return pnodes_[fONegMatchQual].dparams_.size() > dpIdx::maxdr ? pnodes_[fONegMatchQual].dparams_[dpIdx::maxdr] : 0.; }
  void setFwdPosSingleMatchQualLUTMaxDR (double maxDR) { pnodes_[fwdPosSingleMatchQual].dparams_.push_back(maxDR); }
  void setFwdNegSingleMatchQualLUTMaxDR (double maxDR) { pnodes_[fwdNegSingleMatchQual].dparams_.push_back(maxDR); }
  void setOvlPosSingleMatchQualLUTMaxDR (double maxDR) { pnodes_[ovlPosSingleMatchQual].dparams_.push_back(maxDR); }
  void setOvlNegSingleMatchQualLUTMaxDR (double maxDR) { pnodes_[ovlNegSingleMatchQual].dparams_.push_back(maxDR); }
  void setBOPosMatchQualLUTMaxDR        (double maxDR, double maxDREtaFine) { pnodes_[bOPosMatchQual].dparams_.push_back(maxDR); pnodes_[bOPosMatchQual].dparams_.push_back(maxDREtaFine); }
  void setBONegMatchQualLUTMaxDR        (double maxDR, double maxDREtaFine) { pnodes_[bONegMatchQual].dparams_.push_back(maxDR); pnodes_[bONegMatchQual].dparams_.push_back(maxDREtaFine); }
  void setFOPosMatchQualLUTMaxDR        (double maxDR) { pnodes_[fOPosMatchQual].dparams_.push_back(maxDR); }
  void setFONegMatchQualLUTMaxDR        (double maxDR) { pnodes_[fONegMatchQual].dparams_.push_back(maxDR); }

  // Sort rank LUT factors for pT and quality
  unsigned sortRankLUTPtFactor() const   { return pnodes_[sortRank].uparams_.size() > upIdx::ptFactor ? pnodes_[sortRank].uparams_[upIdx::ptFactor] : 0; }
  unsigned sortRankLUTQualFactor() const { return pnodes_[sortRank].uparams_.size() > upIdx::qualFactor ? pnodes_[sortRank].uparams_[upIdx::qualFactor] : 0; }
  void setSortRankLUTFactors(unsigned ptFactor, unsigned qualFactor) { pnodes_[sortRank].uparams_.push_back(ptFactor); pnodes_[sortRank].uparams_.push_back(qualFactor); }

  // print parameters to stream:
  void print(std::ostream&) const;
  friend std::ostream& operator<<(std::ostream& o, const L1TMuonGlobalParamsHelper & p) { p.print(o); return o; }

private:

  std::bitset<6> eomtfInputsToDisable(const int &startLink);
  void setEOmtfInputsToDisable(const size_t &startIdx, const int &tfIdx, const std::bitset<6> &disables);

};
#endif
