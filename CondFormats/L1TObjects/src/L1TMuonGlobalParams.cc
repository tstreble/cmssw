#include "CondFormats/L1TObjects/interface/L1TMuonGlobalParams.h"

std::bitset<72> L1TMuonGlobalParams::inputsToDisable() const
{
  std::bitset<72> inputsToDisable;
  if (pnodes_[INPUTS_TO_DISABLE].uparams_.size() < 4) {
    return inputsToDisable;
  }

  for (size_t i = 0; i < 28; ++i) {
    inputsToDisable[CALOLINK1 + i] = ((pnodes_[INPUTS_TO_DISABLE].uparams_[CALOINPUTS_TO_DISABLE] >> i) & 0x1);
    if (i < CALOLINK1) {
      // disable unused inputs
      inputsToDisable[i] = 0x1;
    }
    if (i < 12) {
      inputsToDisable[BMTFLINK1 + i] = ((pnodes_[INPUTS_TO_DISABLE].uparams_[BMTFINPUTS_TO_DISABLE] >> i) & 0x1);
      if (i < 6) {
        inputsToDisable[EMTFPLINK1 + i] = ((pnodes_[INPUTS_TO_DISABLE].uparams_[EMTFINPUTS_TO_DISABLE] >> i) & 0x1);
        inputsToDisable[OMTFPLINK1 + i] = ((pnodes_[INPUTS_TO_DISABLE].uparams_[OMTFINPUTS_TO_DISABLE] >> i) & 0x1);
        inputsToDisable[OMTFNLINK1 + i] = ((pnodes_[INPUTS_TO_DISABLE].uparams_[OMTFINPUTS_TO_DISABLE] >> (i + 6)) & 0x1);
        inputsToDisable[EMTFNLINK1 + i] = ((pnodes_[INPUTS_TO_DISABLE].uparams_[EMTFINPUTS_TO_DISABLE] >> (i + 6)) & 0x1);
      }
    }
  }
  return inputsToDisable;
}


void L1TMuonGlobalParams::setFwVersion(unsigned fwVersion)
{
  pnodes_[FWVERSION].uparams_.resize(1);
  pnodes_[FWVERSION].uparams_[FWVERSION_IDX] = fwVersion;
}


void L1TMuonGlobalParams::setBxMin(int bxMin)
{
  pnodes_[BXRANGE].iparams_.resize(2);
  pnodes_[BXRANGE].iparams_[BXMIN] = bxMin;
}


void L1TMuonGlobalParams::setBxMax(int bxMax)
{
  pnodes_[BXRANGE].iparams_.resize(2);
  pnodes_[BXRANGE].iparams_[BXMAX] = bxMax;
}


void L1TMuonGlobalParams::setInputsToDisable(const std::bitset<72> &inputsToDisable)
{
  pnodes_[INPUTS_TO_DISABLE].uparams_.resize(4);
  for (size_t i = 0; i < 28; ++i) {
    pnodes_[INPUTS_TO_DISABLE].uparams_[CALOINPUTS_TO_DISABLE] += (inputsToDisable.test(CALOLINK1 + i) << i);
    if (i < 12) {
      pnodes_[INPUTS_TO_DISABLE].uparams_[BMTFINPUTS_TO_DISABLE] += (inputsToDisable.test(BMTFLINK1 + i) << i);
      if (i < 6) {
        pnodes_[INPUTS_TO_DISABLE].uparams_[OMTFINPUTS_TO_DISABLE] += (inputsToDisable.test(OMTFPLINK1 + i) << i);
        pnodes_[INPUTS_TO_DISABLE].uparams_[OMTFINPUTS_TO_DISABLE] += (inputsToDisable.test(OMTFNLINK1 + i) << (i + 6));
        pnodes_[INPUTS_TO_DISABLE].uparams_[EMTFINPUTS_TO_DISABLE] += (inputsToDisable.test(EMTFPLINK1 + i) << i);
        pnodes_[INPUTS_TO_DISABLE].uparams_[EMTFINPUTS_TO_DISABLE] += (inputsToDisable.test(EMTFNLINK1 + i) << (i + 6));
      }
    }
  }
}


void L1TMuonGlobalParams::setCaloInputsToDisable(const std::bitset<28> &inputsToDisable)
{
  pnodes_[INPUTS_TO_DISABLE].uparams_.resize(4);
  for (size_t i = 0; i < 28; ++i) {
    pnodes_[INPUTS_TO_DISABLE].uparams_[CALOINPUTS_TO_DISABLE] += (inputsToDisable.test(i) << i);
  }
}


void L1TMuonGlobalParams::setEOmtfInputsToDisable(const size_t &startIdx, const int &tfIdx, const std::bitset<6> &inputsToDisable)
{
  pnodes_[INPUTS_TO_DISABLE].uparams_.resize(4);
  for (size_t i = 0; i < 6; ++i) {
    pnodes_[INPUTS_TO_DISABLE].uparams_[tfIdx] += (inputsToDisable.test(i) << (i + startIdx));
  }
}


void L1TMuonGlobalParams::setBmtfInputsToDisable(const std::bitset<12> &inputsToDisable)
{
  pnodes_[INPUTS_TO_DISABLE].uparams_.resize(4);
  for (size_t i = 0; i < 12; ++i) {
    pnodes_[INPUTS_TO_DISABLE].uparams_[BMTFINPUTS_TO_DISABLE] += (inputsToDisable.test(i) << i);
  }
}


void L1TMuonGlobalParams::print(std::ostream& out) const {

  out << "L1 MicroGMT Parameters" << std::endl;

  out << "Firmware version: " << this->fwVersion() << std::endl;

  out << "Output BX range from " << this->bxMin() << " to " << this->bxMax() << std::endl;

  out << "InputsToDisable: " << this->inputsToDisable() << std::endl;
  out << "                 EMTF-|OMTF-|   BMTF    |OMTF+|EMTF+|            CALO           |  res  0" << std::endl;

  out << "LUT paths (LUTs are generated analytically if path is empty)" << std::endl;
  out << " Abs isolation checkMem LUT path: "        << this->absIsoCheckMemLUTPath() << std::endl;
  out << " Rel isolation checkMem LUT path: "        << this->relIsoCheckMemLUTPath() << std::endl;
  out << " Index selMem phi LUT path: "              << this->idxSelMemPhiLUTPath() << std::endl;
  out << " Index selMem eta LUT path: "              << this->idxSelMemEtaLUTPath() << std::endl;
  out << " Forward pos MatchQual LUT path: "         << this->fwdPosSingleMatchQualLUTPath() << ", max dR (Used when LUT path empty): " << this->fwdPosSingleMatchQualLUTMaxDR() << std::endl;
  out << " Forward neg MatchQual LUT path: "         << this->fwdNegSingleMatchQualLUTPath() << ", max dR (Used when LUT path empty): " << this->fwdNegSingleMatchQualLUTMaxDR() << std::endl;
  out << " Overlap pos MatchQual LUT path: "         << this->ovlPosSingleMatchQualLUTPath() << ", max dR (Used when LUT path empty): " << this->ovlPosSingleMatchQualLUTMaxDR() << std::endl;
  out << " Overlap neg MatchQual LUT path: "         << this->ovlNegSingleMatchQualLUTPath() << ", max dR (Used when LUT path empty): " << this->ovlNegSingleMatchQualLUTMaxDR() << std::endl;
  out << " Barrel-Overlap pos MatchQual LUT path: "  << this->bOPosMatchQualLUTPath() << ", max dR (Used when LUT path empty): " << this->bOPosMatchQualLUTMaxDR() << ", max dR when eta-fine bit set: " << this->bOPosMatchQualLUTMaxDREtaFine() << std::endl;
  out << " Barrel-Overlap neg MatchQual LUT path: "  << this->bONegMatchQualLUTPath() << ", max dR (Used when LUT path empty): " << this->bONegMatchQualLUTMaxDR() << ", max dR when eta-fine bit set: " << this->bONegMatchQualLUTMaxDREtaFine() << std::endl;
  out << " Forward-Overlap pos MatchQual LUT path: " << this->fOPosMatchQualLUTPath() << ", max dR (Used when LUT path empty): " << this->fOPosMatchQualLUTMaxDR() << std::endl;
  out << " Forward-Overlap neg MatchQual LUT path: " << this->fONegMatchQualLUTPath() << ", max dR (Used when LUT path empty): " << this->fONegMatchQualLUTMaxDR() << std::endl;
  out << " Barrel phi extrapolation LUT path: "      << this->bPhiExtrapolationLUTPath() << std::endl;
  out << " Overlap phi extrapolation LUT path: "     << this->oPhiExtrapolationLUTPath() << std::endl;
  out << " Forward phi extrapolation LUT path: "     << this->fPhiExtrapolationLUTPath() << std::endl;
  out << " Barrel eta extrapolation LUT path: "      << this->bEtaExtrapolationLUTPath() << std::endl;
  out << " Overlap eta extrapolation LUT path: "     << this->oEtaExtrapolationLUTPath() << std::endl;
  out << " Forward eta extrapolation LUT path: "     << this->fEtaExtrapolationLUTPath() << std::endl;
  out << " Sort rank LUT path: "                     << this->sortRankLUTPath() << ", pT and quality factors (Used when LUT path empty): pT factor: " << this->sortRankLUTPtFactor() << ", quality factor: " << this->sortRankLUTQualFactor() << std::endl;
}
