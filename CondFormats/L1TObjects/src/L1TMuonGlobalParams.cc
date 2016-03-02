#include "CondFormats/L1TObjects/interface/L1TMuonGlobalParams.h"

std::bitset<32> L1TMuonGlobalParams::caloInputEnables()
{
  std::bitset<32> subset;
  size_t sum = 0;
  for (int link = CALOLINK1; link < CALOLINK1 + 32; ++link, ++sum) {
    subset.set(sum, inputEnables_.test(link-1));
  }
  return subset;
}


std::bitset<6> L1TMuonGlobalParams::eomtfInputEnables(const int &startLink)
{
  std::bitset<6> subset;
  size_t proc = 0;
  for (int link = startLink; link < startLink + 6; ++link, ++proc) {
    subset.set(proc, inputEnables_.test(link-1));
  }
  return subset;
}


std::bitset<12> L1TMuonGlobalParams::bmtfInputEnables()
{
  std::bitset<12> subset;
  size_t proc = 0;
  for (int link = BMTFLINK1; link < BMTFLINK1 + 12; ++link, ++proc) {
    subset.set(proc, inputEnables_.test(link-1));
  }
  return subset;
}


void L1TMuonGlobalParams::setCaloInputEnables(const std::bitset<32> &enables)
{
  for (size_t sum = 0; sum < 32; ++sum) {
    inputEnables_.set(CALOLINK1 + sum, enables[sum]);
  }
}


void L1TMuonGlobalParams::setEOmtfInputEnables(const int &startLink, const std::bitset<6> &enables)
{
  for (size_t proc = 0; proc < 6; ++proc) {
    inputEnables_.set(startLink + proc, enables[proc]);
  }
}


void L1TMuonGlobalParams::setBmtfInputEnables(const std::bitset<12> &enables)
{
  for (size_t proc = 0; proc < 12; ++proc) {
    inputEnables_.set(BMTFLINK1 + proc, enables[proc]);
  }
}


void L1TMuonGlobalParams::print(std::ostream& out) const {

  out << "L1 MicroGMT Parameters" << std::endl;

  out << "Firmware version: " << fwVersion_ << std::endl;

  out << "Output BX range from " << bxMin_ << " to " << bxMax_ << std::endl;

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
