#include "TString.h"
#include "TFile.h"
#include "TList.h"

#include "Lattice.h"
#include "AnalysisLattice.h"

ClassImp(AnalysisLattice)

AnalysisLattice::AnalysisLattice() {}

AnalysisLattice::AnalysisLattice(const TFile& f) : file_name(f.GetName()) {
  TList * li = f.GetListOfKeys();
  for(int i = 0; i < li->GetEntries(); i++){
    l.push_back(li->At(i)->GetName());
  }
}

AnalysisLattice::AnalysisLattice(const TString& fname) : file_name(fname) {
  TFile f(file_name);
  TList * li = f.GetListOfKeys();
  for(int i = 0; i < li->GetEntries(); i++){
    l.push_back(li->At(i)->GetName());
  }
}

TString AnalysisLattice::getFileName() const { return file_name; }

vector<string> AnalysisLattice::getList() const { return l; }
