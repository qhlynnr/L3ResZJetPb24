#include "TSystem.h"
#include <iostream>

void compile(bool clean = true) {

  if (clean) {
    std::cout << "Cleaning pre-existing compilation artifacts..." << std::endl;
    
    // Remove histograms.C artifacts
    gSystem->Unlink("histograms_C.so");
    gSystem->Unlink("histograms_C.d");
    gSystem->Unlink("histograms_C_ACLiC_dict_rdict.pcm");
    
    // Remove eventhistograms.C artifacts
    gSystem->Unlink("eventhistograms_C.so");
    gSystem->Unlink("eventhistograms_C.d");
    gSystem->Unlink("eventhistograms_C_ACLiC_dict_rdict.pcm");

    // Remove analyse_PhotonJet.cc artifacts. This one embeds a dictionary
    // that references histograms/eventhistograms classes, so its cached .so
    // must always be rebuilt alongside them - otherwise it goes stale
    // relative to the new histograms_C.so/eventhistograms_C.so and every
    // job loading it crashes (cling dictionary redefinition -> DeclUnloader
    // segfault, root exit status 129).
    gSystem->Unlink("analyse_PhotonJet_cc.so");
    gSystem->Unlink("analyse_PhotonJet_cc.d");
    gSystem->Unlink("analyse_PhotonJet_cc_ACLiC_dict_rdict.pcm");

    std::cout << "Cleanup complete." << std::endl;
  }

  std::cout << "Compiling histograms.C..." << std::endl;
  gROOT->ProcessLine(".L histograms.C+");

  std::cout << "Compiling eventhistograms.C..." << std::endl;
  gROOT->ProcessLine(".L eventhistograms.C+");

  // Must come after histograms.C/eventhistograms.C so its dictionary is
  // built against their current versions.
  std::cout << "Compiling analyse_PhotonJet.cc..." << std::endl;
  gROOT->ProcessLine(".L analyse_PhotonJet.cc+");

  std::cout << "Compilation complete!" << std::endl;
  gROOT->ProcessLine(".q");

}
