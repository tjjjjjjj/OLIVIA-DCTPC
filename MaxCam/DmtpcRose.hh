#ifndef DMTPC_ROSE_HH
#define DMTPC_ROSE_HH

#include "TROOT.h"

// Plotting rose histograms
// thanks to Asher

class TObjArray;
class TH1;
class TString;

/** Way to make rose histograms for plotting data \author Asher Kaboth*/

namespace DmtpcRose {
   /**
      Draws a rose histogram. Options are SAME = draw on the same plot.
      \param slices a TObjArray from makeRoseSlices containing the information for the rose plot
      \param opt a TString for the options
   */
  void DrawRose(TObjArray* slices, TString opt="");
   /**
      Makes an array of TCrowns in preparation for plotting a rose histogram
      \param h histogram IN RADIANS to be plotted.
   */
  TObjArray* makeRoseSlices(TH1* h);

};

#endif

