{
  //   gSystem->Setenv("MCTABLES","/Users/ataralas/Documents/ReferenceDocuments/GradSchool/FisherResearch/MonteCarlo/projects/DarkMatter/MaxCam/tables");
   gSystem->Setenv("MCTABLES","../../MaxCam/tables");

  gSystem->Load("$ROOTSYS/lib/libPhysics.so");
  gSystem->Load("$ROOTSYS/lib/libGraf.so");

// do we need this? 
//  gSystem->Load("/usr/local/lib/libcfitsio.so");

// Load MaxCam libraries 
// ======================
// To use production copy of MaxCam libs 
//  gSystem->Load("/usr/local/lib/MaxCam.so"); // production copy 

// To use your own copy of MaxCam libs 
  gSystem->Load("../../MaxCam/MaxCam_macosx.so");

// Load RooFit
  gSystem->Load("libRooFit.so"); 

  gStyle->SetPalette(1);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleColor(0);
  gStyle->SetTitleFont(42);
  gStyle->SetStatColor(0);
  //gStyle->SetFillColor(0);
  gStyle->SetTitleColor(1);

  // set the paper & margin sizes
  //gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.2); // 0.2
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadLeftMargin(0.2);

  gStyle->SetNdivisions(505,"x");
  gStyle->SetNdivisions(505,"y");

  // use large Times-Roman fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.08);

  gStyle->SetLabelFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleFont(42,"z");

  gStyle->SetTitleOffset(1.25, "x");
  gStyle->SetTitleOffset(1.25, "y");
  gStyle->SetTitleOffset(1.25, "z");

  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.06,"z");

  // stat box
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.95);
  gStyle->SetStatW(.200);
  gStyle->SetStatH(.125);
  gStyle->SetStatColor(0);
  gStyle->SetStatStyle(0);
  gStyle->SetStatFont(42);

  // title
  gStyle->SetTitleX(0.3);
  gStyle->SetTitleW(0.5);


  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetHistLineWidth(1.85);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);


  gStyle->SetPalette(1);
}

