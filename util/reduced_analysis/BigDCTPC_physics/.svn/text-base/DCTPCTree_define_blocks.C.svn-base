#define DCTPC_runtree_cxx
#include "../input_files/DCTPC_runtree.h"
#include "../input_files/Calibration_constants_BigDCTPC.h";
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "TVirtualFitter.h"
#include <vector>

//This macro is to define blocks of 50-99 consecutive runs within a sequence, minimizing leftover runs

void DCTPC_runtree::Loop()
{

TStopwatch timer;
gROOT->SetStyle("Plain");
gStyle->SetEndErrorSize(3);
gStyle->SetPalette(1,0);
gStyle->SetLineWidth(2);
gStyle->SetHistLineWidth(2);
gStyle->SetOptStat(kFALSE);
//gStyle->SetOptFit(kFALSE);
gStyle->SetOptFit(1);
TH1::AddDirectory(false);
  
TFile *outtree = new TFile("$BigDCTPC_physics_input");
TTree *dctreepc = (TTree*)outtree->Get("dctpc_runinfo");
DCTPC_runtree aStep(dctreepc);
Long64_t nentries = dctreepc->GetEntries();

////////////////////////////////////////////////////////////////
//cout<<"Determine the first and last run and sequence"<<endl;
////////////////////////////////////////////////////////////////

int first_run;
int first_seq;
int last_seq;

for (int event = 0; event<nentries; event++){
	aStep.GetEntry(event); 
	if(event == 0){
		first_run = aStep.RunNum;
		first_seq = aStep.SequenceNum;
		last_seq = aStep.SequenceNum;
		}
	if (last_seq != aStep.SequenceNum){
		last_seq = aStep.SequenceNum;
		}
	}

////////////////////////////////////////////////////////////////
//cout<<"Count the number of runs in each sequence"<<endl;
////////////////////////////////////////////////////////////////

int runs_this_seq = 1;
int curr_run = first_run;
int curr_seq = first_seq;
std::vector<int> runs_in_seq(last_seq);

for (int event = 0; event<nentries; event++){
	aStep.GetEntry(event);
	if (curr_seq == aStep.SequenceNum && curr_run != aStep.RunNum){
		runs_this_seq++;
		curr_run = aStep.RunNum;
		}
	if (curr_seq != aStep.SequenceNum){
		runs_in_seq[curr_seq] = runs_this_seq;
		curr_seq = aStep.SequenceNum;
		runs_this_seq = 1;
		curr_run = aStep.RunNum;
		}
	}
runs_in_seq[curr_seq] = runs_this_seq;	

// for (curr_seq = 0; curr_seq<last_seq+1; curr_seq++){
// 	cout<<curr_seq<<" "<<runs_in_seq[curr_seq]<<endl;
// 	}



////////////////////////////////////////////////////////////////
// cout<<"Optimize the number of runs per block for each sequence"<<endl;
////////////////////////////////////////////////////////////////

int curr_remainder;
int min_remainder;
int best_j;
int boolean_done;
std::vector<int> runs_per_block(last_seq);

for (int i=first_seq; i<last_seq+1; i++){
	boolean_done = 0;
	for (int j=runs_in_seq[i]/50; j>0; j--){
		if (runs_in_seq[i]/j>100){
			boolean_done = 1;
			}
		if (boolean_done == 0){
			curr_remainder = runs_in_seq[i]%j;
			if (curr_remainder == 0){
				boolean_done = 1;
				}
			if (min_remainder > curr_remainder || j == runs_in_seq[i]/50){
				min_remainder = curr_remainder;
				best_j = j;
				}	
			}
		}
	runs_per_block[i] = runs_in_seq[i]/best_j;	
  	
  	cout<<"Seq. "<<i<<"\t\t"<<runs_in_seq[i]<<" runs\t"<<runs_in_seq[i]/runs_per_block[i]<<" blocks\t"<<runs_per_block[i]<<" runs/block\t"<<runs_in_seq[i]%runs_per_block[i]<<" leftover runs"<<endl;
	}		

////////////////////////////////////////////////////////////////
// cout<<"Determine the first and last run for each block"<<endl;
////////////////////////////////////////////////////////////////

int dummy_run;
int run_counter = 1;
int boolean_newblock = 0;
std::vector<int> first_run_in_block;
std::vector<int> last_run_in_block;

//.push_back 
//.size()

curr_run = first_run;
curr_seq = first_seq;
first_run_in_block.push_back(first_run);

// cout<<"nentries: "<<nentries<<endl;

for (int event = 0; event<nentries; event++){
	aStep.GetEntry(event);
	
	if (boolean_newblock == 1 && curr_seq == aStep.SequenceNum && curr_run != aStep.RunNum){
		curr_run = aStep.RunNum;
		run_counter++;
		first_run_in_block.push_back(curr_run);
		boolean_newblock = 0;			
		}		
	if (curr_seq == aStep.SequenceNum && curr_run != aStep.RunNum){
		run_counter++;		
		curr_run = aStep.RunNum;
		if (run_counter%runs_per_block[curr_seq] == 0){
			last_run_in_block.push_back(curr_run);
			boolean_newblock = 1;
			}
		}
	if (boolean_newblock == 0 && curr_seq != aStep.SequenceNum){
		curr_seq = aStep.SequenceNum;
		curr_run = aStep.RunNum;
		first_run_in_block[first_run_in_block.size()-1] = curr_run;
		run_counter = 1;	
		}
	if (boolean_newblock == 1 && curr_seq != aStep.SequenceNum){
		curr_seq = aStep.SequenceNum;
		curr_run = aStep.RunNum;
		first_run_in_block.push_back(curr_run);
		run_counter = 1;	
		boolean_newblock = 0;
		}		
		
	}	
	
cout<<endl;	
cout<<last_run_in_block.size()<<" blocks total"<<endl;
cout<<endl;	
cout<<"Block\tFirst run\tLast run"<<endl;

for (int i = 0; i<last_run_in_block.size(); i++){
 	cout<<i<<"\t"<<first_run_in_block[i]<<"\t"<<last_run_in_block[i]<<endl;
	}

for (int event = 0; event<nentries; event++){	
 	aStep.GetEntry(event);
	for (int i = 0; i<last_run_in_block.size(); i++){		
		if (aStep.RunNum >= first_run_in_block[i] && aStep.RunNum <= last_run_in_block[i]){
//  	 		cout<<aStep.SequenceNum<<"\t"<<aStep.RunNum<<"\t"<<i<<endl;
			}
		}
	}
}