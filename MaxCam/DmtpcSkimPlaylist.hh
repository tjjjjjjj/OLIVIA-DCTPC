#ifndef DMTPC_SKIM_PLAYLIST_HH
#define DMTPC_SKIM_PLAYLIST_HH

#include "DmtpcSkimDataset.hh"
#include "DmtpcSkimEvent.hh"
#include <vector>

class DmtpcSkimPlaylist : public TObject
{
 
  public:
    DmtpcSkimPlaylist(); 
    DmtpcSkimPlaylist(const char * file); 
    ~DmtpcSkimPlaylist();

    DmtpcSkimDataset * getDataset(){ return cur_dataset; }
    DmtpcSkimEvent * getEvent() { return cur_dataset == NULL ? NULL : cur_dataset->event();} 
    int getCam() {return cams[i];} 
    int getTrack() {return tracks[i];} 

    void next(){go(i+1);} 
    void previous(){go(i-1);}
    int index(){ return i;}
    void go(int index); 
    int n(){ return run_nums.size();}


    void add(const char * det_tag, int run, int ev, int cam, int track); 
    void save(const char * file); 
    void open(const char * file); 

  private:
    
    int i; //!
    string current_det; //!
    int current_run;  //!
    DmtpcSkimDataset * cur_dataset; //! 

    TTree * _tree; //! 
    vector<string> det_tags; //!
    vector<int> run_nums;    //!
    vector<int> events;      //! 
    vector<int> cams;        //! 
    vector<int> tracks;      //! 

    ClassDef(DmtpcSkimPlaylist,0); 

};
#endif
