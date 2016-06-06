#include "SPEFile.hh"

SPEFile::SPEFile (const char* filename) {
  myfile.open(filename, ios::in | ios::binary);
  myfile.seekg(42);
  myfile.read((char*)&myheader, sizeof(myheader));
}

int
SPEFile::Getxpix(){
  return (myheader.xpix);
}

int
SPEFile::Getypix(){
  return (myheader.ypix);
}

int
SPEFile::Getdata_type(){
  return (myheader.data_type);
}

long
SPEFile::Getnimages(){
  return (myheader.nimages);
}

TH2F*
SPEFile::GetImagetoROOT(int a){

  int bytesPerPixel;
  const int constxpix=myheader.xpix;
  const int constypix=myheader.ypix;

  double imagearray[constxpix][constypix];
  

 
 
  if (myheader.data_type == 0){
    // float (4 bytes)
    bytesPerPixel = 4;
int nbytesperim=bytesPerPixel*myheader.xpix*myheader.ypix;
    float tempimage[constxpix][constypix];
    myfile.seekg(4100+a*nbytesperim);
    myfile.read((char*)&tempimage, sizeof(tempimage));
    // this compiles
    //**imagearray = *(double*)tempimage;
    for (int ii=0; ii<myheader.xpix; ii++) {
      for (int jj=0; jj<myheader.xpix; jj++) {
    	imagearray[ii][jj] = (double)tempimage[ii][jj];
	
      }
    }
  }
  else if (myheader.data_type == 1){
   //long (4 bytes)
   bytesPerPixel = 4;
  int nbytesperim=bytesPerPixel*myheader.xpix*myheader.ypix;
   long tempimage[constxpix][constypix];
   myfile.seekg(4100+a*nbytesperim);
   myfile.read((char*)&tempimage, sizeof(tempimage));
   for (int ii=0; ii<myheader.xpix; ii++) {
     for (int jj=0; jj<myheader.xpix; jj++) {
   	imagearray[ii][jj] = (double)tempimage[ii][jj];
     }
   }
  }
  else if (myheader.data_type == 2) {
   //short (2 bytes)
   bytesPerPixel = 2;
  int nbytesperim=bytesPerPixel*myheader.xpix*myheader.ypix;
   short tempimage[constxpix][constypix];
   myfile.seekg(4100+a*nbytesperim);
   myfile.read((char*)&tempimage, sizeof(tempimage));
   for (int ii=0; ii<myheader.xpix; ii++) {
     for (int jj=0; jj<myheader.xpix; jj++) {
   	imagearray[ii][jj] = (double)tempimage[ii][jj];
     }
   }
  }
  else if (myheader.data_type == 3){
    //unsigned short (2 bytes)
    bytesPerPixel = 2;
int nbytesperim=bytesPerPixel*myheader.xpix*myheader.ypix;
    unsigned short tempimage[constxpix][constypix];
    myfile.seekg(4100+a*nbytesperim);
    myfile.read((char*)&tempimage, sizeof(tempimage));
 for (int ii=0; ii<myheader.xpix; ii++) {
     for (int jj=0; jj<myheader.xpix; jj++) {
   	imagearray[ii][jj] = (double)tempimage[ii][jj];
     }
   }
  }

  TH2F *image = new TH2F("image", "image", myheader.xpix, 0, myheader.xpix, myheader.ypix, 0, myheader.ypix);
  for (int ii=0; ii<myheader.xpix; ii++){
    for (int jj=0; jj<myheader.ypix; jj++){
      image->SetBinContent(jj, ii, float(imagearray[myheader.xpix-1-ii][jj]));
      
    }
  }
  
  return image;
}






