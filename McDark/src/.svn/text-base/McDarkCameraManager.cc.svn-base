#include "McDarkCameraManager.hh"



vector<McDarkCamera*> McDarkCameraManager::listOfCameras;



void 
McDarkCameraManager::addCamera(McDarkCamera *newCam) {

  unsigned int ncam=listOfCameras.size();
  unsigned int icam=0;
  for (; icam<ncam; icam++) {  
    McDarkCamera *cam=listOfCameras[icam];
    if (newCam==cam) break;
  }

  if (icam<ncam) return;

  G4cout << "McDarkCameraManager: adding new camera to digitizer modules!" << G4endl;
  listOfCameras.push_back(newCam);

}

