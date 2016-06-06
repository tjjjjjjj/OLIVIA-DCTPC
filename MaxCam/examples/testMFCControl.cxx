// code to test the digital control of the MFC purge mode

{

  MaxCamChannel* mfcpurge = new MaxCamChannel("mfc_purge");
  mfcpurge->setDio(1);
  gSystem->Sleep(3000);
  mfcpurge->setDio(0);
  gSystem->Sleep(3000);
  mfcpurge->setDio(1);
  gSystem->Sleep(3000);
  mfcpurge->setDio(0);

}
