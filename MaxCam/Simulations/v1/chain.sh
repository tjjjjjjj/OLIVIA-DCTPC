rm output/data/*
rm skim/*
./scripts/runFromFile runParameters/LittleDCTPC_far/fileManager.temp
cp ./output/data/dmtpc_mc_00001.root ./skim
cp ./output/data/dmtpc_mc_00001.root ./
./cleanSkim -c ./seed_tj_mc.cfg ./tmp/fileinfo.txt ./skim/dirinfo.txt 
#from a machine that one can run an X window from:
#./DmtpcSkimViewer /net/hisrv0001/home/spitzj/DCTPC_soft/MaxCam/Simulations/v1/dmtpc_mc_00001.root