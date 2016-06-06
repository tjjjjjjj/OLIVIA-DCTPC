
CC=g++
CFLAGS=-Os -g `root-config --cflags`
LFLAGS=`root-config --libs` -L../../ -lMaxCam


all: DeadTime DeadTimeOldData

clean: 
	rm DeadTime
	rm DeadTimeOldData


DeadTime: DeadTime.cc
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)

DeadTimeOldData: DeadTimeOldData.cc
	$(CC) $(CFLAGS) $? -o $@  $(LFLAGS)
