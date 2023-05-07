CXX=g++

all: run clean

run: tx.exe
	tx.exe
.phony: run

tx.exe: tx.o dsp.o ../librtaudio.a
	$(CXX) $^ -static -lksuser -lmfplat -lmfuuid -lwmcodecdspuuid -lwinmm -lole32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32 -lpthread -o $@

tx.o: tx.cpp ../RtAudio.h ../fifo.hpp ../dsp.hpp
	$(CXX) $< -Wall -c -std=c++11 -o $@ -I../

dsp.o: ../dsp.cpp ../dsp.hpp
	$(CXX) $< -O3 -Wall -c -std=c++11 -mavx2 -mfma -o $@ -I../

clean:
	-del *.o *.exe
.phony: clean
