CXX=g++

all: run clean

run: rx.exe
	rx.exe
.phony: run

rx.exe: rx.o ../librtaudio.a
	$(CXX) $^ -lksuser -lmfplat -lmfuuid -lwmcodecdspuuid -lwinmm -lole32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32 -lpthread -o $@

rx.o: rx.cpp ../RtAudio.h ../ThreadSafeQueue.hpp ../dsp.hpp
	$(CXX) $< -c -std=c++11 -o $@ -I../

clean:
	-rm *.o *.exe
.phony: clean
