CXX=clang++

all: run clean

run: tx.exe
	cp tx.exe ..
	../tx.exe
.phony: run

tx.exe: tx.o ../librtaudio.dll.a
	$(CXX) $^ -lksuser -lmfplat -lmfuuid -lwmcodecdspuuid -lwinmm -lole32 -lkernel32 -luser32 -lgdi32 -lwinspool -lshell32 -lole32 -loleaut32 -luuid -lcomdlg32 -ladvapi32 -lpthread -o $@

../librtaudio.dll.a: ../librtaudio.dll

tx.o: tx.cpp ../RtAudio.h ../ThreadSafeQueue.hpp
	$(CXX) $< -c -std=c++11 -o $@ -I../

clean:
	-rm -f *.o *.exe
	-rm -f ../*.exe
.phony: clean
