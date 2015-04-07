all:
	$(CXX) -g -c ../../src/LCD_TimeWaveform.cpp -std=c++11 -I ../../include/
	$(CXX) -g -c ../../src/LCD1D_FDM1DSolver.cpp -std=c++11 -I ../../include/
	$(CXX) -g *.o testLCSolver.cpp -std=c++11 -I ../../include/
clean:
	rm *.o a.out
