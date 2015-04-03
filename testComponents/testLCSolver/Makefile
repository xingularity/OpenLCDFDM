all:
	g++-4.8 -O3  -c ../../src/LCD_TimeWaveform.cpp -std=c++11 -I ../../include/
	g++-4.8 -O3  -c ../../src/LCD1D_FDM1DSolver.cpp -std=c++11 -I ../../include/
	g++-4.8 -O3  *.o testLCSolver.cpp -std=c++11 -I ../../include/
clean: 
	rm *.o a.out
