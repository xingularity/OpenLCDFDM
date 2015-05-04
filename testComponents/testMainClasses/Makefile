all:
	$(CXX) -O3 -fopenmp -c ../../src/LCD_TimeWaveform.cpp -std=c++11 -I ../../include/
	$(CXX) -O3 -fopenmp -c ../../src/LCD1D_FDM1DSolver.cpp -std=c++11 -I ../../include/
	$(CXX) -O3 -fopenmp -c ../../src/LCD_Optics2x2.cpp -std=c++11 -I ../../include/
	$(CXX) -O3 -fopenmp -c ../../src/LCD_SpectrumInterpolator.cpp -std=c++11 -I ../../include/
	$(CXX) -O3 -fopenmp -c ../../src/LCD1D_ExtendedJones.cpp -std=c++11 -I ../../include/
	$(CXX) -O3 -fopenmp -c ../../src/LCD1D_LCD1DMain.cpp -std=c++11 -I ../../include/
	$(CXX) -O3 -fopenmp *.o testMainClasses.cpp -std=c++11 -I ../../include/

clean:
	rm *.o a.out
