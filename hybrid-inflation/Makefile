MODEL = STOLAS
NOISE = noisemap

CXX := g++

CXXFLAGS := -std=c++14 -O3 -Xpreprocessor -fopenmp -march=native -ffast-math -flto -funroll-loops
CPPFLAGS := -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/boost/include

LDFLAGS := -L/opt/homebrew/opt/libomp/lib
LDLIBS := -lfftw3_threads -lfftw3 -lm -lomp

all: $(NOISE) $(MODEL)

$(NOISE): $(NOISE).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(MODEL): $(MODEL).o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)


$(NOISE).o: parameters.hpp src/bias_noise.hpp src/vec_op.hpp src/fft.hpp
$(MODEL).o: src/vec_op.hpp STOLAS.hpp parameters.hpp model.hpp src/zoom.hpp src/output.hpp


clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o
	$(RM) $(NOISE)
	$(RM) $(MODEL)
