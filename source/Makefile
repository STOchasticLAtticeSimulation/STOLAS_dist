NOISE = noisemap
BIAS = biasmap
PARAM = parameters

#CXX := g++-13
CXX := g++
#CXX := icpc

#CXXFLAGS := -std=c++11 -O2 -fopenmp
CXXFLAGS := -std=c++11 -O2 -Xpreprocessor -fopenmp -lomp
#CXXFLAGS := -std=c++11 -O2 -qopenmp

LDFLAGS := -L/opt/homebrew/opt/libomp/lib
CPPFLAGS := -I/opt/homebrew/opt/libomp/include

all: $(NOISE) $(BIAS)
$(NOISE): $(NOISE).o
	$(CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $(NOISE) $(NOISE).o
$(BIAS): $(BIAS).o
	$(CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $(BIAS) $(BIAS).o

$(NOISE).o: parameters.hpp bias_noise.hpp vec_op.hpp fft.hpp
$(BIAS).o: parameters.hpp bias_noise.hpp vec_op.hpp fft.hpp

clean:
	$(RM) *.o
	$(RM) $(NOISE)
	$(RM) $(BIAS)
