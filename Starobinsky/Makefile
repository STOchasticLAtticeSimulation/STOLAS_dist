MODEL = Starobinsky
SOURCE = ../source

#CXX := g++-13
CXX := g++
#CXX := icpc

#CXXFLAGS := -std=c++17 -O2 -fopenmp
CXXFLAGS := -std=c++11 -O2 -Xpreprocessor -fopenmp -lomp
#CXXFLAGS := -std=c++11 -O2 -qopenmp

LDFLAGS := -L/opt/homebrew/opt/libomp/lib
#CPPFLAGS := -I/opt/homebrew/opt/libomp/include

$(MODEL): $(MODEL).o $(SOURCE)/STOLAS.o
	$(CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(SOURCE)/STOLAS.o

$(MODEL).o: $(SOURCE)/STOLAS.hpp $(SOURCE)/parameters.hpp
STOLAS.o: $(SOURCE)/STOLAS.hpp $(SOURCE)/parameters.hpp $(SOURCE)/vec_op.hpp

clean:
	$(RM) *.o
	$(RM) $(SOURCE)/*.o
	$(RM) $(MODEL)

# MODEL = Starobinsky
# SOURCE = ../source
# DATA = ./data

# #CXX := g++-13
# CXX := g++
# #CXX := icpc

# #CXXFLAGS := -std=c++11 -O2 -fopenmp
# CXXFLAGS := -std=c++11 -O2 -Xpreprocessor -fopenmp -lomp
# #CXXFLAGS := -std=c++11 -O2 -qopenmp

# LDFLAGS := -L/opt/homebrew/opt/libomp/lib
# CPPFLAGS := -I/opt/homebrew/opt/libomp/include

# $(MODEL): $(MODEL).o $(SOURCE)/STOLAS_Sta.o
# 	$(CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $(MODEL) $(MODEL).o $(SOURCE)/STOLAS_Sta.o

# $(MODEL).o: $(SOURCE)/STOLAS_Sta.hpp
# STOLAS_Sta.o: $(SOURCE)/STOLAS_Sta.hpp $(SOURCE)/RK4.hpp $(SOURCE)/vec_op.hpp

# clean:
# 	$(RM) *.o
# 	$(RM) $(SOURCE)/*.o
# 	$(RM) $(MODEL)
