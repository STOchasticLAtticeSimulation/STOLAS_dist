# Build STOLAS (main simulation) and functions (post-processing).
#
#   make              build both
#   make STOLAS       build only the simulation
#   make functions    build only the post-processing tool
#   make clean
#   make test         regression test of the noise generation against the original implementation
#
# Choose the model from the command line (overrides the default in model.hpp):
#   make MODEL=2              0: chaotic, 1: Starobinsky, 2: USR, 3: hybrid
#   make MODEL=3 NFIELDS=2    number of waterfall fields (hybrid only)
#
# Requirements: a C++20 compiler with OpenMP, FFTW3 (with threads), Boost (header only).
# Non-standard install locations can be given on the command line, e.g.
#   make FFTW_ROOT=/path/to/fftw BOOST_ROOT=/path/to/boost
# The compiler can be changed with CXX=..., e.g.  make CXX=g++-14

SIM   = STOLAS
FUNC  = functions

ifeq ($(origin CXX),default)
CXX := g++
endif

CXXFLAGS := -std=c++20 -O3 -march=native -ffast-math -flto -funroll-loops
CPPFLAGS :=
LDFLAGS  :=
LDLIBS   := -lfftw3_threads -lfftw3 -lm

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  # macOS + Homebrew (brew install libomp fftw boost)
  BREW_PREFIX := $(shell brew --prefix 2>/dev/null || echo /opt/homebrew)
  CXXFLAGS += -Xpreprocessor -fopenmp
  CPPFLAGS += -I$(BREW_PREFIX)/opt/libomp/include -I$(BREW_PREFIX)/include
  LDFLAGS  += -L$(BREW_PREFIX)/opt/libomp/lib -L$(BREW_PREFIX)/lib
  LDLIBS   += -lomp
else
  CXXFLAGS += -fopenmp
endif

# model selection -> -D flags for the compiler
MODELDEFS :=
ifdef MODEL
  ifeq ($(filter 0 1 2 3,$(MODEL)),)
    $(error MODEL must be 0, 1, 2 or 3 (got '$(MODEL)'))
  endif
  MODELDEFS += -DMODEL=$(MODEL)
endif
ifdef NFIELDS
  MODELDEFS += -DNFIELDS=$(NFIELDS)
endif
# options of the regression tests (tests/)
#   REFNOISE=1                      build with the original noise generation
#   FFTW_PLANNER=FFTW_ESTIMATE      reproducible FFTW plans (default: FFTW_MEASURE)
ifdef REFNOISE
  MODELDEFS += -DSTOLAS_REF_NOISE
endif
ifdef FFTW_PLANNER
  MODELDEFS += -DSTOLAS_FFTW_PLANNER=$(FFTW_PLANNER)
endif
CPPFLAGS += $(MODELDEFS)

# Rebuild everything when the selected model changes (the stamp file is
# rewritten only when MODELDEFS differs from the previous build).
STAMP := .build_flags
$(shell printf '%s\n' '$(MODELDEFS)' | cmp -s - $(STAMP) || printf '%s\n' '$(MODELDEFS)' > $(STAMP))

ifdef FFTW_ROOT
  CPPFLAGS += -I$(FFTW_ROOT)/include
  LDFLAGS  += -L$(FFTW_ROOT)/lib -Wl,-rpath,$(FFTW_ROOT)/lib
endif
ifdef BOOST_ROOT
  CPPFLAGS += -I$(BOOST_ROOT)/include
endif

HEADERS := $(wildcard *.hpp src/*.hpp tests/*.hpp)

all: $(SIM) $(FUNC)

$(SIM): $(SIM).o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(FUNC): $(FUNC).o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp $(HEADERS) $(STAMP)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# regression test of the noise generation (see tests/)
tests/test_noise: tests/test_noise.cpp tests/noise_bias_ref.hpp $(HEADERS) $(STAMP)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

test: tests/test_noise
	SIGMA=$${SIGMA:-100.} ./tests/test_noise

clean:
	$(RM) *.o $(SIM) $(FUNC) $(STAMP) tests/test_noise
	$(RM) -r tests/work

.PHONY: all clean test
