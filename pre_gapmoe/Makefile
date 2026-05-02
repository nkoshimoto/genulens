CXX         ?= g++
GSL_CONFIG  ?= gsl-config
GSL_CFLAGS  ?= $(shell $(GSL_CONFIG) --cflags)
GSL_LIBS    ?= $(shell $(GSL_CONFIG) --libs)
GSL_PREFIX  ?= $(shell $(GSL_CONFIG) --prefix)
GSL_RPATH   ?= $(GSL_PREFIX)/lib
CXXFLAGS    ?= -O3 -std=c++11
CPPFLAGS    += -I. $(GSL_CFLAGS)
LDFLAGS     += -Wl,-rpath,$(GSL_RPATH)
LDLIBS      += $(GSL_LIBS)

SHARED_OBJ = galactic_model.o galactic_kinematics.o option.o
TARGETS    = calc_rho_profile calc_mass_dist calc_murel_dist

.PHONY: all clean

all: $(TARGETS)

# ---- Shared objects ----
galactic_model.o: galactic_model.cpp galactic_model.h option.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

option.o: option.cpp option.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# ---- Tool 1: density profile ----
calc_rho_profile: calc_rho_profile.o $(SHARED_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

calc_rho_profile.o: calc_rho_profile.cpp galactic_model.h option.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# ---- Tool 2: mass distribution per component ----
calc_mass_dist: calc_mass_dist.o $(SHARED_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

calc_mass_dist.o: calc_mass_dist.cpp galactic_model.h option.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

galactic_kinematics.o: galactic_kinematics.cpp galactic_model.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# ---- Tool 3: murel distribution ----
calc_murel_dist: calc_murel_dist.o $(SHARED_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

calc_murel_dist.o: calc_murel_dist.cpp galactic_model.h option.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGETS)
