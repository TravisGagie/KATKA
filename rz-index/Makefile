DEPS   ?= ../deps
SDSL   ?= $(DEPS)/sdsl
RINDEX ?= $(DEPS)/r-index-pfp/internal
CXX    ?= g++
CXXFLAGS ?= -std=c++17 -O3 -march=native -DNDEBUG -w
INC    = -I$(SDSL)/include -I$(RINDEX)/..
LIB    = -L$(SDSL)/lib -lsdsl -ldivsufsort -ldivsufsort64 -lz

all: rz-prep rz-lz77 rz-build rz-bench rz-genpat rz-prof rz-mvbuild rz-mvtest sr-build rz-aux rz-dgcheck rz-classify rz-csabuild rz-vfybuild rz-tagbuild rz-tagruns
test-tools: mbwt ri-build-from

rz-prep: rz_prep.cpp digest.hpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lz
rz-lz77: rz_lz77.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB) -pthread
rz-build: rz_build.cpp rz_index.hpp digest.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB) -pthread
rz-bench: rz_bench.cpp rz_index.hpp digest.hpp rlcsa.hpp sr.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
mbwt: mbwt.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
ri-build-from: ri_build.cpp
	$(CXX) $(CXXFLAGS) $(INC) -I$(RINDEX) -o $@ $< $(LIB)
rz-genpat: rz_genpat.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f rz-prep rz-lz77 rz-build rz-bench rz-genpat rz-prof rz-mvbuild rz-mvtest sr-build rz-aux rz-dgcheck rz-classify rz-csabuild rz-vfybuild rz-tagbuild rz-tagruns mbwt ri-build-from
rz-prof: rz_prof.cpp rz_index.hpp digest.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
rz-mvbuild: rz_mvbuild.cpp mv.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
rz-mvtest: rz_mvtest.cpp rz_index.hpp digest.hpp mv.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
sr-build: sr_build.cpp sr.hpp rz_index.hpp digest.hpp mv.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
rz-aux: rz_aux.cpp rz_index.hpp digest.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
rz-dgcheck: rz_dgcheck.cpp rz_index.hpp digest.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
rz-classify: rz_classify.cpp rz_index.hpp digest.hpp rlcsa.hpp sr.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
rz-csabuild: rz_csabuild.cpp rz_index.hpp rlcsa.hpp runrank.hpp efrank.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)

rz-vfybuild: rz_vfybuild.cpp vfy.hpp digest.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)

rz-tagbuild: rz_tagbuild.cpp tag.hpp rz_index.hpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)

rz-tagruns: rz_tagruns.cpp
	$(CXX) $(CXXFLAGS) $(INC) -o $@ $< $(LIB)
