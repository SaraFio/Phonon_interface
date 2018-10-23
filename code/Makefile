CPP=g++-7.2.0
CPPFLAGS=-I./include -g -Wall

LDLIBS=
LDFLAGS=

EXE := main_try

SRCS=BucNamesFile.cpp MaterialStructureFile.cpp PoscarFile.cpp ForceConstFile.cpp FindNN.cpp get_k.cpp find_vec.cpp getHtot.cpp get_dispersion_direct.cpp cut_FC.cpp sym_FC.cpp check_slabs.cpp extract_submatrices.cpp getHtotOMEN.cpp  asr_sym_scale.cpp get_dispersion.cpp main_phonon_interface.cpp

OBJS = $(SRCS:.cpp=.o)


.PHONY: all clean run

all: $(EXE)

%.o: %.cpp
	$(CPP) -c $(CPPFLAGS) $< -o $@

$(EXE): $(OBJS)
	$(CPP) $(LDFLAGS) $^ $(LDLIBS) -o $@

run:
	./$(EXE)

clean:
	rm -f $(EXE) $(OBJS)
