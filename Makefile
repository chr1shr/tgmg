#Load the common configuration file
include ../config.mk

#List of the common source files
objs=poisson_fd.o poisson_fem.o conj_grad.o mpcg_fem.o vpoiss_fem.o
src=$(patsubst %.o,%.cc,$(objs))
tgmg_src=tgmg_config.hh tgmg.hh tgmg.cc
execs=pfd_test pfd_time pfem_test pfem_sym vpfem_test mpcg_test

#Makefile rules
all: $(execs)

lib: $(tgmg_src)

depend: $(src)
	$(cxx) $(cflags) -MM $(src) >Makefile.dep

-include Makefile.dep

pfd_test: pfd_test.cc poisson_fd.o
	$(cxx) $(cflags) -o $@ $^

pfd_time: pfd_time.cc poisson_fd.o
	$(cxx) $(cflags) -o $@ $^

pfem_test: pfem_test.cc poisson_fem.o
	$(cxx) $(cflags) -o $@ $^

vpfem_test: vpfem_test.cc vpoiss_fem.o
	$(cxx) $(cflags) -o $@ $^

pfem_sym: pfem_sym.cc vpoiss_fem.o
	$(cxx) $(cflags) -o $@ $^

mpcg_test: mpcg_test.cc conj_grad.o mpcg_fem.o vpoiss_fem.o
	$(cxx) $(cflags) -o $@ $^

%.o: %.cc
	$(cxx) $(cflags) -c $<

tgmg.cc: tgmg.in.cc tgmg_compile.pl
	perl tgmg_compile.pl

clean:
	rm -f $(execs) $(objs)

.PHONY: all clean lib depend
