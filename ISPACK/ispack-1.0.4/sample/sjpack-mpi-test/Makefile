SHELL = /bin/sh

include ../../Mkinclude

LIBDIR = ../../
LIB = isp

clean:
	rm -f *.out
	rm -f *~

%.o : %.f
	$(MPIFC) $(MPIFFLAGS) -c $<

%-align.out: %-align.o %-align-sub.o
	$(MPIFC) $(MPIFFLAGS) $^ -L$(LIBDIR) -l$(LIB) -o $@

%.out: %.o
	$(MPIFC) $(MPIFFLAGS) $< -L$(LIBDIR) -l$(LIB) -o $@

