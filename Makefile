
DEPENDALL=staticparameters.f90
include Makefile.inc

ifeq ($(NODEVICE),1)
ifeq ($(MPARALLEL),1)
 HDIR=matrix_parallel
 DIRS=matrix_parallel
 DLIB=
 HLIB=matrix_parallel/libmparallel.a
 EXE=bfss_mparallel
else
 HDIR=hostcode
 DIRS=hostcode
 DLIB=
 HLIB=hostcode/libhost.a
 EXE=bfss_host_serial
endif
else
 HDIR=hostcode
 DIRS=devicecode hostcode
 DLIB=devicecode/libdevice.a
 HLIB=hostcode/libhost.a
 EXE=testprogram test_timing bfss_host_serial bfss_device_serial
endif

all: $(EXE)

${DLIB}: force_look
	echo looking into devicecode : $(MAKE) $(MFLAGS)
	cd devicecode; $(MAKE) $(MFLAGS)

${HLIB}: force_look
	echo looking into $(HDIR) : $(MAKE) $(MFLAGS)
	cd $(HDIR); $(MAKE) $(MFLAGS)

compare_host_device.mod: compare_host_device.f90
	${FC} -c ${FCFLAGS} $<

bfss_host_serial: ${HLIB} main_host_serial.o 
	${FC} ${FCFLAGS} main_host_serial.o ${HLIB} -o $@

bfss_device_serial: ${HLIB} ${DLIB} main_device_serial.o 
	${FC} ${FCFLAGS} main_device_serial.o ${DLIB} ${HLIB}  -o $@

testprogram:  ${HLIB} ${DLIB} compare_host_device.mod main_testprogram.o compare_host_device.o
	${FC} ${FCFLAGS}  main_testprogram.o compare_host_device.o ${HLIB} ${DLIB}  -o $@

test_timing: main_time.o ${HLIB} ${DLIB}  
	${FC} ${FCFLAGS} main_time.o  ${HLIB} ${DLIB} -o $@

bfss_mparallel: ${HLIB} main_parallel.o
	${FC} ${FCFLAGS} main_parallel.o ${DLIB} ${HLIB}  -o $@

# main_parallel.f90 does not include staticparameters.f90 (DEPENDALL) but these:
ifeq ($(MPARALLEL),1)
main_parallel.o: matrix_parallel/size_parallel.h matrix_parallel/include_parallel.h \
	Fourier.inc unit_number.inc matrix_parallel/remez_md.dat matrix_parallel/remez_pf.dat \
	matrix_parallel/mt19937.f90
endif

.PHONY: clean

clean:
	echo cleaning up in .
	$(RM) -f $(EXE) *.o *.mod
	for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true


