CC = gcc -I /home/chenxi/local/include -L /home/chenxi/local/lib
CXX = g++ -I /home/chenxi/local/include -L /home/chenxi/local/lib

SWA_FLAGS = -DDEB=0 -DRDT=0 -DMAXI=0 -DNEW=1 -DSORT_PAIRS=0
MEM_FLAGS = -DPAIRED_END=1 -DMAINY=0 -DSAIS=1
CPPFLAGS = -DENABLE_PREFETCH $(MEM_FLAGS) $(SWA_FLAGS)
LIBS = -lhts -lbz2 -llzma -lpthread -lm -lz

CFLAGS = -O0 -std=gnu99 -g -march=native -DBGZF_MT
CXXFLAGS = -O0 -std=c++0x -g -march=native -fpermissive -DBGZF_MT

COBJ = utils.o \
       xth.o \
       sam_bio.o \
       bedidx.o \
       sample.o \
       bio.o \
       tmp.o \
       str.o \
       aux.o \
       cigar.o \
       hash.o \
       read.o \
			 simd_check.o \
       hash_func.o

CXXOBJ = main.o \
         ebam.o \
         sam_ext.o \
         mem2_main.o \
				 build_main.o \
				 sort_and_mark.o \
         mem2/bandedSWA.o \
         mem2/bntseq.o \
         mem2/bwa.o \
         mem2/bwamem.o \
         mem2/bwamem_extra.o \
         mem2/bwamem_pair.o \
         mem2/bwtbuild.o \
         mem2/bwtindex.o \
         mem2/fastmap.o \
         mem2/FMI_search.o \
         mem2/kstring.o \
         mem2/profiling.o \
         mem2/ksw.o \
         mem2/kswv.o \
         mem2/kthread.o \
         mem2/read_index_ele.o \
         mem2/utils.o

wgs_pipe: $(COBJ) $(CXXOBJ)
	$(CXX) -o fast_aligner $(CXXFLAGS) $(COBJ) $(CXXOBJ) $(LIBS)

.PHONY: clean
clean:
	rm -f */*.o *.o fast_aligner
