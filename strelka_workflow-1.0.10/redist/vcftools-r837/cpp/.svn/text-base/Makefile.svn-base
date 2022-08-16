# Make file for vcftools
# Author: Adam Auton
# ($Revision: 230 $)

# Compiler
CPP = g++
# Output executable
EXECUTABLE = vcftools
# Flag used to turn on compilation of PCA routines
ifndef VCFTOOLS_PCA
	VCFTOOLS_PCA = 0
endif
# Compiler flags
CPPFLAGS = -O2 -Wall -Wextra -D_FILE_OFFSET_BITS=64
#CPPFLAGS = -g 
# Included libraries (zlib)
LIB = -lz 
#LIB = -lz -I/opt/local/include/ -L/opt/local/lib/

OBJS = vcftools.o vcf_file.o vcf_entry.o \
		vcf_entry_getters.o vcf_entry_setters.o \
		vcf_file_filters.o vcf_file_output.o \
		vcf_file_format_convert.o \
		vcf_file_diff.o parameters.o \
		vcf_file_index.o \
		output_log.o

ifeq ($(VCFTOOLS_PCA), 1)
	# Define flag for PCA routine compilation
	CPPFLAGS += -DVCFTOOLS_PCA
	# Add LAPACK library
	LIB += -llapack	
	# Add PCA source code
	OBJS+= dgeev.o
endif

vcftools: $(OBJS)
	$(CPP) $(OBJS) -o vcftools $(LIB)
ifdef BINDIR
	cp $(CURDIR)/$@ $(BINDIR)/$@
endif

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

%.o: %.cpp
	$(CPP) -c $(CPPFLAGS) $*.cpp -o $*.o
	$(CPP) -MM $(CPPFLAGS) $*.cpp > $*.d

# remove compilation products
clean:
	@rm -f vcftools *.o *.d
	@rm -f $(BINDIR)/vcftools
