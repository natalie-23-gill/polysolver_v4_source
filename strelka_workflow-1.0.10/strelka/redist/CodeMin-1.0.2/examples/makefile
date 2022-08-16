
CM_DIR = ../include

EX1D = codemin_1d_example
EX1D_OBJS = $(EX1D).o

EXCJ = codemin_conj_dir_example
EXCJ_OBJS = $(EXCJ).o

EXCG = codemin_conj_grad_example
EXCG_OBJS = $(EXCG).o

CXXFLAGS = -O3 -Wall -I$(CM_DIR)
#CXXFLAGS += -DCODEMIN_DEBUG
#CXXFLAGS += -DUSE_CBLAS
LDFLAGS = $(CXXFLAGS)

all : $(EX1D) $(EXCJ) $(EXCG)

$(EX1D): $(EX1D_OBJS) 
	$(CXX) $(LDFLAGS) $(EX1D_OBJS) $(LDLIBS) -o $@

$(EXCJ): $(EXCJ_OBJS) 
	$(CXX) $(LDFLAGS) $(EXCJ_OBJS) $(LDLIBS) -o $@

$(EXCG): $(EXCG_OBJS) 
	$(CXX) $(LDFLAGS) $(EXCG_OBJS) $(LDLIBS) -o $@

clean:
	$(RM) $(EX1D_OBJS) $(EXCJ_OBJS) $(EXCG_OBJS)

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<

