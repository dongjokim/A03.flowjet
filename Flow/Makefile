PROGRAM         = flow_jet

LIBDIRARCH      = lib

CXX             = g++
CXXFLAGS        = -Wall -fPIC

LD              = g++
SOFLAGS         = -shared -Wno-deprecated
CXXFLAGS        += $(shell root-config --cflags)
LIBS            = $(shell root-config --libs)


LIBS	        += -L$(PYTHIA8)/$(LIBDIRARCH) -lpythia8 
INCS            = -I$(PYTHIA8)/include

HDRS		= 
			 
SRCS 		= $(HDRS:.h=.C)
OBJS		= $(HDRS:.h=.o)

%.o: %.C %.h
	$(COMPILE.cc) $(OUTPUT_OPTION) $(INCS) -c $<

#%.o:	$(SRCS) $(HDRS)           
#		$(CXX) $(CXXFLAGS) $(INCS) -c $< -o $@

$(PROGRAM):     $(OBJS) $(HDRS) $(PROGRAM).C 
		$(CXX) $(CXXFLAGS) $(OBJS) $(PROGRAM).C $(LIBS) $(INCS) -o $(PROGRAM)
		@echo "$(PROGRAM) done"

.PHONY : clean debug
clean	:
	@echo cleaning up
	rm -f $(OBJS) core *Dict*  $(PROGRAM).o $(PROGRAM)

debug:
	echo $(OBJS)
