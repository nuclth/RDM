#SHELL=/bin/sh
 

# NOTE: FOLDERS WILL BE CREATED IF NOT EXISTENT 
OBJ_DIR  := obj_files		# where to put object files
OUTPUT   := output		# create output directory

MAKENAME := makefile		# The name of the makefile
TARGET   := run_rdm		# program command
SRCS     := $(wildcard *.cpp)	# program source files
HDRS     := $(wildcard *.h)	# program header files
OBJS     := $(patsubst %.cpp,%.o,$(SRCS)) 	# program object files
#OBJS     := $(SRCS:.cpp=.o) 	# program object files


###########################################################################
# Commands and options for compiling
########################################################################### 

CC        = g++ -std=c++11
CFLAGS    = -g -O0 -fopenmp 
LDFLAGS   = -g
WARNFLAGS = -Wall -Wextra -Wshadow -fno-common -Wno-unused-parameter
MOREFLAGS = -ansi -pedantic -Wpointer-arith -Wcast-qual -Wcast-align \
           -Wwrite-strings -fshort-enums 
LDLIBS    = -lgsl -lgslcblas -lgomp
 
###########################################################################
# Instructions to compile and link -- allow for different dependencies
########################################################################### 

all: $(TARGET)

$(TARGET): $(OBJS) # $(HDRS) $(MAKENAME)# | $(OBJ_DIR) $(OUTPUT) 
	$(CC) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LDLIBS)

%.o: %.cpp  
	$(CC) $(CFLAGS) -c $< -o $@ $(WARNFLAGS) 


#$(OBJ_DIR)/%.o: %.cpp $(HDRS) 
#	$(CC) $(CFLAGS) -c $< -o $@ 


##########################################################################
# Additional tasks 
##########################################################################

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(OUTPUT):
	mkdir -p $(OUTPUT)

.PHONY: all clean

clean:
	rm -f *.o
#	rm -f $(OBJ_DIR)/*.o

zip:
	zip -r $(COMMAND).zip $(MAKEFILE) $(SRCS) $(HDRS)

##########################################################################
# End of makefile 
##########################################################################
