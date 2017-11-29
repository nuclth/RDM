# NOTE: FOLDERS WILL BE CREATED IF NOT EXISTENT 
OBJ_DIR  = obj
SRC_DIR  = src
ME_DIR   = me_files
FLAG_DIR = flag_files
OBME_DIR = obme
TBME_DIR = tbme
REF_DIR  = ref_files

OUTPUT   = sdp_files

BINARY   = run_rdm
SRCS     = $(wildcard $(SRC_DIR)/*.cpp)
OBJS     = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))


###########################################################################
# Commands and options for compiling
########################################################################### 

CC        = g++ -std=c++11
CFLAGS    = -g -O3 -fopenmp 
LDFLAGS   = -g
WARNFLAGS = -Wall -Wextra -Wshadow -fno-common -Wno-unused-parameter
MOREFLAGS = -ansi -pedantic -Wpointer-arith -Wcast-qual -Wcast-align \
           -Wwrite-strings -fshort-enums 
LDLIBS    = -lgsl -lgslcblas -lgomp
 
###########################################################################
# Instructions to compile and link -- allow for different dependencies
########################################################################### 

all: $(BINARY) 

$(BINARY): $(OBJS) | $(OBJ_DIR)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR) 
	$(CC) $(CFLAGS) -c $< -o $@ $(WARNFLAGS) 


#$(OBJ_DIR)/%.o: %.cpp $(HDRS) 
#	$(CC) $(CFLAGS) -c $< -o $@ 


##########################################################################
# Additional tasks 
##########################################################################

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p $(OUTPUT)
	mkdir -p $(FLAG_DIR)
	mkdir -p $(ME_DIR)
	cd $(ME_DIR) && mkdir -p $(OBME_DIR)
	cd $(ME_DIR) && mkdir -p $(TBME_DIR)
	cd $(ME_DIR) && mkdir -p $(REF_DIR)

.PHONY: all clean

clean:
	rm -f $(OBJ_DIR)/*.o

zip:
	zip -r $(COMMAND).zip $(MAKEFILE) $(SRCS) $(HDRS)

##########################################################################
# End of makefile 
##########################################################################
