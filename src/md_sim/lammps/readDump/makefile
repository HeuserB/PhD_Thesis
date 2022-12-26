# the compiler: gcc for C program, define as g++ for C++
CPP = g++

# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -g -Wall
LIBS = -I /usr/local/lib/read-dump
DUMPFILE = /home/benjamin/Documents/strain_expansion/restart_test.dump
# DUMPFILE = /media/benjamin/BenDrive/Nextcloud/Work/W_PhD/W_PhD_Analysis/Simulations/LAMMPS/strain_expansion/edip/grid_3_45/2000K_8x8x8_50ns.dump

# The build target
TARGET = dumpRead_test

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CPP) $(CFLAGS) -o $(TARGET) $(TARGET).cpp $(LIBS)
test:

	$(CPP) $(CFLAGS) -fopenmp -o $(TARGET) $(TARGET).cpp $(LIBS) && ./$(TARGET) $(DUMPFILE)

clean:
	$(RM) $(TARGET)
