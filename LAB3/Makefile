CXX ?= g++
STD = -std=c++23
OMPFLAG = -fopenmp
CXXFLAGS = -O3 -Wall -Wextra -Wpedantic -march=native
TARGET_SEQUENTIAL ?= correlate_matrix_sequential
TARGET_PARALLEL ?= correlate_matrix_parallel
ARGS ?= 1000 1000

ALL_SOURCES = $(wildcard *.cpp)
SOURCES_SEQUENTIAL = $(filter-out main.cpp, $(ALL_SOURCES))
SOURCES_PARALLEL = $(filter-out main_sequential.cpp, $(ALL_SOURCES))
OBJECTS_SEQUENTIAL = $(SOURCES_SEQUENTIAL:.cpp=.o)
OBJECTS_PARALLEL = $(SOURCES_PARALLEL:.cpp=.o)

all: sequential parallel

sequential: $(TARGET_SEQUENTIAL)

$(TARGET_SEQUENTIAL): $(OBJECTS_SEQUENTIAL)
	@echo "Linking $@"
	@$(CXX) $(STD) $(CXXFLAGS) $(OMPFLAG) -o $@ $^

parallel: $(TARGET_PARALLEL)

$(TARGET_PARALLEL): $(OBJECTS_PARALLEL)
	@echo "Linking $@"
	@$(CXX) $(STD) $(CXXFLAGS) $(OMPFLAG) -o $@ $^

%.o: %.cpp
	@echo "Compiling $<"
	@$(CXX) $(STD) $(CXXFLAGS) $(OMPFLAG) -c -o $@ $<

clean:
	@echo "Cleaning..."
	@rm -f *.o $(TARGET_PARALLEL) $(TARGET_SEQUENTIAL)

run: $(TARGET_PARALLEL)
	@./$(TARGET_PARALLEL) $(ARGS)

run-seq: $(TARGET_SEQUENTIAL)
	@./$(TARGET_SEQUENTIAL) $(ARGS)

.PHONY: all sequential parallel clean run run-seq
