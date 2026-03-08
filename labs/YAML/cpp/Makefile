# Compiler
CXX = g++
CXXFLAGS = -O2 -Wall

# Target
TARGET = simulation

# Source
SRC = simulation.cpp

# Libraries
LIBS = -lyaml-cpp

# Default rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LIBS)

# Run example
run: $(TARGET)
	./$(TARGET)

# Clean build artifacts
clean:
	rm -f $(TARGET)
