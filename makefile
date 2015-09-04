CXX=clang++
CFLAGS=-std=c++11
TARGET=ishtypes

all: $(TARGET)

$(TARGET) : $(TARGET).cpp
	$(CXX) $(CFLAGS) -o $(TARGET) $(TARGET).cpp
