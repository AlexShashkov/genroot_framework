# the compiler: gcc for C program, define as g++ for C++
CC = gcc
CXX = g++
# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -std=c++20 
# LINKING =
TARGET = main

all:
	$(CXX) $(CFLAGS) -o out $(TARGET).cpp $(LINKING)
	./out

clean:
	$(RM) $(TARGET)
