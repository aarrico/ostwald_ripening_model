CC = g++
SRC = ./src
BUILD = ./build
OPT = -Wwrite-strings

proj2 : proj2.o Crystal.o
	$(CC) -fopenmp $(BUILD)/proj2.o $(BUILD)/Crystal.o -o $(BUILD)/proj2

proj2.o : $(SRC)/proj2.cpp
	$(CC)  -c $(OPT) $< -o  $(BUILD)/proj2.o

#RK4.o : $(SRC)/RK4.h
#	$(CC) -c $< -o $(BUILD)/RK4.o

Crystal.o : $(SRC)/Crystal.cpp $(SRC)/Crystal.h
	$(CC) -c $< -o $(BUILD)/Crystal.o

clean : 
	rm $(BUILD)/*
