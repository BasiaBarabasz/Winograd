CPP=g++
CFLAGS=-Wall -std=c++11
CONV_TC_1D=conv_Toom-Cook_1D
CONV_TC_2D_HUFF=conv_Toom-Cook_2D_Huffman
WIN_MODIF=winograd_modified
WIN_FULL=winograd_modified_full
C_OBJ=bin/utils.o
C_OBJ2=bin/utils_win.o
OBJS=bin/$(CONV_TC_1D).o bin/$(CONV_TC_2D_HUFF).o bin/$(WIN_MODIF).o bin/$(WIN_FULL).o $(C_OBJ) $(C_OBJ2)

help:
	@echo 'Make targets:'
	@echo '  compile        - compile programs'
	@echo '  build          - build programs'
	@echo '  clean          - remove files created by other targets'
	@echo '  test           - run some tests'

bin/%.o: src/%.cpp
	$(CPP) -c -o $@ $< $(CFLAGS)

compile: $(OBJS)
	

build: compile
	$(CPP) -o bin/$(CONV_TC_1D) bin/$(CONV_TC_1D).o $(C_OBJ)
	$(CPP) -o bin/$(CONV_TC_2D_HUFF) bin/$(CONV_TC_2D_HUFF).o $(C_OBJ)
	$(CPP) -o bin/$(WIN_MODIF) bin/$(WIN_MODIF).o $(C_OBJ2)
	$(CPP) -o bin/$(WIN_FULL) bin/$(WIN_FULL).o $(C_OBJ2)

clean:
	rm -f bin/* matrix_output/*

test:
	bin/winograd_modified_full
	bin/winograd_modified 3 4 0 1 -1 2 -1 1 1 1 2 1
	diff matrix_output/ATMat.txt matrix_input/ATMat.txt
	diff BTMat_Win_6_sym.txt matrix_output/BTMat_Win_6_sym.txt

.PHONY: help compile build clean

.DEFAULT_GOAL := build
