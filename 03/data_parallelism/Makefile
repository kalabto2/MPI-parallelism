PROGRAM=kalabto2
CC=g++
CFLAGS=-std=c++20 -Wall -pedantic -O3 -fopenmp # -Wextra -g -Wno-long-long -O0 -ggdb
#DIRNAME=files2/answers files2/answerSheets files2/exports files2/questions files2/quizzes files2/sheets

## Spuštění programu 'all'
.PHONY: all
all: compile

### Vygenerování dokumentace
#.PHONY: doc
#doc: src/main.cpp src/App.h
#	doxygen Doxyfile

## Kompilace
compile: $(PROGRAM)

## Spuštění
run: compile
	./$(PROGRAM)

## Spuštění programu 'main'
main: build/main
	./build/main

## Odstranění vygenerovatelných souborů
.PHONY: clean
clean:
	rm -rf ./build ./$(PROGRAM)

#.PHONY: dir
#dir:
#	mkdir -p $(DIRNAME)

# PROGRAMY
## Program 'kalabto2'
$(PROGRAM): build/main.o build/Instance_handler.o build/Graph.o
	mkdir -p build
	$(CC) $(CFLAGS) $^ -o $@

# Třídy
## Třída 'main.cpp'
build/main.o: src/main.cpp src/Instance_handler.h
	mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

## Třída 'Instance_handler.cpp'
build/Instance_handler.o: src/Instance_handler.cpp src/Instance_handler.h src/Graph.h
	mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@

## Třída 'Graph.cpp'
build/Graph.o: src/Graph.cpp src/Graph.h
	mkdir -p build
	$(CC) $(CFLAGS) -c $< -o $@


