.PHONY: all test draw
all:
	g++ -O3 -march=native  -fopenmp -ffast-math -DNDEBUG main.cpp ./../functions/*.cpp  -o  main  -lgslcblas -lm -lgsl
	./main
draw:
	g++ draw.cpp ./../functions/*.cpp -O3 -o draw  -lgslcblas -lm -lgsl
	./draw
	python draw.py
test:
	g++ test.cpp ./../functions/*.cpp -O3 -o test  -lgslcblas -lm -lgsl
	./test