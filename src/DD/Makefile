.PHONY: all test draw
all:
	g++ main.cpp ./../functions/*.cpp -O3 -march=native -ffast-math -DNDEBUG -o main  -lgslcblas -lm -lgsl
	./main
draw:
	g++ draw.cpp ./../functions/*.cpp -O3 -o draw  -lgslcblas -lm -lgsl
	./draw
	python draw.py
test:
	g++ test.cpp ./../functions/*.cpp -O3 -o test  -lgslcblas -lm -lgsl
	./test