.PHONY: all ddg test
all:
	clear
	g++ DDGIC.cpp ./functions/*.cpp -O3 -o DDGIC  -lgslcblas -lm -lgsl
	./DDGIC
ddg:
	g++ DDG.cpp ./functions/*.cpp -O3 -o DDG  -lgslcblas -lm -lgsl
	./DDG
test:
	g++ test.cpp ./functions/*.cpp -O3 -o test  -lgslcblas -lm -lgsl
	./test