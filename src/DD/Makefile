.PHONY: all test draw deps
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

req:
	@echo "--- Updating package lists ---"
	sudo apt update
	@echo "--- Installing dependencies ---"
	sudo apt install -y make build-essential libeigen3-dev libgsl-dev
	@echo "--- Creating symbolic link for Eigen ---"
	@if [ ! -L /usr/include/Eigen ]; then \
		sudo ln -s /usr/include/eigen3/Eigen /usr/include/Eigen; \
		echo "Symbolic link for Eigen created successfully."; \
	else \
		echo "Symbolic link for Eigen already exists. Skipping."; \
	fi
	@echo "--- Development environment setup complete! ---"