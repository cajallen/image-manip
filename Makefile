CFLAGS = -g -std=c++17

main:
	make build/main

test:
	make clean
	g++ -fsanitize=address -std=c++14 src/main.cpp src/pixel.cpp src/image.cpp -o test

build/main: build/main.o build/image.o build/pixel.o
	g++ $(CFLAGS) -o $@ $^

build/main.o: src/main.cpp
	@mkdir -p build
	g++ $(CFLAGS) -c -o $@ $^

build/pixel.o: src/pixel.cpp
	@mkdir -p build
	g++ $(CFLAGS) -c -o $@ $^

build/image.o: src/image.cpp
	@mkdir -p build
	g++ $(CFLAGS) -c -o $@ $^

clean:
	-rm -r build