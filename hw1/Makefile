target = hw1
CFLAGS = -O2 -Wall
CPPFLAGS = -O2 -std=c++11 -Wall

# all: $(target)

all: clean hw1

%: %.cpp
	g++ $(CPPFLAGS) $< -o $@

hw1: hw1.cpp
	g++ $(CPPFLAGS) $< -o $@

clean:
	rm -rf $(target)

run:
	./hw1
