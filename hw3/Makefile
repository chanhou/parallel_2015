target = libIDAAOS libIDASOA libIDAFast mainAOS.exe mainSOA.exe
CFLAGS = -O2 -Wall
CPPFLAGS = -O2 -std=c++11 -Wall

all: $(target)

# iCBLAS=`pkg-config --cflags cblas`
# lCBLAS=`pkg-config --libs cblas`

%: %.cpp
	g++ $(CPPFLAGS) $< -o $@

# hw03: hw02-1.cpp stopWatch.o
# 	g++ $(CPPFLAGS) $^ -o $@

libIDAAOS: libIDAAOS.cpp 
	rm libIDAAOS.a
	g++ -c libIDAAOS.cpp
	ar -q libIDAAOS.a libIDAAOS.o
	
libIDASOA: libIDASOA.cpp
	rm libIDASOA.a
	g++ -c libIDASOA.cpp
	ar -q libIDASOA.a libIDASOA.o

libIDAFast: libIDAFast.cpp
	g++ -c libIDAFast.cpp
	ar -q libIDAFast.a libIDAFast.o

mainAOS.exe: main.cpp libIDAAOS.cpp stopWatch.o
	g++ $(CPPFLAGS) main.cpp stopWatch.o -lIDAAOS -L./ -o $@

mainSOA.exe: main.cpp libIDASOA.cpp stopWatch.o
	g++ $(CPPFLAGS) main.cpp stopWatch.o -lIDASOA -L./ -o $@

generate: generate.cpp
	g++ $(CPPFLAGS) $^ -o $@	

# hw02-1_version1: hw02-1_version1.cpp stopWatch.o
# 	g++ $(CPPFLAGS) $^ -o $@

# hw02-1_version_blas: hw02-1_version2.cpp stopWatch.o
# 	g++ $(CPPFLAGS) $^ -o $@ $(iCBLAS) $(lCBLAS)

run: 
	# ./mainAOS.exe ./test1.dat test1_AOS.out
	# ./mainSOA.exe ./test1.dat test1_SOA.out
	./mainAOS.exe ./test2.dat test2_AOS.out
	./mainSOA.exe ./test2.dat test2_SOA.out

run_generate:
	./generate 2 4 3 10 test_2_dim.dat

clean:
	rm -rf $(target)
	rm -rf *.output
