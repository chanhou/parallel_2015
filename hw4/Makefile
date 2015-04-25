target = libIDAAOSp libIDASOAp mainAOSp.exe mainSOAp.exe
CFLAGS = -O2 -Wall
CPPFLAGS = -O2 -std=c++11 -Wall

all: $(target)

# iCBLAS=`pkg-config --cflags cblas`
# lCBLAS=`pkg-config --libs cblas`

%: %.cpp
	g++ $(CPPFLAGS) $< -o $@

# hw03: hw02-1.cpp stopWatch.o
# 	g++ $(CPPFLAGS) $^ -o $@

libIDAAOSp: libIDAAOSp.cpp 
	rm libIDAAOSp.a
	g++ -c -fopenmp libIDAAOSp.cpp
	ar -q libIDAAOSp.a libIDAAOSp.o
	
libIDASOAp: libIDASOAp.cpp
	rm libIDASOAp.a
	g++ -c -fopenmp libIDASOAp.cpp
	ar -q libIDASOAp.a libIDASOAp.o

# libIDAFast: libIDAFast.cpp
# 	g++ -c libIDAFast.cpp
# 	ar -q libIDAFast.a libIDAFast.o

mainAOSp.exe: main.cpp libIDAAOSp.cpp stopWatch.o
	g++ $(CPPFLAGS) -fopenmp main.cpp stopWatch.o -lIDAAOSp -L./ -o $@

mainSOAp.exe: main.cpp libIDASOAp.cpp stopWatch.o
	g++ $(CPPFLAGS) -fopenmp main.cpp stopWatch.o -lIDASOAp -L./ -o $@

# mainAOS.exe: main.cpp libIDAAOS.cpp stopWatch.o
# 	g++ $(CPPFLAGS) main.cpp stopWatch.o -lIDAAOS -L./ -o $@

# mainSOA.exe: main.cpp libIDASOAp.cpp stopWatch.o
# 	g++ $(CPPFLAGS) main.cpp stopWatch.o -lIDASOA -L./ -o $@

generate: generate.cpp
	g++ $(CPPFLAGS) $^ -o $@	

# hw02-1_version1: hw02-1_version1.cpp stopWatch.o
# 	g++ $(CPPFLAGS) $^ -o $@

# hw02-1_version_blas: hw02-1_version2.cpp stopWatch.o
# 	g++ $(CPPFLAGS) $^ -o $@ $(iCBLAS) $(lCBLAS)

run: 
	./mainAOSp.exe ./test1.dat test1_AOS.out
	./mainSOAp.exe ./test1.dat test1_SOA.out
	./mainAOSp.exe ./test2.dat test2_AOS.out
	./mainSOAp.exe ./test2.dat test2_SOA.out

run_generate:
	./generate 2 4 3 10 test_2_dim.dat

clean:
	rm -rf $(target)
	rm -rf *.output
