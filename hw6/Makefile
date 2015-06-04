# target = libIDAAOSp libIDASOAp mainAOSp.exe mainSOAp.exe
target = mainAOScl.exe mainSOAcl.exe
CFLAGS = -O2 -Wall
CPPFLAGS = -Wall -O3 -std=c++11 -march=native -funroll-loops -fomit-frame-pointer -fopenmp

all: $(target)

# iCBLAS=`pkg-config --cflags cblas`
# lCBLAS=`pkg-config --libs cblas`

iOCL=-I. -I/opt/AMDAPP/current/include -I/opt/AMDAPP/current/include/CL
lOCL=-L/opt/AMDAPP/current/lib/x86_64 -lOpenCL

# %: %.cpp
# 	g++ $(CPPFLAGS) $< -o $@



# hw03: hw02-1.cpp stopWatch.o
# 	g++ $(CPPFLAGS) $^ -o $@

YoUtil.o: YoUtil.cpp
	g++ $(CPPFLAGS) -c $< -o $@

# libIDAAOScl: libIDAAOScl.cpp YoUtil.o
# 	# rm libIDAAOScl.a
# 	# g++ -c -fopenmp libIDAAOScl.cpp
# 	g++ $(CPPFLAGS) $< -o $@ $(iOCL) $(lOCL) -lrt YoUtil.o
# 	ar -q libIDAAOScl.a libIDAAOScl.o
	
# libIDASOAcl: libIDASOAcl.cpp
# 	rm libIDASOAcl.a
# 	g++ -c -fopenmp libIDASOAcl.cpp
# 	ar -q libIDASOAcl.a libIDASOAcl.o

# libIDAFast: libIDAFast.cpp
# 	rm libIDAFast.a
# 	g++ -c -fopenmp libIDAFast.cpp
# 	ar -q libIDAFast.a libIDAFast.o

# libIDAFastcl: libIDAFastcl.cpp
#         rm libIDAFastcl.a
# 	g++ -c -fopenmp libIDAFastcl.cpp
# 	ar -q libIDAFastcl.a libIDAFastcl.o

mainAOScl.exe: main.cpp libIDAAOScl.cpp YoUtil.o
	g++ $(CPPFLAGS) main.cpp libIDAAOScl.cpp -o $@ $(iOCL) $(lOCL) -lrt YoUtil.o

mainSOAcl.exe: main.cpp libIDASOAcl.cpp YoUtil.o
	g++ $(CPPFLAGS) main.cpp libIDASOAcl.cpp -o $@ $(iOCL) $(lOCL) -lrt YoUtil.o

# mainSOAcl.exe: main.cpp libIDASOAcl.cpp stopWatch.o
# 	g++ $(CPPFLAGS) -fopenmp main.cpp stopWatch.o -lIDASOAp -L./ -o $@

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
	./mainAOScl.exe ./test1.dat test1_AOS.out
	./mainSOAcl.exe ./test1.dat test1_SOA.out
	./mainAOScl.exe ./test2.dat test2_AOS.out
	./mainSOAcl.exe ./test2.dat test2_SOA.out

run_generate:
	./generate 2 4 3 10 test_2_dim.dat

clean:
	rm -rf $(target)
	rm -rf *.output
