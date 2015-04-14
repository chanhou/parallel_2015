target = hw02-1 hw02-2 hw02-1_version1 hw02-1_version_blas
CFLAGS = -O2 -Wall
CPPFLAGS = -O2 -std=c++11 -Wall

all: $(target)

iCBLAS=`pkg-config --cflags cblas`
lCBLAS=`pkg-config --libs cblas`

%: %.cpp
	g++ $(CPPFLAGS) $< -o $@

hw02-1: hw02-1.cpp stopWatch.o
	g++ $(CPPFLAGS) $^ -o $@

hw02-2: hw02-2.cpp
	g++ $(CPPFLAGS) $^ -o $@

hw02-1_version1: hw02-1_version1.cpp stopWatch.o
	g++ $(CPPFLAGS) $^ -o $@

hw02-1_version_blas: hw02-1_version2.cpp stopWatch.o
	g++ $(CPPFLAGS) $^ -o $@ $(iCBLAS) $(lCBLAS)

clean:
	rm -rf $(target)
	rm -rf *.output


generate: run_hw02-2_vv run_hw02-2_mv run_hw02-2_mm

run: run_hw02-1-v run_hw02-1-mv run_hw02-1-mm

run_version1: run_hw02-1-version1_v run_hw02-1-version1_mv run_hw02-1-version1_mm

run_blas: 
	./hw02-1_version_blas vv-1.dat >> v_version_blas.output
	./hw02-1_version_blas vv-2.dat >> v_version_blas.output
	./hw02-1_version_blas vv-3.dat >> v_version_blas.output
	./hw02-1_version_blas vv-4.dat >> v_version_blas.output
	./hw02-1_version_blas vv-5.dat >> v_version_blas.output
	./hw02-1_version_blas vv-6.dat >> v_version_blas.output
	./hw02-1_version_blas vv-7.dat >> v_version_blas.output
	./hw02-1_version_blas vv-8.dat >> v_version_blas.output
	./hw02-1_version_blas vv-9.dat >> v_version_blas.output
	./hw02-1_version_blas vv-10.dat >> v_version_blas.output
	./hw02-1_version_blas vv-11.dat >> v_version_blas.output
	./hw02-1_version_blas vv-12.dat >> v_version_blas.output
	./hw02-1_version_blas vv-13.dat >> v_version_blas.output
	# ./hw02-1_version_blas vv-14.dat >> v_version_blas.output
	# ./hw02-1_version_blas vv-15.dat >> v_version_blas.output
	# ./hw02-1_version_blas vv-16.dat >> v_version_blas.output
	# ./hw02-1_version_blas vv-17.dat >> v_version_blas.output
	# ./hw02-1_version_blas vv-18.dat >> v_version_blas.output
	# ./hw02-1_version_blas vv-19.dat >> v_version_blas.output
	# ./hw02-1_version_blas vv-20.dat >> v_version_blas.output
	./hw02-1_version_blas mv-1.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-2.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-3.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-4.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-5.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-6.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-7.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-8.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-9.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-10.dat >> mv_version_blas.output
	./hw02-1_version_blas mv-11.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-12.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-13.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-14.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-15.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-16.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-17.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-18.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-19.dat >> mv_version_blas.output
	# ./hw02-1_version_blas mv-20.dat >> mv_version_blas.output
	./hw02-1_version_blas mm-1.dat >> mm_version_blas.output
	./hw02-1_version_blas mm-2.dat >> mm_version_blas.output
	./hw02-1_version_blas mm-3.dat >> mm_version_blas.output
	./hw02-1_version_blas mm-4.dat >> mm_version_blas.output
	./hw02-1_version_blas mm-5.dat >> mm_version_blas.output
	./hw02-1_version_blas mm-6.dat >> mm_version_blas.output
	./hw02-1_version_blas mm-7.dat >> mm_version_blas.output
	./hw02-1_version_blas mm-8.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-9.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-10.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-11.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-12.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-13.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-14.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-15.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-16.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-17.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-18.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-19.dat >> mm_version_blas.output
	# ./hw02-1_version_blas mm-20.dat >> mm_version_blas.output

run_hw02-1-version1_v:
	# rm v_version1.output
	./hw02-1_version1 vv-1.dat >> v_version1.output
	./hw02-1_version1 vv-2.dat >> v_version1.output
	./hw02-1_version1 vv-3.dat >> v_version1.output
	./hw02-1_version1 vv-4.dat >> v_version1.output
	./hw02-1_version1 vv-5.dat >> v_version1.output
	./hw02-1_version1 vv-6.dat >> v_version1.output
	./hw02-1_version1 vv-7.dat >> v_version1.output
	./hw02-1_version1 vv-8.dat >> v_version1.output
	./hw02-1_version1 vv-9.dat >> v_version1.output
	./hw02-1_version1 vv-10.dat >> v_version1.output
	./hw02-1_version1 vv-11.dat >> v_version1.output
	./hw02-1_version1 vv-12.dat >> v_version1.output
	./hw02-1_version1 vv-13.dat >> v_version1.output
	./hw02-1_version1 vv-14.dat >> v_version1.output
	./hw02-1_version1 vv-15.dat >> v_version1.output
	./hw02-1_version1 vv-16.dat >> v_version1.output
	./hw02-1_version1 vv-17.dat >> v_version1.output
	./hw02-1_version1 vv-18.dat >> v_version1.output
	./hw02-1_version1 vv-19.dat >> v_version1.output
	./hw02-1_version1 vv-20.dat >> v_version1.output

run_hw02-1-version1_mv:
	# rm mv.output
	./hw02-1_version1 mv-1.dat >> mv_version1.output
	./hw02-1_version1 mv-2.dat >> mv_version1.output
	./hw02-1_version1 mv-3.dat >> mv_version1.output
	./hw02-1_version1 mv-4.dat >> mv_version1.output
	./hw02-1_version1 mv-5.dat >> mv_version1.output
	./hw02-1_version1 mv-6.dat >> mv_version1.output
	./hw02-1_version1 mv-7.dat >> mv_version1.output
	./hw02-1_version1 mv-8.dat >> mv_version1.output
	./hw02-1_version1 mv-9.dat >> mv_version1.output
	./hw02-1_version1 mv-10.dat >> mv_version1.output
	./hw02-1_version1 mv-11.dat >> mv_version1.output
	./hw02-1_version1 mv-12.dat >> mv_version1.output
	./hw02-1_version1 mv-13.dat >> mv_version1.output
	./hw02-1_version1 mv-14.dat >> mv_version1.output
	./hw02-1_version1 mv-15.dat >> mv_version1.output
	./hw02-1_version1 mv-16.dat >> mv_version1.output
	./hw02-1_version1 mv-17.dat >> mv_version1.output
	./hw02-1_version1 mv-18.dat >> mv_version1.output
	./hw02-1_version1 mv-19.dat >> mv_version1.output
	./hw02-1_version1 mv-20.dat >> mv_version1.output

run_hw02-1-version1_mm:
	# rm mm_version1.output
	./hw02-1_version1 mm-1.dat >> mm_version1.output
	./hw02-1_version1 mm-2.dat >> mm_version1.output
	./hw02-1_version1 mm-3.dat >> mm_version1.output
	./hw02-1_version1 mm-4.dat >> mm_version1.output
	./hw02-1_version1 mm-5.dat >> mm_version1.output
	./hw02-1_version1 mm-6.dat >> mm_version1.output
	./hw02-1_version1 mm-7.dat >> mm_version1.output
	./hw02-1_version1 mm-8.dat >> mm_version1.output
	./hw02-1_version1 mm-9.dat >> mm_version1.output
	./hw02-1_version1 mm-10.dat >> mm_version1.output
	./hw02-1_version1 mm-11.dat >> mm_version1.output
	./hw02-1_version1 mm-12.dat >> mm_version1.output
	./hw02-1_version1 mm-13.dat >> mm_version1.output
	./hw02-1_version1 mm-14.dat >> mm_version1.output
	./hw02-1_version1 mm-15.dat >> mm_version1.output
	./hw02-1_version1 mm-16.dat >> mm_version1.output
	./hw02-1_version1 mm-17.dat >> mm_version1.output
	./hw02-1_version1 mm-18.dat >> mm_version1.output
	./hw02-1_version1 mm-19.dat >> mm_version1.output
	./hw02-1_version1 mm-20.dat >> mm_version1.output

run_hw02-1-v:
	rm v.output
	./hw02-1 vv-1.dat >> v.output
	./hw02-1 vv-2.dat >> v.output
	./hw02-1 vv-3.dat >> v.output
	./hw02-1 vv-4.dat >> v.output
	./hw02-1 vv-5.dat >> v.output
	./hw02-1 vv-6.dat >> v.output
	./hw02-1 vv-7.dat >> v.output
	./hw02-1 vv-8.dat >> v.output
	./hw02-1 vv-9.dat >> v.output
	./hw02-1 vv-10.dat >> v.output
	./hw02-1 vv-11.dat >> v.output
	./hw02-1 vv-12.dat >> v.output
	./hw02-1 vv-13.dat >> v.output
	./hw02-1 vv-14.dat >> v.output
	./hw02-1 vv-15.dat >> v.output
	./hw02-1 vv-16.dat >> v.output
	./hw02-1 vv-17.dat >> v.output
	./hw02-1 vv-18.dat >> v.output
	./hw02-1 vv-19.dat >> v.output
	./hw02-1 vv-20.dat >> v.output
	
run_hw02-1-mv:
	rm mv.output
	./hw02-1 mv-1.dat >> mv.output
	./hw02-1 mv-2.dat >> mv.output
	./hw02-1 mv-3.dat >> mv.output
	./hw02-1 mv-4.dat >> mv.output
	./hw02-1 mv-5.dat >> mv.output
	./hw02-1 mv-6.dat >> mv.output
	./hw02-1 mv-7.dat >> mv.output
	./hw02-1 mv-8.dat >> mv.output
	./hw02-1 mv-9.dat >> mv.output
	./hw02-1 mv-10.dat >> mv.output
	./hw02-1 mv-11.dat >> mv.output
	./hw02-1 mv-12.dat >> mv.output
	./hw02-1 mv-13.dat >> mv.output
	./hw02-1 mv-14.dat >> mv.output
	./hw02-1 mv-15.dat >> mv.output
	./hw02-1 mv-16.dat >> mv.output
	./hw02-1 mv-17.dat >> mv.output
	./hw02-1 mv-18.dat >> mv.output
	./hw02-1 mv-19.dat >> mv.output
	./hw02-1 mv-20.dat >> mv.output


run_hw02-1-mm:
	rm mm.output
	./hw02-1 mm-1.dat >> mm.output
	./hw02-1 mm-2.dat >> mm.output
	./hw02-1 mm-3.dat >> mm.output
	./hw02-1 mm-4.dat >> mm.output
	./hw02-1 mm-5.dat >> mm.output
	./hw02-1 mm-6.dat >> mm.output
	./hw02-1 mm-7.dat >> mm.output
	./hw02-1 mm-8.dat >> mm.output
	./hw02-1 mm-9.dat >> mm.output
	./hw02-1 mm-10.dat >> mm.output
	./hw02-1 mm-11.dat >> mm.output
	./hw02-1 mm-12.dat >> mm.output
	./hw02-1 mm-13.dat >> mm.output
	./hw02-1 mm-14.dat >> mm.output
	./hw02-1 mm-15.dat >> mm.output
	./hw02-1 mm-16.dat >> mm.output
	./hw02-1 mm-17.dat >> mm.output
	./hw02-1 mm-18.dat >> mm.output
	./hw02-1 mm-19.dat >> mm.output
	./hw02-1 mm-20.dat >> mm.output


run_hw02-2_vv:
	rm vv*.dat
	./hw02-2 1 64 64 1 vv-1.dat
	./hw02-2 1 128 128 1 vv-2.dat
	./hw02-2 1 256 256 1 vv-3.dat
	./hw02-2 1 512 512 1 vv-4.dat
	./hw02-2 1 1024 1024 1 vv-5.dat
	./hw02-2 1 2048 2048 1 vv-6.dat
	./hw02-2 1 4096 4096 1 vv-7.dat
	./hw02-2 1 8192 8192 1 vv-8.dat
	./hw02-2 1 16384 16384 1 vv-9.dat
	./hw02-2 1 32768 32768 1 vv-10.dat
	./hw02-2 1 65536 65536 1 vv-11.dat
	./hw02-2 1 131072 131072 1 vv-12.dat
	./hw02-2 1 262144 262144 1 vv-13.dat
	./hw02-2 1 524288 524288 1 vv-14.dat
	./hw02-2 1 1048576 1048576 1 vv-15.dat
	./hw02-2 1 2097152 2097152 1 vv-16.dat
	./hw02-2 1 4194304 4194304 1 vv-17.dat
	./hw02-2 1 8388608 8388608 1 vv-18.dat
	./hw02-2 1 16777216 16777216 1 vv-19.dat
	./hw02-2 1 33554432 33554432 1 vv-20.dat

run_hw02-2_mv:
	rm mv*.dat
	./hw02-2 2 32 32 1 mv-1.dat
	./hw02-2 4 64 64 1 mv-2.dat
	./hw02-2 6 128 128 1 mv-3.dat
	./hw02-2 8 256 256 1 mv-4.dat
	./hw02-2 10 512 512 1 mv-5.dat
	./hw02-2 12 1024 1024 1 mv-6.dat
	./hw02-2 14 2048 2048 1 mv-7.dat
	./hw02-2 16 4096 4096 1 mv-8.dat
	./hw02-2 18 8192 8192 1 mv-9.dat
	./hw02-2 20 16384 16384 1 mv-10.dat
	./hw02-2 22 32768 32768 1 mv-11.dat
	./hw02-2 24 65536 65536 1 mv-12.dat
	./hw02-2 26 131072 131072 1 mv-13.dat
	./hw02-2 28 262144 262144 1 mv-14.dat
	./hw02-2 30 524288 524288 1 mv-15.dat
	./hw02-2 32 1048576 1048576 1 mv-16.dat
	./hw02-2 34 1011111 1011111 1 mv-17.dat
	./hw02-2 36 2345123 2345123 1 mv-18.dat
	./hw02-2 38 2800000 2800000 1 mv-19.dat
	./hw02-2 40 3100000 3100000 1 mv-20.dat


run_hw02-2_mm:
	rm mm*.dat
	./hw02-2 2 32 32 2 mm-1.dat
	./hw02-2 4 64 64 4 mm-2.dat
	./hw02-2 8 128 128 8 mm-3.dat
	./hw02-2 16 256 256 16 mm-4.dat
	./hw02-2 32 700 700 32 mm-5.dat
	./hw02-2 56 901 901 42 mm-6.dat
	./hw02-2 60 2500 2500 60 mm-7.dat
	./hw02-2 66 7000 7000 66 mm-8.dat
	./hw02-2 77 12000 12000 77 mm-9.dat
	./hw02-2 80 14231 14231 82 mm-10.dat
	./hw02-2 85 23124 23124 87 mm-11.dat
	./hw02-2 90 25341 25341 89 mm-12.dat
	./hw02-2 101 45678 45678 101 mm-13.dat
	./hw02-2 123 56789 56789 123 mm-14.dat
	./hw02-2 234 99999 99999 231 mm-15.dat
	./hw02-2 256 100000 100000 256 mm-16.dat
	./hw02-2 268 100231 100231 268 mm-17.dat
	./hw02-2 300 110000 110000 300 mm-18.dat
	./hw02-2 350 120000 120000 350 mm-19.dat
	./hw02-2 400 130000 130000 400 mm-20.dat
	