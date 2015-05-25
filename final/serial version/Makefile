target = dnn

CPPFLAGS = -O2 -std=c++11 
# -Wall
all: $(target)

%: %.cpp
	g++ $(CPPFLAGS) $< -o $@

dnn: dnn_main.cpp dnn.cpp dnn.h
	g++ $(CPPFLAGS) $^ -o $@

test: test.cpp
	g++ $(CPPFLAGS) $^ -o $@

run_dnn:
	./dnn ./homework4-forests/heart.train
	

# run_test:
# 	./test ./homework4-forests/heart.train