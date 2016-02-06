CXX ?= gcc
CFLAGS = -Wall -Wconversion -fPIC -w -fpermissive -g

all: flow_final

flow_final: flow_final.c svmpredict.c svm.o
	$(CXX) $(CFLAGS) flow_final.c svmpredict.c noiseremoval.c svm.o -o flow_final -lm 

clean:
	rm -f *~ flow_final
