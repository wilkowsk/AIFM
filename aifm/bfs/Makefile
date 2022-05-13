#AIFM runtime configuration
AIFM_PATH=..
SHENANGO_PATH=$(AIFM_PATH)/../shenango
include $(SHENANGO_PATH)/shared.mk

librt_libs = $(SHENANGO_PATH)/bindings/cc/librt++.a
libaifm = $(AIFM_PATH)/libaifm.a
INC += -I$(SHENANGO_PATH)/bindings/cc -I$(AIFM_PATH)/inc -I$(SHENANGO_PATH)/ksched
LDFLAGS += -lmlx5 -lpthread

CXXFLAGS := $(filter-out -std=gnu++17,$(CXXFLAGS))
override CXXFLAGS += -std=gnu++2a -fconcepts -Wno-unused-function -mcmodel=medium

defule:
	g++ -fopenmp -msse2 -msse3 -lm -g -O2 -o bfs main.cpp bfs.cpp graph.cpp 
test:
	./bfs graphs/tiny.graph
clean:
	rm -rf bfs *.o *~ *.*~
