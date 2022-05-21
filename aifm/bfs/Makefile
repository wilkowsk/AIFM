#AIFM runtime configuration
#AIFM_PATH=..
#SHENANGO_PATH=$(AIFM_PATH)/../shenango
#include $(SHENANGO_PATH)/shared.mk
#
#librt_libs = $(SHENANGO_PATH)/bindings/cc/librt++.a
#libaifm = $(AIFM_PATH)/libaifm.a
#INC += -I$(SHENANGO_PATH)/bindings/cc -I$(SHENANGO_PATH)/ksched -I$(AIFM_PATH)/inc -I$(AIFM_PATH)/DataFrame/AIFM/include/
#LDFLAGS += -lmlx5 -lpthread
#
#CXXFLAGS := $(filter-out -std=gnu++17,$(CXXFLAGS))
#override CXXFLAGS += -std=gnu++2a -fconcepts -Wno-unused-function -mcmodel=medium
#
#defule:
#	g++ $(INC) $(LDFLAGS) -fopenmp -msse2 -msse3 -lm -g -O2 -o bfs_up_down_sequential main.cpp bfs_up_down_sequential.cpp graph.cpp 
#test:
#	./bfs_up_down_sequential graphs/tiny.graph
#clean:
#	rm -rf bfs_up_down_sequential *.o *~ *.*~
#
AIFM_PATH=..
SHENANGO_PATH=$(AIFM_PATH)/../shenango
include $(SHENANGO_PATH)/shared.mk

librt_libs = $(SHENANGO_PATH)/bindings/cc/librt++.a
INC += -I$(SHENANGO_PATH)/bindings/cc -I$(SHENANGO_PATH)/ksched -I../inc \
       -I../DataFrame/AIFM/include/

test_bfs_up_down_sequential_src = bfs_up_down_sequential.cpp
test_bfs_up_down_sequential_obj = $(test_bfs_up_down_sequential_src:.cpp=.o)

test_bfs_src = bfs.cpp
test_bfs_obj = $(test_bfs_src:.cpp=.o)

test_graph_src = graph.cpp
test_graph_obj = $(test_graph_src:.cpp=.o)

test_main_src = main.cpp
test_main_obj = $(test_main_src:.cpp=.o)

lib_src = $(wildcard src/*.cpp)
lib_src := $(filter-out src/tcp_device_server.cpp,$(lib_src))
lib_obj = $(lib_src:.cpp=.o)

test_src = $(test_bfs_up_down_sequential_src) $(test_bfs_src) $(test_graph_src) $(test_main_src)
test_obj = $(test_src:.cpp=.o)

src = $(lib_src) $(test_src)
obj = $(src:.cpp=.o)
dep = $(obj:.o=.d)

tcp_device_server_src = ../src/tcp_device_server.cpp
tcp_device_server_obj = $(tcp_device_server_src:.cpp=.o)

override CXXFLAGS += -std=gnu++2a -fconcepts -Wno-unused-function
CXXFLAGS := $(filter-out -std=gnu++17,$(CXXFLAGS))
override LDFLAGS += -lnuma

all: bfs_up_down_sequential graph main libaifm.a

graph: $(test_graph_obj) $(test_main_obj) $(test_bfs_up_down_sequential_obj) $(librt_libs) $(RUNTIME_DEPS) $(lib_obj)
	$(LDXX) -o $@ $^ $(lib_obj) $(librt_libs) $(RUNTIME_LIBS) $(LDFLAGS)

main: $(test_graph_obj) $(test_main_obj) $(test_bfs_up_down_sequential_obj) $(librt_libs) $(RUNTIME_DEPS) $(lib_obj)
	$(LDXX) -o $@ $^ $(lib_obj) $(librt_libs) $(RUNTIME_LIBS) $(LDFLAGS)

bfs_up_down_sequential: $(test_graph_obj) $(test_main_obj) $(test_bfs_up_down_sequential_obj) $(librt_libs) $(RUNTIME_DEPS) $(lib_obj)
	$(LDXX) -o $@ $^ $(lib_obj) $(librt_libs) $(RUNTIME_LIBS) $(LDFLAGS)

$(tcp_device_server_obj): $(tcp_device_server_src)
	$(CXX) $(CXXFLAGS) -c $< -o $@

libaifm.a: $(lib_obj)
	$(AR) rcs $@ $^

#rule to generate a dep file by using the C preprocessor
#(see man cpp for details on the - MM and - MT options)
%.d: %.cpp
	@$(CXX) $(CXXFLAGS) $< -MM -MT $(@:.d=.o) >$@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

ifneq ($(MAKECMDGOALS),clean)
-include $(dep)   # include all dep files in the makefile
endif

.PHONY: clean
clean:
	rm -f src/*.o test/*.o $(dep) bin/* lib*.a
