OPSYS=$(shell uname)

PLATFORM=$(shell uname -p)
ARCH=.$(PLATFORM)

CXX:=ccache $(CXX)

CXXFLAGS+=-std=c++17 -DROOT_DIR=\"$(shell pwd ..)/\"

CPPFLAGS+=$(LOCAL_CFLAGS)
LDFLAGS+=$(LOCAL_LDFLAGS)


# CPPFLAGS+=-I./include 
CPPFLAGS+=-I./include -I /usr/include/eigen3 
LDFLAGS+=-lpthread -lyaml-cpp 

# CXXFLAGS+= -g -march=native
CXXFLAGS+= -O3 -march=native

