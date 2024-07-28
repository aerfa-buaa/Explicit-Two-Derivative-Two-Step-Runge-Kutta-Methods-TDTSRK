# 指定编译器
CXX = g++

# 指定编译选项
CXXFLAGS = -std=c++11 -fopenmp

# 指定目标文件名
TARGET = run.out

SRCS=$(wildcard *.cpp)
HEAD=$(wildcard *.hpp)

# print:
#     $(info $(HEAD))

# 指定源文件
# SRCS = main.cpp FluidSimulation.cpp

# 生成目标可执行文件
$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRCS)


# Input.cpp
# Variable_Cons_Dest.cpp
# Driver.cpp
# Burgers.cpp
# Variable.cpp 
# Integration.cpp 
# Initialization.cpp 
# PartDiffeQ.cpp


# Burgers.hpp
# PartDiffEq.hpp 
# Variable.hpp 
# Primitive.hpp

# 生成对象文件
# PartDiffeQ.o: PartDiffeQ.cpp PartDiffeQ.hpp
# 	$(CXX) $(CXXFLAGS) -c PartDiffeQ.cpp

.PHONY: clean
clean:
	rm -f $(TARGET) *.o *.dat
