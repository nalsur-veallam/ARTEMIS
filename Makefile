# Script for GNU Make to build the PARENT_GPU program suite
# Copyright (C) 2020  Markus Fleck

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 3 as 
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

CXX = g++
CXXFLAGS = -O3 -Wall

all : get_map_gpu get_map_cpu

get_map_gpu: bin/get_map_gpu

get_map_cpu: bin/get_map_cpu

clean :
	- rm -r bin
	- rm -r obj

obj :
	mkdir obj

bin :
	mkdir bin

# All packages for gpu
obj/Entropy_Matrix.o: src/cpp/gpu/util/Entropy_Matrix.cpp src/cpp/gpu/util/Entropy_Matrix.h | obj    
	$(CXX) -std=c++11 -c src/cpp/gpu/util/Entropy_Matrix.cpp -o obj/Entropy_Matrix.o $(CXXFLAGS)
    
obj/Arg_Parser.o: src/cpp/gpu/util/Arg_Parser.cpp src/cpp/gpu/util/Arg_Parser.h | obj
	$(CXX) -c src/cpp/gpu/util/Arg_Parser.cpp -o obj/Arg_Parser.o $(CXXFLAGS)
    
obj/util_gpu.o: src/cpp/gpu/util/util.cpp src/cpp/gpu/util/util.h | obj
	$(CXX) -c src/cpp/gpu/util/util.cpp -o obj/util_gpu.o $(CXXFLAGS)
    
obj/Residue_Representation.o: src/cpp/gpu/util/Residue_Representation.cpp src/cpp/gpu/util/Residue_Representation.h | obj
	$(CXX) --std=c++11 -c src/cpp/gpu/util/Residue_Representation.cpp -o obj/Residue_Representation.o $(CXXFLAGS)
    
bin/get_map_gpu: src/cpp/gpu/get_map_gpu.cpp obj/Residue_Representation.o obj/Entropy_Matrix.o obj/Arg_Parser.o obj/util_gpu.o | bin
	$(CXX) --std=c++11 -O3 src/cpp/gpu/get_map_gpu.cpp obj/Residue_Representation.o obj/Entropy_Matrix.o obj/Arg_Parser.o obj/util_gpu.o -o bin/get_map_gpu $(CXXFLAGS)

# All packages for cpu
obj/io_binary.o : src/cpp/cpu/util/io_binary.cpp src/cpp/cpu/util/io_binary.h | obj
	$(CXX) -c src/cpp/cpu/util/io_binary.cpp -o obj/io_binary.o $(CXXFLAGS)
    
obj/io_text.o : src/cpp/cpu/util/io_text.cpp src/cpp/cpu/util/io_text.h | obj
	$(CXX) -c src/cpp/cpu/util/io_text.cpp -o obj/io_text.o $(CXXFLAGS)
    
obj/io.o : obj/io_binary.o obj/io_text.o | obj
	ld -r obj/io_binary.o obj/io_text.o -o obj/io.o
    
obj/util_cpu.o : src/cpp/cpu/util/util.cpp src/cpp/cpu/util/util.h | obj
	$(CXX) -c src/cpp/cpu/util/util.cpp -o obj/util_cpu.o $(CXXFLAGS)
    
bin/get_map_cpu: src/cpp/cpu/get_map_cpu.cpp obj/io.o obj/util_cpu.o obj/Arg_Parser.o | bin
	$(CXX) src/cpp/cpu/get_map_cpu.cpp obj/io.o obj/util_cpu.o  obj/Arg_Parser.o -o bin/get_map_cpu $(CXXFLAGS)
