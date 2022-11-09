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

all : get_map ac

get_map: bin/get_map

ac: bin/ac

clean :
	- rm -r bin
	- rm -r obj

obj :
	mkdir obj

bin :
	mkdir bin

# All packages for gpu
obj/Entropy_Matrix.o: src/cpp/util/Entropy_Matrix.cpp src/cpp/util/Entropy_Matrix.h | obj    
	$(CXX) -std=c++11 -c src/cpp/util/Entropy_Matrix.cpp -o obj/Entropy_Matrix.o $(CXXFLAGS)
    
obj/Arg_Parser.o: src/cpp/util/Arg_Parser.cpp src/cpp/util/Arg_Parser.h | obj
	$(CXX) -c src/cpp/util/Arg_Parser.cpp -o obj/Arg_Parser.o $(CXXFLAGS)
    
obj/util.o: src/cpp/util/util.cpp src/cpp/util/util.h | obj
	$(CXX) -c src/cpp/util/util.cpp -o obj/util.o $(CXXFLAGS)
    
obj/Residue_Representation.o: src/cpp/util/Residue_Representation.cpp src/cpp/util/Residue_Representation.h | obj
	$(CXX) --std=c++11 -c src/cpp/util/Residue_Representation.cpp -o obj/Residue_Representation.o $(CXXFLAGS)
    
bin/get_map: src/cpp/get_map.cpp obj/Residue_Representation.o obj/Entropy_Matrix.o obj/Arg_Parser.o obj/util.o | bin
	$(CXX) --std=c++11 -O3 src/cpp/get_map.cpp obj/Residue_Representation.o obj/Entropy_Matrix.o obj/Arg_Parser.o obj/util.o -o bin/get_map $(CXXFLAGS)
    
# All packages for autocorrelator
obj/Arg_Parser_ac.o: src/cpp/ac/util/Arg_Parser.cpp src/cpp/ac/util/Arg_Parser.h | obj
	$(CXX) --std=c++17 -c src/cpp/ac/util/Arg_Parser.cpp -o obj/Arg_Parser_ac.o $(CXXFLAGS)
    
obj/io_binary.o : src/cpp/ac/util/io_binary.cpp src/cpp/ac/util/io_binary.h | obj
	$(CXX) --std=c++17 -c src/cpp/ac/util/io_binary.cpp -o obj/io_binary.o $(CXXFLAGS)
    
obj/io_text.o : src/cpp/ac/util/io_text.cpp src/cpp/ac/util/io_text.h | obj
	$(CXX) --std=c++17 -c src/cpp/ac/util/io_text.cpp -o obj/io_text.o $(CXXFLAGS)
    
obj/io.o : obj/io_binary.o obj/io_text.o | obj
	ld -r obj/io_binary.o obj/io_text.o -o obj/io.o
    
obj/util_ac.o : src/cpp/ac/util/util.cpp src/cpp/ac/util/util.h | obj
	$(CXX) --std=c++17 -c src/cpp/ac/util/util.cpp -o obj/util_ac.o $(CXXFLAGS)
    
bin/ac: src/cpp/ac/ac.cpp obj/io.o obj/util_ac.o obj/Arg_Parser_ac.o | bin
	$(CXX) --std=c++17 src/cpp/ac/ac.cpp obj/io.o obj/util_ac.o  obj/Arg_Parser_ac.o -o bin/ac $(CXXFLAGS)
