CXX = g++
CXXFLAGS = -O3 -Wall

all : get_map

get_map: bin/get_map

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
    
