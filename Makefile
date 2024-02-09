CXX = g++
CXXFLAGS = -O3 -Wall

all : get_map denoise

get_map: bin/get_map

denoise: bin/denoise

clean :
	- rm -r bin
	- rm -r obj

obj :
	mkdir obj

bin :
	mkdir bin

# All packages for gpu
obj/Entropy_Matrix.o: artemis/cpp/util/Entropy_Matrix.cpp artemis/cpp/util/Entropy_Matrix.h | obj
	$(CXX) -std=c++11 -c artemis/cpp/util/Entropy_Matrix.cpp -o obj/Entropy_Matrix.o $(CXXFLAGS)
    
obj/Arg_Parser.o: artemis/cpp/util/Arg_Parser.cpp artemis/cpp/util/Arg_Parser.h | obj
	$(CXX) -c artemis/cpp/util/Arg_Parser.cpp -o obj/Arg_Parser.o $(CXXFLAGS)
    
obj/util.o: artemis/cpp/util/util.cpp artemis/cpp/util/util.h | obj
	$(CXX) -c artemis/cpp/util/util.cpp -o obj/util.o $(CXXFLAGS)
    
obj/Residue_Representation.o: artemis/cpp/util/Residue_Representation.cpp artemis/cpp/util/Residue_Representation.h | obj
	$(CXX) --std=c++11 -c artemis/cpp/util/Residue_Representation.cpp -o obj/Residue_Representation.o $(CXXFLAGS)
    
bin/get_map: artemis/cpp/get_map.cpp obj/Residue_Representation.o obj/Entropy_Matrix.o obj/Arg_Parser.o obj/util.o | bin
	$(CXX) --std=c++11 -O3 artemis/cpp/get_map.cpp obj/Residue_Representation.o obj/Entropy_Matrix.o obj/Arg_Parser.o obj/util.o -o bin/get_map $(CXXFLAGS)
    
bin/denoise: artemis/cpp/denoise.cpp obj/Residue_Representation.o obj/Entropy_Matrix.o obj/Arg_Parser.o obj/util.o | bin
	$(CXX) --std=c++11 -O3 artemis/cpp/denoise.cpp obj/Residue_Representation.o obj/Entropy_Matrix.o obj/Arg_Parser.o obj/util.o -o bin/denoise $(CXXFLAGS)
