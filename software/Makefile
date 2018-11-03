
all: darwin_ref_guided read_overlaps 

ntcoding.o: ntcoding.h ntcoding.cpp
	g++ -O4  -c ntcoding.cpp -g

fasta.o: fasta.h fasta.cpp 
	g++ -O4  -c fasta.cpp -g

seed_pos_table.o: ntcoding.o seed_pos_table.h seed_pos_table.cpp
	g++ -O4  -c seed_pos_table.cpp -g

align.o: align.h align.cpp 
	g++ -O4  -c align.cpp -g

gact.o: gact.h gact.cpp 
	g++ -O4  -c gact.cpp -g

darwin_ref_guided.o: gact.h darwin_ref_guided.cpp
	g++ -std=c++11 -O4 -Wall -c darwin_ref_guided.cpp -g

read_overlaps.o: gact.h read_overlaps.cpp
	g++ -std=c++11 -O4 -Wall -c read_overlaps.cpp -g

chameleon.o: Chameleon.cpp
	g++ -std=c++11 -O4 -Wall -o chameleon.o -c Chameleon.cpp -g

config_file.o: ConfigFile.cpp
	g++ -std=c++11 -O4 -Wall -o config_file.o -c ConfigFile.cpp -g

darwin_ref_guided: ntcoding.o fasta.o seed_pos_table.o align.o chameleon.o config_file.o gact.o darwin_ref_guided.o 
	g++ -std=c++11 -O4 -Wall -fopenmp fasta.o ntcoding.o seed_pos_table.o align.o chameleon.o config_file.o gact.o darwin_ref_guided.o -o darwin_ref_guided -pthread 

read_overlaps: ntcoding.o fasta.o seed_pos_table.o align.o chameleon.o config_file.o gact.o read_overlaps.o 
	g++ -std=c++11 -O4 -Wall -fopenmp fasta.o ntcoding.o seed_pos_table.o align.o chameleon.o config_file.o gact.o read_overlaps.o -o read_overlaps -pthread 

clean:
	rm -rf *.o *.out darwin_ref_guided read_overlaps 

