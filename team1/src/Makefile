CXX=g++
CXXFLAGS=-O3 -Iteam01/srilm_include

all: disambig_my.cpp
	g++ -Iteam01/srilm_include -o disambig_my disambig_my.cpp -Lteam01 -loolm -ldstruct -lmisc

run:
	./disambig_my team01/ZhuYin-Big5.map team01/ZhuYin team01/bigram.lm > ans.txt 
clean:
	rm -f disambig_my ans.txt
