#include<iostream>
#include<cstdlib>
#include<sstream>
#include<string>
#include<vector>

using namespace std;

class Vertex{
	static int graphSize;
	static Vertex **idIndex; //pointer array to all vertices

	int id;
	vector<Vertex*> neighbors;
	//for DFS
	bool searched;
public:
	static void init(int); //initialize the static variables
	static int getNumOfVertices();
	static Vertex* getVertexWithId(int);

	Vertex(int);
	int getId();
	void addNeighbor(Vertex*);
	int getDegree();
	Vertex* getNeighbor(int);
	//for DFS
	void setSearched(bool);
	bool isSearched();
};

int Vertex::graphSize;
Vertex **Vertex::idIndex;

void Vertex::init(int graphSize){
	Vertex::graphSize = graphSize;
	idIndex = new Vertex*[graphSize];
	for(int i = 0; i < graphSize; i++){
		idIndex[i] = new Vertex(i);
	}
}

int Vertex::getNumOfVertices(){
	return graphSize;
}

Vertex* Vertex::getVertexWithId(int id){
	return idIndex[id];
}

Vertex::Vertex(int id){
	this->id = id;
	searched = false;
}

int Vertex::getId(){
	return id;
}

void Vertex::addNeighbor(Vertex *neighbor){
	neighbors.push_back(neighbor);
}

int Vertex::getDegree(){
	return (int)neighbors.size();
}

Vertex* Vertex::getNeighbor(int index){
	return neighbors.at(index);
}

void Vertex::setSearched(bool value){
	searched = value;
}

bool Vertex::isSearched(){
	return searched;
}

//main program start from here

void init();
vector<string> split(string, char);

int main(int argc, char *argv[]){
	init(); //read the input graph and build up the structure

	int componentCount = 0;
	//DFS without function recursion
	vector<Vertex*> searchBuffer;
	for(int i = 0; i < Vertex::getNumOfVertices(); i++){
		Vertex *startVertex = Vertex::getVertexWithId(i);
		if(startVertex->isSearched() == false){
			componentCount++;
			startVertex->setSearched(true);
			searchBuffer.push_back(startVertex);
		}
		while(searchBuffer.empty() == false){
			Vertex *vertex = searchBuffer.back();
			searchBuffer.pop_back();
			for(int j = 0; j < vertex->getDegree(); j++){
				Vertex *neighbor = vertex->getNeighbor(j);
				if(neighbor->isSearched() == false){
					neighbor->setSearched(true);
					searchBuffer.push_back(neighbor);
				}
			}
		}
	}
	cout << componentCount << endl;

	return 0;
}

void init(){
	string line;
	getline(cin, line);
	Vertex::init(atoi(line.c_str()));

	while(getline(cin, line)){
		if(line.find("->") != string::npos){
			size_t found = line.find("->");
			Vertex *current = Vertex::getVertexWithId(atoi(line.substr(0, found).c_str()));
			line = line.substr(found + 2);
			vector<string> neighbors = split(line, ',');
			for(int i = 0; i < (int)neighbors.size(); i++){
				int id = atoi(neighbors.at(i).c_str());
				current->addNeighbor(Vertex::getVertexWithId(id));
			}
		}
	}
}

vector<string> split(string str, char delim){
	vector<string> items;
	string item;
	size_t found;
	while((found = str.find(delim)) != string::npos){
		items.push_back(str.substr(0, found));
		str = str.substr(found + 1);
		found = str.find(delim);
	}
	items.push_back(str);
	return items;
}
