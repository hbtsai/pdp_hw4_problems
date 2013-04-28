#include <iostream>
#include <vector>
#include <queue>
#include <climits>

using namespace std;

typedef struct position
{
	unsigned row;
	unsigned col;
	position() : row(0), col(0) {}
	position(unsigned row, unsigned col): row(row), col(col) {}
	position(const position& pos): row(pos.row), col(pos.col) {}
	~position() {}
}position;

typedef struct pos_distance
{
	position pos;
	unsigned distance;
	pos_distance(position pos, unsigned distance): pos(pos), distance(distance) {}
	pos_distance(): pos(), distance(0) {}
	pos_distance(const pos_distance& posDis): pos(posDis.pos), distance(posDis.distance) {}
	~pos_distance() {}
}pos_distance;

void findAllShortest(int**, unsigned, const vector<position>&, unsigned**);
unsigned TSP(unsigned**, unsigned);
void TSP_recursive(unsigned**, unsigned, unsigned, unsigned, unsigned&, bool*, unsigned);


int main()
{
	unsigned width;
	cin >> width;
	
	int **maze = new int*[width];
	vector<position> dots;	// store the positions of food and snake

	string input;
	for(int i = 0;cin >> input;i++)
	{
		maze[i] = new int[width];

		// '1': wall, '0': road, '2': snake, '3': food
		for(string::size_type index = 0;index != input.size();index++)
		{
			maze[i][index] = input[index] - '0';
			if(input[index] == '3')
				dots.push_back(position(i, index));
			else if(input[index] == '2')
				dots.insert(dots.begin(), position(i, index)); // insert the position of snake at first element
		}
	}

	unsigned **shortest = new unsigned* [dots.size()];
	for(vector<position>::size_type index = 0;index != dots.size();index++)
		shortest[index] = new unsigned [dots.size()];

	findAllShortest(maze, width, dots, shortest);

	cout << TSP(shortest, dots.size()) << endl;
	
	for(unsigned i = 0;i < width;i++)
		delete [] maze[i];
	delete [] maze;
	for(vector<position>::size_type i = 0;i != dots.size();i++)
		delete [] shortest[i];
	delete [] shortest;
	return 0;
}


void findAllShortest(int **maze, unsigned width, const vector<position> &dots, unsigned **shortest)
{
	int posMapToIndex[width][width];
	for(unsigned i = 0 ;i != width;i++)
		for(unsigned j = 0;j != width;j++)
			posMapToIndex[i][j] = -1;
	for(vector<position>::const_iterator dot = dots.begin();dot != dots.end();dot++)
		posMapToIndex[dot->row][dot->col] = dot - dots.begin();

	for(vector<position>::const_iterator dot = dots.begin();dot != dots.end();dot++)
	{
		unsigned tmpMaze[width][width];
		for(unsigned i = 0;i < width;i++)
			for(unsigned j = 0;j < width;j++)
				tmpMaze[i][j] = maze[i][j];
		
		queue<pos_distance> dotsQueue;
		tmpMaze[dot->row][dot->col] = 1;
		dotsQueue.push(pos_distance(*dot, 0));

		while(!dotsQueue.empty())
		{
			const pos_distance tmpDot = dotsQueue.front();
			unsigned row = tmpDot.pos.row, col = tmpDot.pos.col;
			unsigned distance = tmpDot.distance;

			if(posMapToIndex[row][col] != -1) 
			{
				shortest[dot - dots.begin()][posMapToIndex[row][col]] = distance;
			}

			if(((int)row - 1) >= 0 && tmpMaze[row - 1][col] != 1)
			{
				dotsQueue.push(pos_distance(position(row - 1, col), distance + 1));
				tmpMaze[row - 1][col] = 1;
			}
			if((row + 1) < width && tmpMaze[row + 1][col] != 1)
			{
				dotsQueue.push(pos_distance(position(row + 1, col), distance + 1));
				tmpMaze[row + 1][col] = 1;
			}
			if(((int)col - 1) >= 0 && tmpMaze[row][col - 1] != 1)
			{
				dotsQueue.push(pos_distance(position(row, col - 1), distance + 1));
				tmpMaze[row][col - 1] = 1;
			}
			if((col + 1) < width && tmpMaze[row ][col + 1] != 1)
			{
				dotsQueue.push(pos_distance(position(row, col + 1), distance + 1));
				tmpMaze[row][col + 1] = 1;
			}
		
			dotsQueue.pop();
		}
	}

}

unsigned TSP(unsigned **shortest, unsigned dotsNum)
{
	unsigned bestDistance = UINT_MAX;
	bool *visited = new bool [dotsNum];
	visited[0] = true;
	for(unsigned i = 1;i != dotsNum;i++)
		visited[i] = false;
	
	TSP_recursive(shortest, 0, dotsNum, 0, bestDistance, visited, 1);

	delete [] visited;

	return bestDistance;
}


void TSP_recursive(unsigned **shortest, unsigned nowDot, unsigned dotsNum, unsigned nowDistance, unsigned &bestDistance, bool *visited, unsigned visitedCount)
{
	if((visitedCount == dotsNum) && (nowDistance < bestDistance))
	{
		bestDistance = nowDistance;
	}

	for(unsigned i = 0;i != dotsNum;i++)
	{
		if(!visited[i] && (nowDistance + shortest[nowDot][i] < bestDistance))
		{
			visited[i] = true;
			TSP_recursive(shortest, i, dotsNum, nowDistance + shortest[nowDot][i], bestDistance, visited, visitedCount + 1);
			visited[i] = false;
		}
	}
}
