#include <cstdio>
#include <iostream>
#include <vector>
#include <queue>
#include <climits>
#include <omp.h>

#if defined(_DEBUG)
#define dprintf(fmt, ...) printf("%s():%d "fmt,__func__,__LINE__,##__VA_ARGS__)
#else
#define dprintf(fmt, ...)
#endif

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

	string input[width];
	int i = 0;
	for(i = 0; i<width; i++)
	{
		cin>>input[i];
		maze[i]= new int[width];
	}

#pragma omp parallel for collapse(2) firstprivate(maze, input)
	for(int i = 0;i<width;i++)
	{
		//maze[i] = new int[width];

		// '1': wall, '0': road, '2': snake, '3': food
		for(string::size_type index = 0;index < width;index++)
		{
			maze[i][index] = input[i][index] - '0';
			if(input[i][index] == '3')
#pragma omp critical
				dots.push_back(position(i, index));
			else if(input[i][index] == '2')
#pragma omp critical
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
#pragma omp parallel for collapse(2) 
	for(unsigned i = 0 ;i < width;i++)
		for(unsigned j = 0;j < width;j++)
			posMapToIndex[i][j] = -1;

	int dit= 0;
	position dot;
	unsigned int dotSize = (unsigned int)dots.size();
	//vector<position>::const_iterator dot = dots.begin();
#pragma omp parallel for private(dot)
	for(dit = 0; dit<dotSize;dit++)
	{
		dot=dots.at(dit);
		posMapToIndex[dot.row][dot.col] = dit;
	}

//	for(vector<position>::const_iterator dot = dots.begin();dot != dots.end();dot++)
#pragma omp parallel for private(dot) firstprivate(maze)
	for(dit=0; dit<dotSize; dit++)
	{
		dot=dots.at(dit);
		unsigned tmpMaze[width][width];
#pragma omp parallel for collapse(2) firstprivate(maze)
		for(unsigned i = 0;i < width;i++)
			for(unsigned j = 0;j < width;j++)
				tmpMaze[i][j] = maze[i][j];
		
		queue<pos_distance> dotsQueue;
		tmpMaze[dot.row][dot.col] = 1;
		dotsQueue.push(pos_distance(dot, 0));

		while(!dotsQueue.empty())
		{
			const pos_distance tmpDot = dotsQueue.front();
			unsigned row = tmpDot.pos.row, col = tmpDot.pos.col;
			unsigned distance = tmpDot.distance;

			if(posMapToIndex[row][col] != -1) 
			{
				shortest[dit][posMapToIndex[row][col]] = distance;
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

int minValue(unsigned int* numList, unsigned int size)
{
	int i=0, minNum=UINT_MAX;
	for(i=0; i<size; i++)
	{
		if(numList[i]==0)
			break;
		if(numList[i]<minNum )
		{
			minNum=numList[i];
		}
	}
	return (minNum>0)?minNum:0;
}

unsigned TSP(unsigned **shortest, unsigned dotsNum)
{
	unsigned bestDistance = UINT_MAX;
	bool *visited = new bool [dotsNum];
	visited[0] = true;
	for(unsigned i = 1;i < dotsNum;i++)
		visited[i] = false;
	
	TSP_recursive(shortest, 0, dotsNum, 0, bestDistance, visited, 1);

	/*

	int i= 0, j=0, sum=0;
	for(i = 0; i<dotsNum; i++)
	{
		for(j=0; j<dotsNum;j++)
		{
			printf(" %d ", shortest[i][j]);
			//printf("%d", minValue(shortest[i], dotsNum));
		}
printf("min: %d \n", minValue(shortest[i], dotsNum));
	sum+=minValue(shortest[i], dotsNum);
		printf("\n");
	}
	
	dprintf("sum=%d\n", sum);
	*/

	delete [] visited;

	return bestDistance;
}


void TSP_recursive(unsigned **shortest, unsigned nowDot, unsigned dotsNum, unsigned nowDistance, unsigned &bestDistance, bool *visited, unsigned visitedCount)
{
	if((visitedCount == dotsNum) && (nowDistance < bestDistance))
	{
		bestDistance = nowDistance;
	}

	for(unsigned i = 0;i < dotsNum;i++)
	{
		if(!visited[i] && (nowDistance + shortest[nowDot][i] < bestDistance))
		{
			visited[i] = true;
			TSP_recursive(shortest, i, dotsNum, nowDistance + shortest[nowDot][i], bestDistance, visited, visitedCount + 1);
			visited[i] = false;
		}
	}
}
