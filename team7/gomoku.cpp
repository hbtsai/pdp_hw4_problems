#include<iostream>


using namespace std;
	


int check_win(int size, int** board, int i, int j)
{
	int who = board[i][j];
	int count1, count2;

	// up + down
	count1 = 0;
	count2 = 0;
	for(int k=1;k<=4;k++){
		if(j+k >= size || board[i][j+k] != who) break;
		count1++;
	}
	for(int k=1;k<=4;k++){
		if(j-k < 0 || board[i][j-k] != who) break;
		count2++;
	}
	if(count1 + count2 >=4) return 1;

	// left + right
	count1 = 0;
	count2 = 0;
	for(int k=1;k<=4;k++){
		if(i+k >= size || board[i+k][j] != who) break;
		count1++;
	}
	for(int k=1;k<=4;k++){
		if(i-k < 0 || board[i-k][j] != who) break;
		count2++;
	}
	if(count1 + count2 >=4) return 1;

	// upleft + downright
	count1 = 0;
	count2 = 0;
	for(int k=1;k<=4;k++){
		if(j+k >= size || i-k < 0 || board[i-k][j+k] != who) break;
		count1++;
	}
	for(int k=1;k<=4;k++){
		if(j-k < 0 || i+k >= size ||  board[i+k][j-k] != who) break;
		count2++;
	}
	if(count1 + count2 >=4) return 1;


	// upright + downleft
	count1 = 0;
	count2 = 0;
	for(int k=1;k<=4;k++){
		if(i+k >= size || j+k>=size || board[i+k][j+k] != who) break;
		count1++;
	}
	for(int k=1;k<=4;k++){
		if(i-k < 0 || j-k<0 || board[i-k][j-k] != who) break;
		count2++;
	}
	if(count1 + count2 >=4) return 1;


	return 0;
}

// NegaMax algorithm
int nega_max(int who, int** board, int size, int depth, int max_depth, int alpha, int beta)
{
	int val = alpha;
	int tmp;
	int status;
	// iterate through every child
	for(int i = 0; i < size; i++){
		for(int j = 0; j < size; j++){
			if( board[i][j] != 0) continue;	
			board[i][j] = who;

			status = check_win( size, board, i, j);
			if( status == 1 || depth == max_depth-1 ) tmp = status;
			else tmp = -nega_max( -who, board, size, depth+1, max_depth, -beta, -val );

			board[i][j] = 0;
			if( tmp >= beta ) return tmp;
			if( tmp > val ) val = tmp;
		}
	}
	return val;

}



int main()
{

	
	int num_given_white;
	int num_given_black;
	int size;				// board_size, N x N
	int depth_to_search;	// maximum search depth
	int x;
	int y;
	int** board;			// white = 1, black = -1, empty = 0;

	cin >> size;
	cin >> depth_to_search;
	board = new int*[size];
	for( int i = 0; i < size; i++){
		board[i] = new int[size];
	}
	for( int i = 0; i < size; i++){
		for( int j = 0; j < size; j++){
			board[i][j] = 0; 
		}
	}
	
	cin >> num_given_white;
	for(int i = 0; i < num_given_white; i++){
		cin >> x;
		cin >> y;
		board[x][y] = 1;
	}


	cin >> num_given_black;
	for( int i = 0; i < num_given_black; i++){
		cin >> x;
		cin >> y;
		board[x][y] = -1;
	}

	for( int k = 1; k <= depth_to_search; k++){
		int result;
		result = nega_max( 1, board, size, 0 , k, -1, 1 );
		
		if( result == -1 ) {
			cout << k << " " << "black";
			break;
		}
		else if( result == 1) {
			cout << k << " " << "white";
			break;
		}
		else if( k == depth_to_search )
			cout << "draw";
		
	}
	
	return 0;
}
