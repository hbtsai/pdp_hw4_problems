#include <iostream>
#include <string.h>
#include <math.h>
#define N 3
using namespace std;

int G[9][9]={0};
bool ans=0;

bool check();
void backtrace(int,int);
void output();

int main()
{
    bool sw=0;
	ans = 0;
	memset(G,0,sizeof(G));
	for(int i=0;i<N*N;i++)
		for(int j=0;j<N*N;j++)
			cin >> G[i][j];
	if(sw){  cout << endl;}
	else{    sw = 1;}
	backtrace(0,0);
	if(!ans)
            cout << "NO SOLUTION" << endl;
    return 0;
}

bool check()
{
    bool sw=1;
    for(int i = 0;i<N*N;i++)
    {
        int tmp1[10] = {0};
        int tmp2[10] = {0};
        int tmp3[10] = {0};
        for(int j=0;j<N*N;j++)
        {
            tmp1[ G[i][j] ]++;
            tmp2[ G[j][i] ]++;
        }
        for(int j=1;j<=N*N;j++)
            if(tmp1[j] > 1 || tmp2[j] > 1)
            {
                sw=0;
                break;
            }
        if(!sw) break;
    }
    if(sw)
    {
        int tmp[3][3][10]={0};
        for(int i=0;i<N*N;i++)
            for(int j=0;j<N*N;j++)
                tmp[i/N][j/N][ G[i][j] ]++;
        for(int i=0;i<N && sw;i++)
            for(int j=0;j<N && sw;j++)
                for(int k=1;k<=N*N;k++)
                    if(tmp[i][j][k] > 1)
                    {
                        sw=0;
                        break;
                    }
    }
    return sw;
}

void backtrace(int x,int y)
{
    if(ans)
        return;
    if(y == N*N && x == N*N-1)
    {
        ans = 1;
        output();
        return;
    }
    if( y == N*N )
        return backtrace(x+1,0);
    if(G[x][y])
        return backtrace(x,y+1);
    for(int i=1;i<=N*N;i++)
    {
        G[x][y] = i;
        if(check())
        {
            backtrace(x,y+1);
            G[x][y] = 0;
        }
    }
    G[x][y] = 0;
}

void output()
{
    for(int i=0,sw=0;i<N*N;i++,cout << endl,sw=0)
        for(int j=0;j<N*N;j++)
        {
            if(sw) cout << " ";
            else sw=1;
            cout << G[i][j];
        }
}
