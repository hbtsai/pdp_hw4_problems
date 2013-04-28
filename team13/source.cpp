#include <cstdio>
#include <cstring>
#include <vector>
#include <queue>

using namespace std;

#define COIN_MAX 50
#define CITY_MAX 1001
#define FEE_MAX 1000001
#define INF 0x7F7F7F7F

struct PATH{
    int des, fee;
};

int start, dest, coin, coin_value[COIN_MAX], coin_dp[FEE_MAX], cost[CITY_MAX];
vector<PATH> path[CITY_MAX];

void read_input(void){
    scanf("%d", &coin);
    for(int i=0; i<coin; ++i)
        scanf("%d", &coin_value[i]);
    int m, n;
    scanf("%d %d", &m, &n);
    for(int i=1; i<=m; ++i)
        path[i].clear();
    for(int i=0; i<n; ++i){
        int from, to;   PATH tmp;
        scanf("%d %d %d", &from, &to, &tmp.fee);
        // bidirectional path
        tmp.des = to;   path[from].push_back(tmp);
        tmp.des = from; path[to].push_back(tmp);
    }
    scanf("%d %d", &start, &dest);
}

int coin_change(void){
    // dp initialization
    memset(coin_dp, 127, sizeof(coin_dp));
    coin_dp[0] = 0;
    // dynamic programming
    for(int i=0; i<coin; ++i){
        for(int j=coin_value[i]; j<FEE_MAX; ++j){
            int tmp = coin_dp[j-coin_value[i]] + 1;
            if( coin_dp[j] > tmp )
                coin_dp[j] = tmp;
        }
    }
}

void bellman_ford(void){
    // cost initialization
    memset(cost, 127, sizeof(cost));
    cost[start] = 0;
    // go shortest path
    queue<int> que; que.push(start);
    while( !que.empty() ){
        int tmp = que.front();  que.pop();
        for(int i=0; i<path[tmp].size(); ++i){
            int new_cost = cost[tmp] + coin_dp[path[tmp][i].fee];
            if( cost[path[tmp][i].des] > new_cost ){
                cost[path[tmp][i].des] = new_cost;
                que.push(path[tmp][i].des);
            }
        }
    }
}

int main(void){
    
    read_input();
    coin_change();
    bellman_ford();
    printf("%d\n", (cost[dest]==INF)?(-1):(cost[dest]));
    
    return 0;
}
