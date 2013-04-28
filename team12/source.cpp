#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include <algorithm>
#define dim 1200
    int source[dim][dim];//matrix used to store graph information
    int ans[10][dim];//matrix used to store answer
    int matrix_n;
    using namespace std;
int Search(int i,int j, int group);
int main(int argc, char** argv)
{ 
  
    int group_n;
    
    cin>>matrix_n;//read dimension of the matrix
    cin>>group_n;//read the number of groups
    
    int group[group_n][10];//matrix used to store group information
    int length[group_n];//number of elements in each group
    
    //read the number of elements in each group
    for (int count = 0; count < group_n; count++)
        cin>>length[count];  
        
      
    //read rule
    for(int i = 0; i< group_n; i++)
      for(int j = 0; j< length[i]; j++)
        cin>>group[i][j];
         
    //read the graph into source[][]
    for(int i=0;i<matrix_n;i++)
      for(int j=0;j<matrix_n;j++)
         cin>>source[i][j];
        


   //clean ans array   
   for(int i=0;i<group_n;i++)
      for(int j=0;j<matrix_n;j++)
        {
            ans[i][j]=0;
        }//end of for j     

    //convert the graph into individual group
    for(int i=0;i<matrix_n;i++)
    {
        for(int j=0;j<matrix_n;j++)
        {
            int t = 0;
            for(int u=0;(t==0)&(u<group_n); u++)
            {
                for(int v=0;(t==0)&(v < length[u]);v++)
                {
                    if(source[i][j] == group[u][v]) 
                    {
                       source[i][j] = u;//make the bit be the group number
                       t = t+1;
                    }//end of if
                }//end of for u
            }//end of for v
        }//end of for j 
    }

    /* For you to check
   //used to print the component
   for(int i=0;i<matrix_n;i++)
   { for(int j=0;j<matrix_n;j++)
        {
        printf("%d ",source[i][j]);
        }//end of for j
       printf("\n");
   }//end of for
   */
   
    for(int i=0;i<group_n;i++)
        length[i]=0;
        
    //use recursive method to solve the problem
    for(int i=0;i<matrix_n;i++)
      for(int j=0;j<matrix_n;j++)
        {   
            //if pixel haven't been processed
            if (source[i][j]>=0) 
            {
            int group = source[i][j];
            int count = Search( i, j, group);
            
            //record the final result
            int r = length[group];//used to store length of the answer
            ans[group][r] = count;
            length[group]=r+1;
            }//end of if
        }

    
    for(int u=0; u<group_n; u++)
    {
       printf("group %d has followings:\n",u);
       sort(ans[u],ans[u]+length[u]);
       for(int v=0; v<length[u]; v++)
        { 
          printf("%d ",ans[u][v]);
        }//end of for(v)
    
       printf("\n");
    }//end of for u
    
    return 0;
    
}//end of main


int Search(int i,int j,int group)
{
    //used for count
    int count=0;
    
    //check whether the index is reasonable or not
    if( i<0 || i>= matrix_n ) return 0;
    if( j<0 || j>= matrix_n ) return 0;

    //if group number is not what we want
    if(group != source[i][j] ) 
    {return 0;}
    
    //if source[i][j]<0, meaning that pixel has been processed
    if(source[i][j]<0) 
    {return 0;
    }//else, this pixel has not been processed
    else if(group == source[i][j])
    {
    count = count + 1;
    source[i][j] = -1;//processed
    count = count + Search( i-1, j-1,  group);
    count = count + Search( i-1, j,  group);
    count = count + Search( i-1, j+1,  group);
    count = count + Search( i, j-1,  group);
    count = count + Search( i, j+1,  group);
    count = count + Search( i+1, j-1,  group);
    count = count + Search( i+1, j,  group);
    count = count + Search( i+1, j+1,  group);
    return count;
    }else{return 0;}
    
}//end of search
