#include<iostream>
#include<fstream>
#include<vector>
#include<queue>
using namespace std;

//#define nodeNo 5
#define inf 999999999
void printMatrix(int**,int);




int computeRow(int** Matrix,int dim )
{
	//Find min in each row
	int* min=new int[dim];

	for(int i=0;i<dim;i++)
	{
		min[i]=inf;
		for(int j=0;j<dim;j++)
		{
			if(Matrix[i][j]<min[i])
				min[i]=Matrix[i][j];
		}
	

	}

	int sumMin=0;
	//Substract zero in each row
	for(int i=0;i<dim;i++)
	{
		for(int j=0;j<dim;j++)
		{
			if(Matrix[i][j]!=inf)
			Matrix[i][j]=Matrix[i][j]-min[i];
		}
		
		if(min[i]!=inf)
		sumMin+=min[i];
	}

	return sumMin;
}


int computeCol(int** Matrix,int dim)
{
	//Find min in each row
	int* min=new int[dim];
	for(int i=0;i<dim;i++)
	{
	min[i]=inf;
		for(int j=0;j<dim;j++)
		{
			if(Matrix[j][i]<min[i])
				min[i]=Matrix[j][i];

		}
	}


	int sumMin=0;
	//Substract zero in each row
	for(int i=0;i<dim;i++)
	{

		for(int j=0;j<dim;j++)
			if(Matrix[j][i]!=inf)
			   Matrix[j][i]=Matrix[j][i]-min[i];
		
		if(min[i]!=inf)
	    	sumMin+=min[i];
	}


	return sumMin;
}





int computeLbound(int **matrix,int dim)
{
	 
	int colVal=computeCol(matrix,dim);
	//cout<<"col"<<endl;
	//printMatrix(matrix,dim);


	int rowVal=computeRow(matrix,dim);
	//cout<<"row"<<endl;
	//printMatrix(matrix,dim);

	return (rowVal+colVal);
}



typedef struct
{
	int **Matrix;
	int currendDim;
	int *solution;
	int fromNode;
	int solutionDimAdd;
	int currentBound;
	int toNode;
	int *linkPath;
	vector<int> remainedSolution;
	vector<int> remainedRow;

} Record;


void printMatrix(int**M, int n)
{
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
			 //myfile << M[i][j]<<" ";
			cout<<M[i][j]<<" ";
		cout<<endl;
	}


//system("Pause");
}


int **deleteRowColumn(int**M,int delete_row,int delete_col,int n)
{
	int **newMatrix;
			
	newMatrix=new int*[n-1];
	int count1,count2;
	count1=0;count2=0;

	for(int i=0;i<n-1;i++)
	  newMatrix[i]=new int[n-1];

	for(int i=0;i<n;i++)
	{
	
		if(i!=delete_row)
		{
			for(int j=0;j<n;j++)
			{
				if(j!=delete_col)
				{
 					newMatrix[count1][count2]=M[i][j];
					count2++;
				}
			}
		
	    count1++;

		}
		count2=0;
	}
	return newMatrix;

}


int **branchLchild(int**M,int row,int col,int n)
{
	int **LMatrix;
	int countRow,countCol;
	int originalTour=M[col][row];

	LMatrix=new int*[n-1];

	for(int i=0;i<n-1;i++)
		LMatrix[i]=new int[n-1];

	
	//printMatrix(M,n);
//	system("Pause");

	//逆的路徑標示不能走
	M[col][row]=inf;
	
	countRow=0;
	for(int i=0;i<n;i++)
		if(i!=row)
		{
			countCol=0;
			for(int j=0;j<n;j++)
				if(j!=col)
				{
					LMatrix[countRow][countCol]=M[i][j];
					countCol++;
				}			
			countRow++;
		}
	

//	printMatrix(LMatrix,n-1);
//	system("Pause");

	M[col][row]=originalTour;
	return LMatrix;
}



int **branchRchild(int**M,int row,int col, int n)
{
	int **RMatrix;
	
	RMatrix=new int*[n];
	for(int i=0;i<n;i++)
		RMatrix[i]=new int[n];

	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			RMatrix[i][j]=M[i][j];

	RMatrix[row][col]=inf;

	return RMatrix;
}

int computeCost(int **M,int* solution,int n)
{
	int sum=0;
	for(int i=0;i<n;i++)
		sum+=M[solution[i]][solution[i+1]];	
	 
	return sum;
}


int computePenalty(int **M,int fromNode,int toNode,int n)
{
	//int penaltyMatrix;
	int maxCol=inf;
	int maxRow=inf;

			for(int j=0;j<n;j++)			
				if(M[j][toNode]<maxCol && j!=fromNode && M[j][toNode]!=inf )
					maxCol=M[j][toNode];					
	
			if(maxCol==inf)
				maxCol=0;

			for(int j=0;j<n;j++)
				if(M[fromNode][j]<maxRow && j!=toNode && M[fromNode][j]!=inf )
					maxRow=M[fromNode][j];	

			if(maxRow==inf)
				maxRow=0;
//	printMatrix(M,n);
//	cout<<fromNode<<","<<toNode<<endl<<maxCol+maxRow;
//	system("Pause");

	return (maxCol+maxRow);
}

bool checkInfeasible(int **M,int fromNode,int n)
{
	int countInf=0;

	for(int i=0;i<n;i++)
		if(M[fromNode][i]==inf)
			countInf++;
	
	if(countInf==n) 
		return 1;
	else
		return 0;


}

bool detectBadCycle(int toNode,int fromNode,int* list,int n)
{
	int count=0;
	int locate;
	bool isEnd=0;
	
 	locate=list[toNode];
	/*
	cout<<"detecting"<<endl;
	cout<<fromNode<<",";
	cout<< toNode<<","<<list[toNode]<<endl;
	for(int i=0;i<n;i++)
		cout<<list[i]<<",";
	cout<<endl<<endl;
	*/


	if(locate==fromNode)
		return 1;


	while(1)
	{
		if(locate==-1)			
			return 0;  //No cycle
		else
		{
			locate=list[locate];
			if(locate==fromNode)
				break;
		}	
	}

	return 1;
}

int detectCycleLength(int toNode,int fromNode,int* list,int n)
{
	int cycleNo=1;

	int count=0;
	int locate;
	bool isEnd=0;
	
 	locate=list[toNode];

	/*
	cout<<"detecting"<<endl;
	cout<<fromNode<<",";
	cout<< toNode<<","<<list[toNode]<<endl;
	for(int i=0;i<n;i++)
		cout<<list[i]<<",";
	cout<<endl<<endl;
	*/


	while(1)
	{
		cycleNo++;
		if(locate==-1)
		{
		return cycleNo;
			return 0;  //No cycle

		}
		else
		{
		
			locate=list[locate];
			if(locate==fromNode)
				break;
		}	
	}

	return cycleNo;
}

int main(int argc,char **argv)
{
	

	Record tempRecord;
	int **Matrix;
	int **penaltyMatrix;
	int *solution, *finalSolution;
	int *linkPath,*finalLinkPath;


	

	int currentLevel=0;
	int nextCandidateNode;
	int fromNode=0;
	int toNode;
	int solutionDimAdd=0;
	int currentBound=0;
	//input
	int nodeNumber;
	bool endSearch=0;
	int upperBound=inf;

	fstream fin;
	

	
	vector<Record> memRecord;
	vector<int> remainedSolution;
	vector<int> remainedRow;
	//vector<Record> memRecord;

	//vector<int**> memMatrix;
	//vector<int> memLbound;
	//vector<int*> memSolution;
	
	if(argc==1)
	{
	
	cin >>nodeNumber;
	}
	else
	{
	cin >>nodeNumber;
	}

	//	myfile.open ("output.txt");
	//讀取檔案(存入節點數量)
 
	
    int currentDim=nodeNumber;

	int** CostMatrix;


	
	linkPath=new int[nodeNumber];
	finalLinkPath=new int[nodeNumber];

	Matrix=new int*[nodeNumber];

    CostMatrix=new int*[nodeNumber];
	for(int i=0;i<nodeNumber;i++)
		CostMatrix[i]=new int[nodeNumber];

	//將檔案寫入costMatrix中
	for(int i=0;i<nodeNumber;i++)
		for(int j=0;j<nodeNumber;j++)
			cin>>CostMatrix[i][j];

	//對costMatrix對角元素設很大的數,防止自己節點走到自己節點
	for(int i=0;i<nodeNumber;i++)
	{
		linkPath[i]=-1;
		CostMatrix[i][i]=inf;
	}
	
	



	solution=new int[nodeNumber+1];
	finalSolution=new int[nodeNumber+1];
	 

	//Matrix<=CostMatrix...計算lower bound與 branching都在Matrix上做
	for(int i=0;i<nodeNumber;i++)
	 {
	  Matrix[i]=new int[nodeNumber];
		 for(int j=0;j<nodeNumber;j++)
			Matrix[i][j]=CostMatrix[i][j];

	 }
	 
	//存入剩下可以選的節點
	 for(int i=0;i<nodeNumber;i++)
	 {
		 remainedSolution.push_back(i);
		 remainedRow.push_back(i);

	 }
	solution=new int[nodeNumber];

	//計算trunk的 lower bound
	currentBound=computeLbound(Matrix,nodeNumber);
	
//	printMatrix(Matrix,nodeNumber);
	//system("Pause");

	/*for(int i=0;i<currentDim;i++)
		{
			for(int j=0;j<currentDim;j++)
				myfile<<Matrix[i][j]<<",";
			myfile<<endl;
		}
		myfile<<endl;
		*/


	//B&B starts

	
	solutionDimAdd=1;
	int colAdj=0;
	bool isTerminated=0;
	bool isInfeasible=0;


	while(isTerminated==false)	
	{
		
		penaltyMatrix=new int*[currentDim];
		
		for(int i=0;i<currentDim;i++)
		{
			penaltyMatrix[i]=new int[currentDim];
			for(int j=0;j<currentDim;j++)
				penaltyMatrix[i][j]=-1;

		}
	
		
		int MaxPenalty=-inf;
		vector<int> recordPenalty;
		vector<int> recordToNode;
		vector<int> recordFromNode;
		bool isAllCycle=1;
	
		if(currentBound<upperBound)
		for(int i=0;i<currentDim;i++)
			for(int j=0;j<currentDim;j++)
			{
												
				if(Matrix[i][j]==0)
				{
					
					if(detectBadCycle(remainedSolution.at(j),remainedRow.at( i),linkPath,nodeNumber) )
					{
						Matrix[i][j]=inf;
						
					}
					else
					{
						isAllCycle=0;

						int penalty=computePenalty(Matrix,i,j,currentDim);
						penaltyMatrix[i][j]=penalty;
						recordPenalty.push_back(penalty);
						recordToNode.push_back(j);
						recordFromNode.push_back(i);
						
						if(penalty>MaxPenalty)
						{
							MaxPenalty=penalty;
							fromNode=i;
							toNode=j;
						}
					}
				}
			}


			if(isAllCycle==1)
			{
			//Infeasible
			//isAllCycle=1;
		//	cout<<"penaltyM"<<endl;
		//	printMatrix(penaltyMatrix,currentDim);
		//	cout<<"Matrix"<<endl;
		//	printMatrix(Matrix,currentDim);
			//system("Pause");
			isInfeasible=1;
			endSearch=1;
			}


	//	printMatrix(penaltyMatrix,currentDim);
		//printMatrix(Matrix,currentDim);
	 //   system("Pause");
		

		/*
		for(int i=0;i<recordPenalty.size();i++)
			if(recordPenalty.at(i)==MaxPenalty && recordToNode.at(i)!=toNode)
			{
								
					Record tempRecord;			
					tempRecord.solutionDimAdd=solutionDimAdd;

					tempRecord.currendDim=currentDim;
					tempRecord.solution=new int[nodeNumber];
					tempRecord.currentBound=currentBound;
					tempRecord.fromNode=recordFromNode.at(i);
					tempRecord.toNode=recordToNode.at(i);
					tempRecord.linkPath=new int[nodeNumber];

					for(int j=0;j<currentDim-1;j++)					
						tempRecord.solution[j]=solution[j];
					
					for(int j=0;j<nodeNumber;j++)					
						tempRecord.linkPath[j]=linkPath[j];


					tempRecord.Matrix=new int*[currentDim];
					tempRecord.fromNode=recordFromNode.at(i);

					for(int j=0;j<currentDim;j++)
					{
						tempRecord.Matrix[j]=new int[currentDim];
						
						for(int k=0;k<currentDim;k++)
							tempRecord.Matrix[j][k]=Matrix[j][k];						
					}

					for(int j=0;j<remainedSolution.size();j++)
						tempRecord.remainedSolution.push_back(remainedSolution.at(j));
					
					for(int j=0;j<remainedRow.size();j++)
						tempRecord.remainedRow.push_back(remainedRow.at(j));

					
					for(int j=0;j<recordToNode.size();j++)
						if(j!=i)
							tempRecord.Matrix[fromNode][recordToNode.at(j)]=inf;
										
					memRecord.push_back(tempRecord);
			}
			
			recordPenalty.clear();
	
		*/



	//	int CN=detectCycleLength(remainedSolution.at(toNode),remainedRow.at(fromNode),linkPath,nodeNumber);
	//	myfile<<endl<<"cycleNo:"<< CN<<endl;
//		myfile<<"cyclic?"<<detectBadCycle(remainedSolution.at(toNode),remainedRow.at(fromNode),linkPath,nodeNumber)<<endl;

		//myfile<<endl<<currentBound<<endl;
		
   		int **LeftChildMatrix,**RightChildMatrix,leftBound,rightBound;

		if(endSearch==false)
		{
		 LeftChildMatrix= branchLchild(Matrix,fromNode,toNode,currentDim);
		 RightChildMatrix= branchRchild(Matrix,fromNode,toNode,currentDim);
	
	//	 cout<<"leftMatrix"<<endl;
	//	 printMatrix(LeftChildMatrix,currentDim-1);
		
		 leftBound=computeLbound(LeftChildMatrix,currentDim-1);
		 rightBound=computeLbound(RightChildMatrix,currentDim);
		}




		
		if(endSearch==false && isInfeasible==false)
 		if(leftBound<=rightBound)
		{ 
			//進入左節點
		//	cout<<endl<<"selected"<<endl;
		//	cout<<remainedRow.at(fromNode)<<","<<remainedSolution.at(toNode)<<endl;
			//system("Pause");

	
			//儲存右節點
			
					Record tempRecord;
					tempRecord.solutionDimAdd=solutionDimAdd;
					tempRecord.currendDim=currentDim;
					tempRecord.solution=new int[nodeNumber];
				
					tempRecord.linkPath=new int[nodeNumber];

					for(int j=0;j<currentDim-1;j++)					
						tempRecord.solution[j]=solution[j];
					
					
					
					for(int j=0;j<nodeNumber;j++)
						tempRecord.linkPath[j]=linkPath[j];


					tempRecord.Matrix=new int*[currentDim];
					
					for(int j=0;j<currentDim;j++)
						tempRecord.Matrix[j]=new int[currentDim];					
					
					tempRecord.fromNode=remainedRow.at(fromNode);
					tempRecord.toNode=remainedSolution.at(toNode);
					tempRecord.currentBound=rightBound+currentBound;
					

					int count1=0;
					int count2=0;
					

					
					for(int j=0;j<currentDim;j++)
						for(int k=0;k<currentDim;k++)													
							tempRecord.Matrix[j][k]=RightChildMatrix[j][k];																							
					
					

					for(int j=0;j<remainedSolution.size();j++)
						tempRecord.remainedSolution.push_back(remainedSolution.at(j));

					for(int j=0;j<remainedRow.size();j++)
						tempRecord.remainedRow.push_back(remainedRow.at(j));
			
					/*
					for(int j=0;j<currentDim;j++)
					{
						for(int k=0;k<currentDim;k++)													
							cout<<tempRecord.Matrix[j][k]<<" ";
					cout<<endl;
					}*/


					memRecord.push_back(tempRecord);
					
			
			

			//進入下個state的參數調整
			currentBound+=leftBound;
			
			solution[solutionDimAdd]= remainedSolution.at(toNode);
			linkPath[remainedRow.at(fromNode)]=remainedSolution.at(toNode);

			remainedSolution.erase(remainedSolution.begin()+toNode);
			remainedRow.erase(remainedRow.begin()+fromNode);
				 	
			
			currentDim--;

			Matrix=LeftChildMatrix;
			
			/*
					for(int j=0;j<currentDim+1;j++)
					{
						for(int k=0;k<currentDim+1;k++)													
							cout<<tempRecord.Matrix[j][k]<<" ";
					cout<<endl;
					}*/
			solutionDimAdd++;

		}
		else
		{
		//進入右節點
		
			//儲存左節點
			Record tempRecord;
			tempRecord.solutionDimAdd=solutionDimAdd+1;
			tempRecord.currendDim=currentDim-1;
			tempRecord.solution=new int[nodeNumber];
			tempRecord.currentBound=currentBound+leftBound;
	

			tempRecord.linkPath=new int[nodeNumber];

			for(int j=0;j<currentDim-1;j++)					
				tempRecord.solution[j]=solution[j];
					
					
					
			for(int j=0;j<currentDim;j++)
				tempRecord.linkPath[j]=linkPath[j];

			tempRecord.linkPath[remainedRow.at(fromNode)]=remainedSolution.at(toNode);

			tempRecord.Matrix=new int*[currentDim-1];
					
			for(int j=0;j<currentDim;j++)
				tempRecord.Matrix[j]=new int[currentDim-1];					
					
				tempRecord.fromNode=remainedRow.at(fromNode);
				tempRecord.toNode=remainedSolution.at(toNode);
				
					

				
					for(int j=0;j<currentDim-1;j++)
						for(int k=0;k<currentDim-1;k++)													
							tempRecord.Matrix[j][k]=LeftChildMatrix[j][k];																							
				
					/*
					cout<<"LeftChildMatrix"<<endl;
					printMatrix(LeftChildMatrix,currentDim-1);
					cout<<"tempChildMatrix"<<endl;
					printMatrix(tempRecord.Matrix,currentDim-1);
				
					system("Pause");
				
					*/
					for(int j=0;j<remainedSolution.size();j++)
						tempRecord.remainedSolution.push_back(remainedSolution.at(j));
					
					tempRecord.remainedSolution.erase(tempRecord.remainedSolution.begin()+toNode);

					
					for(int j=0;j<remainedRow.size();j++)
						tempRecord.remainedRow.push_back(remainedRow.at(j));
			
					tempRecord.remainedRow.erase(tempRecord.remainedRow.begin()+fromNode);



  					memRecord.push_back(tempRecord);

					



			//進入下個state的參數調整
			currentBound+=rightBound;
			Matrix=RightChildMatrix;
		}

//		myfile<<endl<<currentBound<<endl;

/*
		for(int i=0;i<currentDim;i++)
		{
			for(int j=0;j<currentDim;j++)
				myfile<<Matrix[i][j]<<",";
//			myfile<<endl;
		}*/
	    
		if(endSearch==false)
		if(currentBound>=upperBound)
		{
			endSearch=1;

		}
	
		




		if(solutionDimAdd==(nodeNumber))
		{
			endSearch=1;

			

			if(detectCycleLength(remainedSolution.at(0),remainedRow.at(0),linkPath,nodeNumber)< (nodeNumber-1))
			{
				cout<< detectCycleLength(remainedSolution.at(0),remainedRow.at(0),linkPath,nodeNumber);
				printf("\nBad Cycle Occcurred!!!");
			}
			
			
			for(int i=0;i<nodeNumber;i++)
				finalSolution[i]=solution[i];
		
			
			upperBound=currentBound;

			

			linkPath[remainedRow.at(0)]=remainedSolution.at(0);
			

			for(int i=0;i<nodeNumber;i++)
				finalLinkPath[i]=linkPath[i];


			int locateCount=1;
			int locate=0;
			finalSolution[0]=0;
			for(int i=1;i<nodeNumber+1;i++)
			{
				
		       

				finalSolution[locateCount]=linkPath[locate];
				locate=linkPath[locate];
				//finalSolution[locateCount+1]=linkPath[locateCount];
				

				

				locateCount++;
			}

			//for(int i=0;i<nodeNumber+1;i++)
				//cout<<finalSolution[i]+1<<"=>";
			//cout<<endl;

			//需調整
			//isTerminated=1;
		}

//		cout<<"selectedMatrix"<<endl;
//			printMatrix(Matrix,currentDim);

		/*
		if(endSearch==false)
		{
			myfile<<"leftB:"<<leftBound<<" ";
//	myfile<<"rightB:"<<rightBound<<endl;
		}
		
		for(int i=0;i<nodeNumber;i++)
		myfile<<linkPath[i]<<" ";
		*/




	//	myfile<<endl;
		//cout<<currentBound<<endl;
	//	myfile<<endl;



		if(endSearch==1)
		{
   			if(memRecord.size()==0)
				isTerminated=1;
			else
			{
				int memIndex=0;
				
				solutionDimAdd=memRecord.at(memIndex).solutionDimAdd;
				currentDim=memRecord.at(memIndex).currendDim;
				currentBound=memRecord.at(memIndex).currentBound;
				fromNode=memRecord.at(memIndex).fromNode;
				toNode=memRecord.at(memIndex).toNode;

				Matrix=new int*[currentDim];
				for(int i=0;i<currentDim;i++)
					Matrix[i]=new int[currentDim];


				remainedSolution.clear();
				remainedRow.clear();

				for(int i=0;i<memRecord.at(memIndex).remainedSolution.size();i++)
				   remainedSolution.push_back(memRecord.at(memIndex).remainedSolution.at(i));

				for(int i=0;i<memRecord.at(memIndex).remainedRow.size();i++)
				   remainedRow.push_back(memRecord.at(memIndex).remainedRow.at(i));

				//for(int i=0;i<memRecord.at(memIndex).currendDim;i++)
				for(int i=0;i<nodeNumber;i++)
					linkPath[i]=memRecord.at(memIndex).linkPath[i];


	//			myfile<<"recallMatrix"<<endl;
				for(int i=0;i<currentDim;i++)
				{

					//solution[i]=memRecord.at(memIndex).solution[i];
					
					for(int j=0;j<currentDim;j++)
					{
    						Matrix[i][j]=memRecord.at(memIndex).Matrix[i][j];
							


					}
//					myfile<<endl;
				}
		



				memRecord.erase(memRecord.begin()+memIndex);
				endSearch=0;
				isInfeasible=0;



				if(memRecord.size()==0)
					isTerminated=1;

			}
			
			
		}
	

	
	}



	int cost=computeCost(CostMatrix,finalSolution,nodeNumber);

	cout<<"best rout:"<<endl;
	
	for(int i=0;i<nodeNumber+1;i++)
		cout<<finalSolution[i]+1<<" ";


	cout<<endl<<"cost:"<<cost<<endl;
	
//	 myfile.close();
//	system("Pause");



	return 0;
}
