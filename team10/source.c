#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int n,m,k;
char **equation;
int varible[26];

void print_equation(){
	int k_i;
	int sl,sl_i;
	for( k_i = 0 ; k_i < k ; k_i++ ){
		sl = strlen(equation[k_i]);
		for( sl_i = 0 ; sl_i < sl ; sl_i++ ){
			if( equation[k_i][sl_i] >= 'A' && equation[k_i][sl_i] <= 'Z' ){
				printf("(%d)",varible[equation[k_i][sl_i]-'A']);
			}
			else{
				printf("%c",equation[k_i][sl_i]);
			}
		}
		printf("\n");
	}
	return;
}

int check_equation(){
	int k_i;
	int sl,sl_i;
	int left,right,add_temp,mul_temp,temp;
	for( k_i = 0 ; k_i < k ; k_i++ ){
		temp = 0;
		mul_temp = 1;
		add_temp = 0;
		sl = strlen(equation[k_i]);
		for( sl_i = 0 ; sl_i < sl ; sl_i++ ){
			if( equation[k_i][sl_i] >= 'A' && equation[k_i][sl_i] <= 'Z' ){
				temp *= n;
				temp += varible[equation[k_i][sl_i] - 'A'];
			}
			if( equation[k_i][sl_i] == '*' ){
				mul_temp *= temp;
				temp = 0;
			}
			if( equation[k_i][sl_i] == '+' ){
				mul_temp *= temp;
				add_temp += mul_temp;
				temp = 0;
				mul_temp = 1;
			}
			if( equation[k_i][sl_i] == '=' ){
				mul_temp *= temp;
				add_temp += mul_temp;
				left = add_temp;
				temp = 0;
				mul_temp = 1;
				add_temp = 0;
			}
		}
		mul_temp *= temp;
		add_temp += mul_temp;
		right = add_temp;
		if( right != left){
			return 0;
		}
	}

	return 1;
}

void generate_varibles(long long a){
	long long temp_a;
	temp_a = a;
	int m_i,m_j;
	for( m_i = 0 ; m_i < m ; m_i++ ){
		varible[m-1-m_i] = temp_a%(n-m+m_i+1);
		temp_a = temp_a/(n-m+m_i+1);
	}
	for( m_i = m-2 ; m_i >= 0 ; m_i-- ){
		for( m_j = m_i+1 ; m_j < m ; m_j++ ){
			if( varible[m_i] <= varible[m_j] ){
				varible[m_j]++;
			}
		}
	}
	return;
}

void generate_next(long long a){
	int n_i,m_i;
	for( n_i = varible[m-1]+1 ; n_i < n ; n_i++ ){
		for( m_i = 0 ; m_i < m-1 ; m_i++ ){
			if( n_i == varible[m_i] ){
				break;
			}
		}
		if( m_i >= m-1 ){
			varible[m-1] = n_i;
			return;
		}
	}
	if( n_i >= n ){
		generate_varibles(a);
	}
	return;
}

int main(){
	//read input
	{
		scanf("%d %d %d",&n,&m,&k);
		int k_i;
		equation = (char **)malloc(sizeof(long long)*k);
		for( k_i = 0 ; k_i < k ; k_i++ ){
			equation[k_i] = (char *)malloc(sizeof(char)*256);
			scanf("%s",equation[k_i]);
		}
	}

	//compute all case:a
	long long a = 1;
	{
		int m_i;
		for( m_i = n ; m_i > n-m ; m_i-- ){
			a *= m_i;
		}
	}

	//computing
	int s = 0;
	long long a_i,a_j = -100;
	for( a_i = 0 ; a_i < a ; a_i++ ){
		//generate varibles
		if( a_j == a_i-1 ){
			generate_next(a_i);
		}
		else{
			generate_varibles(a_i);
		}
		a_j = a_i;

		//check & record
		//print_equation();
		if( check_equation() ){//1:if right answer,0:if wrong answer
			s++;
		}

	}

	printf("%d",s);

	return 0;
}
