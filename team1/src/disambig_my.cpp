#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Ngram.h"
#define NUMofZhuYin 37
#define MAXWordInLine 100
double getBigramProb(const char *w1, const char *w2, Vocab &voc, Ngram &lm);
void getZhuYin( char ZhuYins[NUMofZhuYin][4], FILE *ZhuYin);
int isZhuYin( char *word, char ZhuYins[NUMofZhuYin][4]);
double P[MAXWordInLine][3000];
int who[MAXWordInLine][3000];
int main( int argc, char *argv[]){
	//intput: (map) (ZhuYin) (bigram.lm)
	FILE *ZhuYin_Big5;
	FILE *ZhuYin;
	char word[MAXWordInLine][3000][10];
	int NUMword[MAXWordInLine];
	int ans_sequence[MAXWordInLine];
	char string[MAXWordInLine][10];
	char *tokenptr;
	char ZhuYins[NUMofZhuYin][4];
	int count;
	double max;
	double max_temp;
	int who_give;
	int i, j, k;
	char temp[100000];

	Vocab voc;
	Ngram lm( voc, 2);
	File lmFile(argv[3], "r");
	lm.read(lmFile);
	lmFile.close();

	ZhuYin_Big5 = fopen( argv[1], "r");
	ZhuYin = fopen( argv[2], "r");
	if(ZhuYin_Big5 == NULL || ZhuYin == NULL){
		printf("file open error\n");
		fclose(ZhuYin_Big5);
		fclose(ZhuYin);
		return 0;
	}
	getZhuYin( ZhuYins, ZhuYin);

	while(fgets( temp, 100000, stdin) != NULL){
		//printf("%s", temp);
		count = 0;
		for( i = 0; i < MAXWordInLine; i++){
			NUMword[i] = 0;
		}
		tokenptr = strtok( temp, " \n");
		while( tokenptr != NULL){
			strcpy(string[count], tokenptr);
			count++;
			tokenptr = strtok( NULL, " \n");
		}

		for( i = 0; i < count; i++){
			if( isZhuYin( string[i], ZhuYins) == 1){
				fseek( ZhuYin_Big5, 0, SEEK_SET);
				while( fgets( temp, 100000, ZhuYin_Big5) != NULL){
					tokenptr = strtok( temp, " \t");
					if(strcmp( string[i], tokenptr) == 0){
						tokenptr = strtok( NULL, " \n");
						while(tokenptr != NULL){
							strcpy( word[i][NUMword[i]], tokenptr);
							NUMword[i]++;
							tokenptr = strtok( NULL, " \n");
						}
						break;
					}
				}
			}
			else{
				strcpy( word[i][0], string[i]);
				NUMword[i]++;
			}
		}
		for( i = 0; i < count; i++){
			for(j = 0; j < NUMword[i]; j++){
				if( i == 0){
					max = 0;
					who_give = 0;
					max_temp = getBigramProb( "<s>", word[i][j], voc, lm);
					if( max == 0){
						max = max_temp;
						who_give = j;
					}
					if( max_temp > max){
						max = max_temp;
						who_give = j;
					}
					P[i][j] = max;
					who[i][j] = who_give;
				}
				else{
					max = 0;
					who_give = 0;
					for( k = 0; k < NUMword[i - 1]; k++){
						max_temp = getBigramProb( word[i - 1][k], word[i][j], voc, lm);
						max_temp = max_temp + P[i - 1][k];
						if( max == 0){
							max = max_temp;
						}
						if( max_temp > max){
							max = max_temp;
							who_give = k;
						}
					}
					P[i][j] = max;
					who[i][j] = who_give;
				} 
			}
		}
		max = 0;
		who_give = 0;
		for( j = 0; j < NUMword[count - 1]; j++){
			if( max == 0){
				max = P[count - 1][j];
				who_give = j;
			}
			else{
				if( P[count - 1][j] > max){
					max = P[count - 1][j];
					who_give = j;
				}
			}
		}
		ans_sequence[count - 1] = who_give; 
		for( i = count - 1; i > 0; i--){
			ans_sequence[i - 1] = who[i][ans_sequence[i]];
		}
		printf("<s> ");
		for( i = 0; i < count; i++){
			printf("%s ", word[i][ans_sequence[i]]);
		}
		printf("</s>\n");
	}
	fclose(ZhuYin_Big5);
	fclose(ZhuYin);
	return 0;
}

//Get P(W2 | W1) -- bigram
double getBigramProb(const char *w1, const char *w2, Vocab &voc, Ngram &lm){
	VocabIndex wid1 = voc.getIndex(w1);
	VocabIndex wid2 = voc.getIndex(w2);
	if(wid1 == Vocab_None)  //OOV
		wid1 = voc.getIndex(Vocab_Unknown);
	if(wid2 == Vocab_None){  //OOV
		wid2 = voc.getIndex(Vocab_Unknown);
		return -20;
	}
	VocabIndex context[] = { wid1, Vocab_None };
	return lm.wordProb( wid2, context);
}

//把所有注音符號填入陣列中
void getZhuYin( char ZhuYins[NUMofZhuYin][4], FILE *ZhuYin){
	int i;
	char temp[4];
	for( i = 0; i < NUMofZhuYin; i++){
		fscanf( ZhuYin, "%s", temp);
		strcpy( ZhuYins[i], temp);
	}
	return;
}

//判斷是不是注音
int isZhuYin( char *word, char ZhuYins[NUMofZhuYin][4]){
	int i;
	for( i = 0; i < NUMofZhuYin; i++){
		if(strcmp( word, ZhuYins[i]) == 0)
			return 1;
	}
	return 0;
}
