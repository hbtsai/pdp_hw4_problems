/*
 * main.c
 *
 *  Created on: 2013/4/18
 *      Author: freeze
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

char redirect_list[30],incremental_list[30];
unsigned char *prime_marker;
long long int boundary;

void init_mapping_indicator(char* array);
void init_increment_helper(char* array);
static void flip_marker(long long int virtual_addr);
static int is_prime(long long int virtual_addr);
static void clr_prime(long long int virtual_addr) ;
void init_prime_list(long long int maximan);
long long int next_prime(long long int current_prime);

int main() {
	unsigned long long int number,fraction,prime,number_sqrt,counter;
	int iterate;
	scanf("%llu", &number);
	number_sqrt=sqrt((double)number)+1;
	init_prime_list(number_sqrt);
	init_mapping_indicator(redirect_list);
	fraction=number;
	prime=next_prime(0);
	counter=0;
	while (prime*prime<=fraction) {
		if(fraction%prime==0){
			do{
				counter++;
				fraction/=prime;
			}while(fraction%prime==0);
			if(counter>1)
				printf("%lld %lld\n",prime,counter);
			else
				printf("%lld\n",prime);
			counter=0;
		}
		prime=next_prime(prime);
		if(prime==0)
			break;
	}
	if(fraction!=1)
		printf("%lld\n",fraction);
	return 0;
}
void init_prime_list(long long int maximan) {
	int remainder;
	long long int sequence, boundary_sqrt, flip_target;
	//modified boundary number
	boundary_sqrt = (long long int) sqrt((double) maximan) + 1;
	remainder = boundary_sqrt % 30;
	if (remainder != 0) {
		boundary_sqrt = boundary_sqrt - remainder + 30;
	}
	boundary = boundary_sqrt * boundary_sqrt;
	init_mapping_indicator(redirect_list);
	init_increment_helper(incremental_list);
	prime_marker = (unsigned char*) calloc((boundary / 30),
			sizeof(unsigned char));
	for (long long int i = 0; i < boundary_sqrt; i++) {
		if (4 * i * i >= boundary)
			break;
		for (long long int j = 0; j < boundary_sqrt; j++) {
			/*case1
			 *if 4*i*i+j*j is legal number, flip prime indicator
			 */
			flip_target = 4 * i * i + j * j;
			remainder = flip_target % 12;
			if (flip_target >= boundary)
				break;
			if (remainder == 1 || remainder == 5) {
				flip_marker(flip_target);
			}
		}
	}
	for (long long int i = 0; i < boundary_sqrt; i++) {
		if (3 * i * i >= boundary)
			break;
		for (long long int j = 0; j < boundary_sqrt; j++) {
			/*case2
			 *if 3*i*i+j*j is legal number, flip prime indicator
			 */
			flip_target = 3 * i * i + j * j;
			if (flip_target >= boundary)
				break;
			remainder = flip_target % 12;
			if (remainder == 7) {
				flip_marker(flip_target);
			}
		}
	}
	for (long long int i = 0; i < boundary_sqrt; i++) {
		if (2 * i * i >= boundary)
			break;
		for (long long int j = i-1; j !=-1; j--) {
			/*case2
			 *if i>j and  3*i*i-j*j is legal number, flip prime indicator
			 */
			flip_target = 3 * i * i - j * j;
			remainder = flip_target % 12;
			if (flip_target >= boundary)
				break;
			if (remainder == 11) {
				flip_marker(flip_target);
			}
		}
	}
	//start to create prime number list
	for (long long int i = 6; i < boundary_sqrt; i+=incremental_list[i%30]) {
		if (is_prime(i)) {
			for (long long int squart = i * i, j = squart; j < boundary; j +=
					squart) {
				clr_prime(j);
			}
		}
	}

}
long long int next_prime(long long int current_prime) {
	switch (current_prime) {
	case 0:
		return 2;
	case 2:
		return 3;
	case 3:
		return 5;
	}
	for (long long int i = current_prime + 1; i < boundary; i+=incremental_list[i%30]) {
		if (is_prime(i)) {
			return i;
		}
	}
	return 0;
}

void init_mapping_indicator(char* array) {
	memset(array, 0, 30);
	for (int i = 0; i < 30; i += 2) {
		array[i] = -1;
	}
	for (int i = 0; i < 30; i += 3) {
		array[i] = -1;
	}
	for (int i = 0; i < 30; i += 5) {
		array[i] = -1;
	}
	for (int i = 0, j = 0; i < 30; i++) {
		if (array[i] == 0)
			array[i] = j++;
	}
}
void init_increment_helper(char* array) {
	memset(array, 1, 30);
	array[1]=6;
	array[7]=4;
	array[11]=2;
	array[13]=4;
	array[17]=2;
	array[19]=4;
	array[23]=6;
	array[29]=2;
}
static void flip_marker(long long int virtual_addr) {
	long long int phys_addr, remainder;
	remainder = virtual_addr % 30;
	if (redirect_list[remainder] == -1) {
		return;
	}
	phys_addr = (virtual_addr / 30);
	prime_marker[phys_addr] ^= 0x01 << redirect_list[remainder];

}

static int is_prime(long long int virtual_addr) {
	long long int phys_addr, remainder;
	remainder = virtual_addr % 30;
	if (redirect_list[remainder] == -1) {
		return 0;
	}
	phys_addr = (virtual_addr / 30);
	return (prime_marker[phys_addr] >> redirect_list[remainder]) & 0x01;
}
static void clr_prime(long long int virtual_addr) {
	long long int phys_addr, remainder;
	remainder = virtual_addr % 30;
	if (redirect_list[remainder] == -1) {
		return;
	}
	phys_addr = (virtual_addr / 30);
	prime_marker[phys_addr] &= ~(0x01 << redirect_list[remainder]);
}
