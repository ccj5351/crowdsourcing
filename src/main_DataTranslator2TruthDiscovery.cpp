/**
* @file: main_DataTranslator2TruthDiscovery.cpp
* @brief:
* @author: Changjiang Cai, ccai1@stevens.edu, caicj5351@gmail.com
* @version: 0.0.1
* @creation date: 17-12-2015
* @last modified: Thu 10 Mar 2016 09:32:02 AM EST
*/

#define _Main_DataTranslator2TruthDiscovery_CPP_
#ifndef _Main_DataTranslator2TruthDiscovery_CPP_


#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

int buff_num[4096], n, a[5];
int task_list[50000][2];
int vote[2000][2000][2];
int vote_n[2000];

int parse(char *buff, int *buff_num){
	int p = 1,i;
	for (i = 0; i < strlen(buff); i++){
		if (buff[i] == '\n')
			return p-1;
		if (buff[i] >= '0' && buff[i] <= '9'){
			buff_num[p] *= 10;
			buff_num[p] += buff[i] - '0';
		}
		else{
			p++;
		}
	}
	return p - 1;
}

int main(int argc, char * argv[]){
	int i, j;
	std::string inFileName =  argv[1];
	std::string outFileName = argv[2];
	std::string s_client_n =  argv[3];
	FILE *f = fopen(inFileName.c_str(), "r");
	if (f == NULL){
		printf("File not found!\n");
		return 0;
	}
	else{
		printf("Data set found!\n");
	}
	int obj_n, task_n, client_n;
	fscanf(f, "%d ", &obj_n);

	char *buff = (char*)malloc(sizeof(char)* 4096);
	client_n = stoi(s_client_n);

	while (1){
		fgets(buff, 4096, f);
		/*if (buff[0] == 'W' && buff[1] == 'o' && buff[2] == 'r')
			client_n++;*/
		if (buff[0] == 'T' && buff[1] == 'h' && buff[2] == 'e' && buff[3] == ' ' && buff[4] == 'r')
			break;
	}

	task_n = 0;

	for (i = 1; i <= client_n; i++){
		memset(buff_num, 0, sizeof(buff_num));
		fgets(buff, 4096, f);
		n = parse(buff, buff_num);
		if (n % 4 != 0){
			printf("Error!\n");
			return 1;
		}
		vote_n[i] = n / 4;
		for (j = 1; j <= vote_n[i]; j++){
			a[1] = buff_num[j * 4 - 3]+1;
			a[2] = buff_num[j * 4 - 2]+1;
			a[3] = buff_num[j * 4 - 1]+1;
			a[4] = buff_num[j * 4]+1;
			if (a[1] > task_n)
				task_n = a[1];
			task_list[a[1]][0] = a[2];
			task_list[a[1]][1] = a[3];
			vote[i][j][0] = a[1];
			vote[i][j][1] = a[4];
			/*if (i == 49 && (j == 1 || j == 16))
				printf("%d, %d, %d, %d", a[1] - 1, a[2] - 1, a[3] - 1, a[4] - 1);
			if (i == 50 && (j == 1 || j == 16))
				printf("%d, %d, %d, %d", a[1] - 1, a[2] - 1, a[3] - 1, a[4] - 1);
				*/
		}
	}
	fclose(f);

	f = fopen(outFileName.c_str(), "w");
	fprintf(f, "%d %d\n", client_n, obj_n);
	for (i = 1; i <= client_n; i++){
		fprintf(f, "%d ", vote_n[i]);
		for (j = 1; j <= vote_n[i]; j++)
			fprintf(f, "%d %d ", vote[i][j][0], vote[i][j][1]);
		fprintf(f, "\n");
	}
	fprintf(f, "%d\n", task_n);
	for (i = 1; i <= task_n; i++)
		fprintf(f, "%d %d\n", task_list[i][0], task_list[i][1]);
	fclose(f);
	return 0;
}
#endif