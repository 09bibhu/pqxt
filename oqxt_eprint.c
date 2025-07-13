#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <openssl/aes.h>
#include <openssl/rand.h>
#include "tools.h"


int main(){

    srand(time(NULL));
    
    signed long time_in_micros_start_c, time_in_micros_end_c, time_in_micros_total_c = 0;

    struct timeval tv_start, tv_end;
    
    int cnt = 0, counter = 0, success = 0, total = 0;

    // uint8_t r[sizewid][16], update_r[sizeid][16], nonce_w[sizew][8], nonce_r[sizewid][8], nonce_id[sizewid][8];
    uint8_t r[16] = {0}, update_r[16] = {0}, upnonce_r[8] = {0}, upnonce_w[8] = {0}, nonce_w[8] = {0}, nonce_id[8] = {0};
    
    // int X_w[sizew][N][N], X_id[sizewid][N][N];
    int A_bar[N][M - W] = {0}, X_w[N][N] = {0}, upX_w[N][N] = {0}, count[sizew] = {0}, X_id[N][N] = {0}, Xt_id[N][N] = {0};
    
    // int TD[sizewid][M - W][W], update_TD[sizeid][M - W][W], Z[sizewid][N][M], update_Z[sizeid][N][M];
    int TD[M - W][W] = {0}, update_TD[M - W][W] = {0}, Z[N][M] = {0}, update_Z[N][M] = {0};

    char DB[sizew][sizeid];
    
    // int xtag[sizewid][N][N] = {0};
    int ***xtag = malloc(sizewid * sizeof(int **));
    if (xtag == NULL) {
        perror("malloc failed for xtag");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sizewid; i++) {
        xtag[i] = malloc(N * sizeof(int *));
        if (xtag[i] == NULL) {
            perror("malloc failed for xtag[i]");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < N; j++) {
            xtag[i][j] = calloc(N, sizeof(int));  // zero-initialized
            if (xtag[i][j] == NULL) {
                perror("calloc failed for xtag[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }

    // int xtag_prime[sizeid][N][N] = {0};
    int ***xtag_prime = malloc(sizeid * sizeof(int **));
    if (xtag_prime == NULL) {
        perror("malloc failed for xtag_prime");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sizeid; i++) {
        xtag_prime[i] = malloc(N * sizeof(int *));
        if (xtag_prime[i] == NULL) {
            perror("malloc failed for xtag_prime[i]");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < N; j++) {
            xtag_prime[i][j] = calloc(N, sizeof(int));  // zero-initialized
            if (xtag_prime[i][j] == NULL) {
                perror("calloc failed for xtag_prime[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }
  
    // int y[sizewid][N][M];
    int ***y = malloc(sizewid * sizeof(int **));
    if (y == NULL) {
        perror("malloc failed for y");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sizewid; i++) {
        y[i] = malloc(N * sizeof(int *));
        if (y[i] == NULL) {
            perror("malloc failed for y[i]");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < N; j++) {
            y[i][j] = calloc(M, sizeof(int));  // zero-initialized
            if (y[i][j] == NULL) {
                perror("calloc failed for y[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }

    // int yt[sizeid][M][N];
    int ***yt = malloc(sizewid * sizeof(int **));
    if (yt == NULL) {
        perror("malloc failed for yt");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sizewid; i++) {
        yt[i] = malloc(M * sizeof(int *));
        if (yt[i] == NULL) {
            perror("malloc failed for yt[i]");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < M; j++) {
            yt[i][j] = calloc(N, sizeof(int));  // zero-initialized
            if (yt[i][j] == NULL) {
                perror("calloc failed for yt[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }
   
    // int xtoken[sizeid][N][M];
    int ***xtoken = malloc(sizeid * sizeof(int **));
    if (xtoken == NULL) {
        perror("malloc failed for xtoken");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sizeid; i++) {
        xtoken[i] = malloc(N * sizeof(int *));
        if (xtoken[i] == NULL) {
            perror("malloc failed for xtoken[i]");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < N; j++) {
            xtoken[i][j] = calloc(M, sizeof(int));  // zero-initialized
            if (xtoken[i][j] == NULL) {
                perror("calloc failed for xtoken[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    // Generate Random Database
    for(int i = 0; i < sizew; i++){
        for(int j = 0; j < sizeid; j++){
            //DB[i][j] = rand() % 2;
            if(i % 4 == 0) {
            	if(j % 50 == 0)
            		DB[i][j] = 1;
            	else 
            		DB[i][j] = 0;
            }
            else if(1) {
            	if(j % 3 != 1)
            		DB[i][j] = 1;
            	else 
            		DB[i][j] = 0;
            }
            /*else if(i % 4 == 2) {
            	if(j % 4 != 2)
            		DB[i][j] = 1;
            	else 
            		DB[i][j] = 0;
            }
            else {
            	if(j % 4 != 3)
            		DB[i][j] = 1;
            	else 
            		DB[i][j] = 0;
            }*/
            //printf("%d ",DB[i][j]);
        }
        //printf("\n");
    }
    
    // SetUp
    
    generate_Random_matrix(A_bar);
    
    tv_start = GetTimeStamp(); // Calculate time
    time_in_micros_start_c = 1000000 * tv_start.tv_sec + tv_start.tv_usec; // Store time in microseconds

    for(int h = 0; h < sizew; h++){
        dectohex(nonce_w, h);

        generate_Random_matrix_X_w(X_w, nonce_w);
        generate_random(r, nonce_w);

        //dectohex(nonce_r, h);
        
        Gen_Trapdoor(Z, A_bar, TD, r);
        
        count[h] = 0;
        
        for(int l = 0; l < sizeid; l++){
            if(DB[h][l] == 1){
                count[h] = count[h] + 1;
                cnt = count[h] - 1;

                dectohex(nonce_id, l);
                generate_Random_matrix_X_id(X_id, nonce_id);

                for(int i = 0; i < N; i++){
                    for(int j = 0; j < N; j++)
                        Xt_id[i][j] = X_id[j][i];
                }

                for(int i = 0; i < N; i++){
                    SampleD_MP(A_bar, TD, Xt_id[i], y[(sizeid * h) + cnt][i]);
                }
                
                for(int i = 0; i < N; i++){
                    for(int j = 0; j < N; j++){
                        for(int k = 0; k < N;k++)
                            xtag[(sizeid * h) + cnt][i][j] = (xtag[(sizeid * h) + cnt][i][j] + (X_w[i][k] * X_id[k][j])) % Q;
                        if( xtag[(sizeid * h) + cnt][i][j] < 0)
                            xtag[(sizeid * h) + cnt][i][j] += Q;
                        xtag[(sizeid * h) + cnt][i][j] = rounding((double)(((double)P/Q) * xtag[(sizeid * h) + cnt][i][j])) % P; 
                    }
                }
            }
        }
    }

    tv_end = GetTimeStamp(); // Calculate time
    time_in_micros_end_c = 1000000 * tv_end.tv_sec + tv_end.tv_usec; // Store time in microseconds

    time_in_micros_total_c += (time_in_micros_end_c - time_in_micros_start_c);
	
    printf("\nElapsed time Setup: %ld microseconds\n", time_in_micros_total_c);
    
    //Search
    
    tv_start = GetTimeStamp(); // Calculate time
    time_in_micros_start_c = 1000000 * tv_start.tv_sec + tv_start.tv_usec; // Store time in microseconds
    
    dectohex(upnonce_r, 0);
    generate_random(update_r, upnonce_r);
    Gen_Trapdoor(update_Z, A_bar, update_TD, update_r);


    dectohex(upnonce_w, 1);
    generate_Random_matrix_X_w(upX_w, upnonce_w);

    for(int h = 0; h < count[0]; h++){
        for(int i = 0; i < N; i++){
            for(int j = 0; j < M; j++){
                for(int k = 0; k < N; k++)
                    xtoken[h][i][j] = (xtoken[h][i][j] + (upX_w[i][k] * update_Z[k][j])) % Q;
                if( xtoken[h][i][j] < 0)
                    xtoken[h][i][j] += Q;
                xtoken[h][i][j] = rounding((double)(((double)O/Q) * xtoken[h][i][j])) % O;
                
            }
        }
    }

    for(int h = 0; h < count[0]; h++){
        for(int i = 0;i < M; i++){
            for(int j = 0; j < N; j++)
                yt[h][i][j] = y[h][j][i];
        }
    
    }

    for(int h = 0; h < count[0]; h++){
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                for(int k = 0; k < M; k++)
                    xtag_prime[h][i][j] = (xtag_prime[h][i][j] + (xtoken[h][i][k] * yt[h][k][j])) % Q;
                if( xtag_prime[h][i][j] < 0)
                    xtag_prime[h][i][j] += Q;
                xtag_prime[h][i][j] = rounding((double)(((double)P/O) * xtag_prime[h][i][j])) % P;       
            }
        }
    }
    
    tv_end = GetTimeStamp(); // Calculate time
    time_in_micros_end_c = 1000000 * tv_end.tv_sec + tv_end.tv_usec; // Store time in microseconds

    time_in_micros_total_c += (time_in_micros_end_c - time_in_micros_start_c);
	
    printf("\nElapsed time Search 1: %ld microseconds\n", time_in_micros_total_c);

    for(int i = 0; i < sizeid; i++){
        if(DB[0][i] == 1){
            for(int j = 0; j < count[1]; j++){
                success = areMatricesEqual(xtag_prime[counter], xtag[sizeid + j]);
                if(success == 1){
                    total++;
                    //printf("id_%d exists in intersection.\n",i + 1);
                    break;
                }
            }
            counter++;
        }
    }
    
    tv_end = GetTimeStamp(); // Calculate time
    time_in_micros_end_c = 1000000 * tv_end.tv_sec + tv_end.tv_usec; // Store time in microseconds

    time_in_micros_total_c += (time_in_micros_end_c - time_in_micros_start_c);
	
    printf("\nElapsed time Search 2: %ld microseconds\n", time_in_micros_total_c);

    printf("Number of identifiers in the intersection is %d.\n",total);
    
    return 0;
    
}

