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

    int cnt = 0, counter = 0, success= 0, total = 0;
    
    int A_bar[N][M - W] = {0}, X_w[N][N] = {0}, upX_w[N][N] = {0}, upXt_w[N][N] = {0}, count[sizew] = {0}, X_id[N][N] = {0};

    int TD[M - W][W] = {0}, Z[N][M] = {0}, update_Z[N][M] = {0}, update_TD[M - W][W] = {0};

    uint8_t nonce_w[8], nonce_r[8], nonce_id[8], upnonce_r[8], upnonce_w[8], r[16], update_r[16];
   
    char DB[sizew][sizeid];

    // Allocate xbind[sizewid][N][M]
    int ***xbind = malloc(sizewid * sizeof(int **));
    if (xbind == NULL) {
        perror("malloc failed for xbind");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sizewid; i++) {
        xbind[i] = malloc(N * sizeof(int *));
        if (xbind[i] == NULL) {
            perror("malloc failed for xbind[i]");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < N; j++) {
            xbind[i][j] = calloc(M, sizeof(int));  // zero-initialized
            if (xbind[i][j] == NULL) {
                perror("calloc failed for xbind[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }

    // Allocate xtag[sizewid][N][N]
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

    // Allocate xtag_prime[sizeid][N][N]
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
    
    // int xtokent[sizeid][M][N];
    int ***xtokent = malloc(sizeid * sizeof(int **));
    if (xtokent == NULL) {
        perror("malloc failed for xtokent");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < sizeid; i++) {
        xtokent[i] = malloc(M * sizeof(int *));
        if (xtokent[i] == NULL) {
            perror("malloc failed for xtokent[i]");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < M; j++) {
            xtokent[i][j] = calloc(N, sizeof(int));  // zero-initialized
            if (xtokent[i][j] == NULL) {
                perror("calloc failed for xtokent[i][j]");
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
    
    tv_start = GetTimeStamp(); // Calculate time
    time_in_micros_start_c = 1000000 * tv_start.tv_sec + tv_start.tv_usec; // Store time in microseconds
            
    // SetUp
    generate_Random_matrix(A_bar);

    for(int h = 0; h < sizew; h++){
        
        dectohex(nonce_w, h);

        generate_Random_matrix_X_w(X_w, nonce_w);
        count[h] = 0;
        for(int l = 0; l < sizeid;l++){
            if(DB[h][l] == 1){
                count[h] = count[h] + 1;
                cnt = count[h] - 1;

                dectohex(nonce_r, (sizeid * h) + cnt);

                generate_random(r,nonce_r);
                
                Gen_Trapdoor(Z, A_bar, TD, r);

                dectohex(nonce_id, l);
                generate_Random_matrix_X_id(X_id,nonce_id);
              
                for(int i = 0; i < N;i++){
                    for(int j = 0; j < M;j++){
                        for(int k = 0; k < N;k++)
                            xbind[(sizeid * h) + cnt][i][j] = (xbind[(sizeid * h) + cnt][i][j] + (X_id[i][k] * Z[k][j])) % Q;
                        if( xbind[(sizeid * h) + cnt][i][j] < 0)
                            xbind[(sizeid * h) + cnt][i][j] += Q;
                        xbind[(sizeid * h) + cnt][i][j] = rounding((double)(((double)O/Q) * xbind[(sizeid * h) + cnt][i][j])) % O;
                        
                    }
                }

                for(int i = 0; i < N;i++){
                    for(int j = 0; j < N;j++){
                        for(int k = 0; k < N;k++)
                            xtag[(sizeid * h) + cnt][i][j] = (xtag[(sizeid * h) + cnt][i][j] + (X_id[i][k] * X_w[k][j])) % Q;
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
    
    // Search
	
    tv_start = GetTimeStamp(); // Calculate time
    time_in_micros_start_c = 1000000 * tv_start.tv_sec + tv_start.tv_usec; // Store time in microseconds
    
    dectohex(upnonce_w, 1);
    generate_Random_matrix_X_w(upX_w, upnonce_w);

    for(int i = 0;i < N;i++){
        for(int j = 0; j < N;j++)
            upXt_w[i][j] = upX_w[j][i];
    }
    
    for(int l = 0;l < count[0];l++){
        dectohex(upnonce_r,l);
        generate_random(update_r,upnonce_r);
        Gen_Trapdoor(update_Z, A_bar, update_TD, update_r);
        for(int j = 0; j < N;j++){
            SampleD_MP(A_bar,update_TD,upXt_w[j],xtoken[l][j]);
        }

    }

    for(int h = 0; h < count[0];h++){
        for(int i = 0;i < M;i++){
            for(int j = 0; j < N;j++)
                xtokent[h][i][j] = xtoken[h][j][i];
        }
    
    }
   
    for(int h = 0; h < count[0];h++){
        for(int i = 0; i < N;i++){
            for(int j = 0; j < N;j++){
                for(int k = 0; k < M;k++)
                    xtag_prime[h][i][j] = (xtag_prime[h][i][j] + (xbind[h][i][k] * xtokent[h][k][j])) % Q;
                if( xtag_prime[h][i][j] < 0)
                    xtag_prime[h][i][j] += Q;
                xtag_prime[h][i][j] = rounding((double)(((double)P/O) * xtag_prime[h][i][j])) % P;
                //if ( xtag_prime )
                
            }
        }
    
    }
    
    tv_end = GetTimeStamp(); // Calculate time
    time_in_micros_end_c = 1000000 * tv_end.tv_sec + tv_end.tv_usec; // Store time in microseconds

    time_in_micros_total_c += (time_in_micros_end_c - time_in_micros_start_c);
	
    printf("\nElapsed time Search: %ld microseconds\n", time_in_micros_total_c);

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
	
    printf("\nElapsed time Search: %ld microseconds\n", time_in_micros_total_c);
	
    printf("Number of identifiers in the intersection is %d.\n",total);
    
    return 0;

}


