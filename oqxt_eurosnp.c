#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <openssl/aes.h>
#include <openssl/rand.h>
#include "tools.h"


int main()
{
    int nMask11, nMask12;
                        	
    nMask11 = (1 << 18) - 1;
    
    int i = 0, j = 0, h = 0, k = 0, l = 0, cnt = 0, counter = 0, success = 0, total = 0;

    uint8_t r[16] = {0}, update_r[16] = {0}, upnonce_r[8] = {0}, upnonce_w[8] = {0};
    uint8_t nonce_w[8] = {0}, nonce_r[8] = {0}, nonce_id[8] = {0};

    int A_bar[N][M-W] = {0}, X_w[N][N] = {0}, upX_w[N][N] = {0},  X_id[N][N] = {0}, Xt_id[N][N] = {0};
    
    int TD[M-W][W] = {0}, update_TD[M-W][W] = {0}, Z[N][M] = {0}, update_Z[N][M] = {0};
    
    int count[sizew] = {0},num[sizew] = {0}, * arr;

    char DB[sizew][sizeid];
    
    //Initializing Empty Database
    for(i = 0; i < sizew; i ++)
    {
    	for(j = 0; j < sizeid; j ++)
    		DB[i][j] = 0;
    }
    
    //Initializing Empty xtag
    int ***xtag = malloc(sizewid * sizeof(int **));
    if (xtag == NULL) {
        perror("malloc failed for xtag");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < sizewid; i++) {
        xtag[i] = malloc(N * sizeof(int *));
        if (xtag[i] == NULL) {
            perror("malloc failed for xtag[i]");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < N; j++) {
            xtag[i][j] = calloc(N, sizeof(int));  // zero-initialized
            if (xtag[i][j] == NULL) {
                perror("calloc failed for xtag[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }

    //Initializing Empty xtag'
    int ***xtag_prime = malloc(sizeid * sizeof(int **));
    if (xtag_prime == NULL) {
        perror("malloc failed for xtag_prime");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < sizeid; i++) {
        xtag_prime[i] = malloc(N * sizeof(int *));
        if (xtag_prime[i] == NULL) {
            perror("malloc failed for xtag_prime[i]");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < N; j++) {
            xtag_prime[i][j] = calloc(N, sizeof(int));  // zero-initialized
            if (xtag_prime[i][j] == NULL) {
                perror("calloc failed for xtag_prime[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }
  
    //Initializing Empty y
    int ***y = malloc(sizewid * sizeof(int **));
    if (y == NULL) {
        perror("malloc failed for y");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < sizewid; i++) {
        y[i] = malloc(N * sizeof(int *));
        if (y[i] == NULL) {
            perror("malloc failed for y[i]");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < N; j++) {
            y[i][j] = calloc(M, sizeof(int));  // zero-initialized
            if (y[i][j] == NULL) {
                perror("calloc failed for y[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }

    //Initializing Empty y^T
    int ***yt = malloc(sizewid * sizeof(int **));
    if (yt == NULL) {
        perror("malloc failed for yt");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < sizewid; i++) {
        yt[i] = malloc(M * sizeof(int *));
        if (yt[i] == NULL) {
            perror("malloc failed for yt[i]");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < M; j++) {
            yt[i][j] = calloc(N, sizeof(int));  // zero-initialized
            if (yt[i][j] == NULL) {
                perror("calloc failed for yt[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }
   
    //Initializing Empty xtoken
    int ***xtoken = malloc(sizeid * sizeof(int **));
    if (xtoken == NULL) {
        perror("malloc failed for xtoken");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < sizeid; i++) {
        xtoken[i] = malloc(N * sizeof(int *));
        if (xtoken[i] == NULL) {
            perror("malloc failed for xtoken[i]");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < N; j++) {
            xtoken[i][j] = calloc(M, sizeof(int));  // zero-initialized
            if (xtoken[i][j] == NULL) {
                perror("calloc failed for xtoken[i][j]");
                exit(EXIT_FAILURE);
            }
        }
    }

    //Taking Frequency From Client
    for(i = 0; i < sizew; i++)
    {
    	printf("\nEnter frequency of keyword w_%d: ", i);
    	scanf("%d", & num[i]);
    }
    
    srand(time(NULL));
    
    // Generate Random Database
    for(i = 0; i < sizew; i++)
    {
    	arr = (int *) malloc (num[i] * sizeof(int));
    	if (arr == NULL) 
    	{
        	printf("Memory allocation failed.\n");
        	return 1;
    	}

    	for (j = 0; j < num[i]; j ++) 
        	arr[j] = j;
    	
    	shuffle(arr, num[i]);
    	
        for(j = 0; j < num[i]; j++)
            DB[i][arr[j]] = 1;
        
        free( arr );
    }
    
    //SetUp
    
    generate_Random_matrix(A_bar);

    for(h = 0; h < sizew; h++)
    {
        dectohex(nonce_w, h);

        generate_Random_matrix_X_w(X_w, nonce_w);

        count[h] = 0;
        
        for(l = 0; l < sizeid; l++)
        {
            if(l < num[h] && DB[h][l] == 1)
            {
            	printf("\nupdating for w_%d, id_%d\n", h, l);
            	
                count[h] = count[h] + 1;
                cnt = count[h] - 1;

                dectohex(nonce_r, (sizeid * h) + cnt);
                generate_random(r,nonce_r);
                
                Gen_Trapdoor(Z, A_bar, TD, r);

                dectohex(nonce_id, l);
                generate_Random_matrix_X_id(X_id, nonce_id);

                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < N;j++)
                        Xt_id[i][j] = X_id[j][i];
                }

                for(i = 0; i < N; i++){
                    SampleD_MP(A_bar, TD, Xt_id[i], y[(sizeid * h) + cnt][i]);
                }
                
                for(i = 0; i < N; i++)
                {
                    for(j = 0; j < N; j++)
                    {
                        for(k = 0; k < N; k++)
                            xtag[(sizeid * h) + cnt][i][j] = (xtag[(sizeid * h) + cnt][i][j] + (X_w[i][k] * X_id[k][j])) % Q;
                            
                        if( xtag[(sizeid * h) + cnt][i][j] >= Q)
                            	xtag[(sizeid * h) + cnt][i][j] -= Q;
                            	
                        {	
                        	nMask12 = nMask11 & xtag[(sizeid * h) + cnt][i][j];
                        	
                        	if (nMask12 > (1 << 17))
	                        	xtag[(sizeid * h) + cnt][i][j] = ((xtag[(sizeid * h) + cnt][i][j] >> 18) + 1) % P;
	                        else
	                        	xtag[(sizeid * h) + cnt][i][j] = (xtag[(sizeid * h) + cnt][i][j] >> 18) % P;                    	 
                        }
                    }
                }
            }
        }
    }
    
    //Search
   
    dectohex(upnonce_w, 1);
    generate_Random_matrix_X_w(upX_w, upnonce_w);

    for(h = 0; h < count[0]; h++)
    {
    	dectohex(upnonce_r, h);
        generate_random(update_r, upnonce_r);
        
        Gen_Trapdoor(update_Z, A_bar, update_TD, update_r);
        
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < M; j++)
            {
                for(k = 0; k < N;k++)
                    xtoken[h][i][j] = (xtoken[h][i][j] + (upX_w[i][k] * update_Z[k][j])) % Q;
                    
                if(xtoken[h][i][j] >= Q)
                    xtoken[h][i][j] -= Q;
                
                {	
                	nMask12 = nMask11 & xtoken[h][i][j];
                        	
                        if (nMask12 > (1 << 17))
	                       	xtoken[h][i][j] = ((xtoken[h][i][j] >> 18) + 1) % P;
	                else
	                       	xtoken[h][i][j] = (xtoken[h][i][j] >> 18) % P;                    	 
                }
            }
        }
    }
   
    for(h = 0; h < count[0];h++)
    {
        for(i = 0;i < M;i++)
        {
            for(j = 0; j < N;j++)
                yt[h][i][j] = y[h][j][i];
        }
    }

    for(h = 0; h < count[0]; h++)
    {
        for(i = 0; i < N; i++)
        {
            for(j = 0; j < N; j++)
            {
                for(k = 0; k < M; k++)
                    xtag_prime[h][i][j] = (xtag_prime[h][i][j] + (xtoken[h][i][k] * yt[h][k][j])) % Q;
                
                if( xtag_prime[h][i][j] > Q)
                    xtag_prime[h][i][j] -= Q;
                
                {	
                	nMask12 = nMask11 & xtag_prime[h][i][j];
                        	
                        if (nMask12 > (1 << 17))
	                       	xtag_prime[h][i][j] = ((xtag_prime[h][i][j] >> 18) + 1) % P;
	                else
	                       	xtag_prime[h][i][j] = (xtag_prime[h][i][j] >> 18) % P;                    	 
                }        
            }
        }
    }
    
    for(i = 0; i < sizeid; i++)
    {
        if(DB[0][i] == 1)
        {
            for(j = 0; j < count[1]; j++)
            {
                success = areMatricesEqual(xtag_prime[counter], xtag[sizeid + j]);
                if(success == 1)
                {
                    total++;
                    printf("\nid_%d exists in intersection.\n",i);
                    break;
                }
            }
            
            counter++;
        }
    }
    
    printf("\nNumber of identifiers in the intersection is %d.\n",total);
    
    return 0;
}
