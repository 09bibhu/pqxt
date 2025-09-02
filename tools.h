#ifndef TOOLS_H

#include <sys/time.h>

#define M 110
#define W 100
#define N 5

#define Q (1<<20)
#define O (1<<17)
#define P (1<<2)
#define K 20
#define L 2
#define c 0.0

#define S 2.0
#define ROWS 2
#define COLS 2

#define sizew 5
#define sizeid 100000
#define sizewid 500000

#define AES_BLOCK_SIZE 16

#define ADD_COL(row, mat, col_i, col_j, scalar_mult)                \
    do                                                              \
    {                                                               \
        for (int l = 0; l < (row); l++)                             \
        {                                                           \
            (mat)[l][(col_i)] -= (scalar_mult) * (mat)[l][(col_j)]; \
        }                                                           \
    } while (0);

#define INNER_PROD(row, basis_mat, GSO_mat, col_i, col_j, inner_prod)              \
    do                                                                             \
    {                                                                              \
        double temp;                                                               \
        (inner_prod) = 0.0;                                                        \
        for (int l = 0; l < (row); l++)                                            \
        {                                                                          \
            temp = (double)(basis_mat)[l][(col_i)];                                \
            (inner_prod) += temp * (GSO_mat)[l][(col_j)];                          \
        }                                                                          \
    } while (0)
		
//----------------------------------------------------------------------------------------------//

//Global Variable
int p1[M];
	
int G[N][W], x1[100][P][L] = {0};

//Time Structure
struct timeval GetTimeStamp()
{
   struct timeval tv;
   
   gettimeofday(&tv,NULL);
   
   return tv;
}

//AES Structure and Variable
typedef struct {
    AES_KEY aes_key;                 
    uint8_t counter[AES_BLOCK_SIZE]; 
    uint8_t stream_block[AES_BLOCK_SIZE]; 
    int stream_index;                
} AES_CTR_PRNG;

AES_CTR_PRNG prng1;

//-----------------------------------Function Declaration------------------------------//

//Initializing Global Variable
void global();

//General Functionalities
void dectohex(uint8_t[8], int);

int rounding(double);

int areMatricesEqual(int **, int **); 

void swap(int *, int *);

void shuffle(int *, int); 

//AES Functionalities
void aes_ctr_prng_init(AES_CTR_PRNG *, uint8_t[16], uint8_t[8]); 

void aes_ctr_generate_block(AES_CTR_PRNG *); 

uint32_t aes_ctr_random_int(AES_CTR_PRNG *); 

void aes_ctr_prng_get_128bit(AES_CTR_PRNG *, uint8_t[16]); 

double aes_ctr_random_real(AES_CTR_PRNG *); 

int aes_ctr_random_int_range(AES_CTR_PRNG *, int); 

//Vector Operations
double norm(int[M]); 

double vector_norm(double[ROWS]);

int product(int[M], int[M]);

double dot_product(double[ROWS], double[ROWS]); 

int init_dot_product(int[ROWS], int[ROWS]);  

void compute_GSO(double[ROWS][COLS], double[ROWS][COLS]);

//Random Vector and Matrix 
void generate_random(uint8_t[16], uint8_t[8]); 

void generate_Random_matrix(int[N][M - W]);

void generate_small_matrix(int[M - W][W], uint8_t[16]);

void generate_Random_matrix_X_w(int[N][N], uint8_t[8]);

void generate_Random_matrix_X_id(int[N][N], uint8_t[8]);

//Trapdoor Matrix
//void generate_gadget_matrix(int[N][W]);

void Gen_Trapdoor(int[N][M], int[N][M - W], int[M - W][W],uint8_t[16]);

//Sampling 
int SampleZ(double, AES_CTR_PRNG *);

void SampleD(double[ROWS][COLS],int[ROWS], AES_CTR_PRNG *);

void SampleD_SL(int[N], int[N * K], AES_CTR_PRNG *);

void SampleD_MP(int[N][M - W], int[M - W][W], int[N], int[M]);

//-----------------------------------Function Definition------------------------------//

void global()
{
    int i, j, p;
    int g[ROWS], v[L], cnt[P] = {0};
    
    double input[ROWS][COLS];
    
    uint8_t key[16] = 
    {
        0x1b, 0x7f, 0x15, 0x16,
        0x08, 0xae, 0xd2, 0xa6,
        0x2b, 0xf7, 0x15, 0x88,
        0x09, 0xcf, 0x4f, 0x3c
    };

    uint8_t nonce[8] = 
    {
    	0x00, 0x00, 0x00, 0x00, 
    	0x00, 0x00, 0x00, 0x00
    };
    	
    aes_ctr_prng_init(&prng1, key, nonce);
    
    for (i = 0; i < N; i++) 
    {
        for(j = 0; j < K; j++)
        	G[i][(i * K) + j] = pow(2, j);
            
        for(j = K; j < W; j++) 
            	G[i][((i * K) + j)% W] = 0;
    }
    
    for(i = 0; i < M; i++)
        p1[i] = SampleZ(6.0,&prng1);
   
    for(i = 0; i < ROWS; i++)
    {
	for(j = 0;j < COLS;j++)
	{
		if (i==j)
			input[i][j] = 1.0;
		else
			input[i][j] = 0.0;
	}
    }

    for(i = 0; i < L; i++)
        g[i] = (1 << i);
    
    for(i = 0; i < 110; i++)
    {
        for(p = 0; p < L; p++)
            v[p] = 0;
            
        SampleD(input,v,&prng1);
    
        j = (init_dot_product(g, v) % P);
        
        if(j < 0)
            j = j + P;
   
        for(p = 0; p < L; p++)
            x1[cnt[j]][j][p] = v[p];
            
        cnt[j]++;
    }
}
    
void dectohex(uint8_t arr[8], int i) 
{
    int j;
    
    for (j = 0; j < 8; j++) 
        arr[7 - j] = (i >> (j * 8)) & 0xFF;

}

int areMatricesEqual(int **mat1, int **mat2) 
{
    int i, j;
    
    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < N; j++) 
        {
        	if (mat1[i][j] != mat2[i][j])
                	return 0; 
        }
    }
    
    return 1; 
}

void swap(int * a, int * b) 
{
    int temp;
    
    temp  = * a;
    * a = * b;
    * b = temp;
}

void shuffle(int *arr, int n)
{
    int i, j; 
    
    for (i = n - 1; i > 0; i --) 
    {
        j = rand() % (i + 1);  
        
        swap(&arr[i], &arr[j]);
    }
}

int rounding(double x)
{
    int i;
    
    i = floor(x);
    
    if((x - i) <= 0.5)
        return i;
    else 
        return (i + 1);
}

int product(int v1[M], int v2[M]) 
{
    int i = 0, sum = 0;
    
    for (i = 0; i < M; i++) 
        sum += v1[i] * v2[i];

    return sum;
}

double dot_product(double v1[ROWS], double v2[ROWS]) 
{
    int i;
    
    double sum = 0.0;
    
    for (i = 0; i < ROWS; i++) 
        sum += v1[i] * v2[i];
    
    return sum;
}

int init_dot_product(int v1[ROWS], int v2[ROWS]) 
{
    int i = 0, sum = 0;
    
    for (i = 0; i < ROWS; i++)
        sum += v1[i] * v2[i];
    
    return sum;
}

double norm(int vector[M]) 
{
    double val;
    
    val = sqrt(product(vector, vector));
    
    return val;
}

double vector_norm(double vector[ROWS]) 
{
    double val;
    
    val = sqrt(dot_product(vector, vector));
    
    return val;
}

void aes_ctr_prng_init(AES_CTR_PRNG *prng, uint8_t key[16], uint8_t nonce[8]) 
{   
    memcpy(prng->counter, nonce, 8); 
    memset(prng->counter + 8, 0, 8); 
    
    AES_set_encrypt_key(key, 128, &prng->aes_key); 
    
    prng->stream_index = AES_BLOCK_SIZE;  
}
    
void aes_ctr_generate_block(AES_CTR_PRNG *prng) 
{
    int i;
    
    AES_encrypt(prng->counter, prng->stream_block, &prng->aes_key);
        
    for (i = 15; i >= 8; i--) 
    {
    	if (++prng->counter[i]) 
        	break;  
    }
        
        prng->stream_index = 0;  
}
    
uint32_t aes_ctr_random_int(AES_CTR_PRNG *prng) 
{
    uint32_t rand_value;
    
    if (prng->stream_index + sizeof(uint32_t) > AES_BLOCK_SIZE) 
    	aes_ctr_generate_block(prng);
    
    memcpy(&rand_value, prng->stream_block + prng->stream_index, sizeof(uint32_t));
    prng->stream_index += sizeof(uint32_t);
        
    return rand_value;
}
    
void aes_ctr_prng_get_128bit(AES_CTR_PRNG *prng, uint8_t out[16]) 
{
    int remaining;
    
    if (prng->stream_index == AES_BLOCK_SIZE) 
        aes_ctr_generate_block(prng);
    
    remaining = AES_BLOCK_SIZE - prng->stream_index;
    
    if (remaining >= 16) 
    {
        memcpy(out, prng->stream_block + prng->stream_index, 16);
        prng->stream_index += 16;
    } 
    else 
    {
        memcpy(out, prng->stream_block + prng->stream_index, remaining);
        aes_ctr_generate_block(prng);
            
        memcpy(out + remaining, prng->stream_block, 16 - remaining);
        prng->stream_index = 16 - remaining;
    }
}

double aes_ctr_random_real(AES_CTR_PRNG *prng) 
{
    return ((double)aes_ctr_random_int(prng) / (double)UINT32_MAX);
}

int aes_ctr_random_int_range(AES_CTR_PRNG *prng, int s) 
{
    int min, max;
    
    uint32_t range, rand_val;
    
    min = c - s;
    max = c + s;
    
    range = (uint32_t)(max - min + 1);
    
    rand_val = aes_ctr_random_int(prng) % range;
    
    return min + rand_val;
}

void compute_GSO(double inp_basis[ROWS][COLS], double GSO_basis[ROWS][COLS])
{
    int i, j, k;
    
    double inner_prod, norm;
    double mu_mat[ROWS][COLS] = {0};
    
    for(i= 0; i < ROWS; i++)                            
        GSO_basis[i][0] = inp_basis[i][0];
        
    for(j = 1; j < COLS; j++)
    {
        for(i = 0; i < ROWS; i++)                      
            GSO_basis[i][j] = inp_basis[i][j];
            
        for(k = j - 1; k > -1; k--)             
        {
            INNER_PROD(ROWS, inp_basis, GSO_basis, j, k, inner_prod);
            INNER_PROD(ROWS, GSO_basis, GSO_basis, k, k, norm);
            
            mu_mat[j][k] = inner_prod / norm;
            
            ADD_COL(ROWS, GSO_basis, j, k, mu_mat[j][k]);
        }
    }
}

void generate_Random_matrix_X_w(int matrix[N][N], uint8_t nonce[8]) 
{
    int i, j; 
    
    uint8_t key_X[16] = {
        0x1b, 0x7f, 0x15, 0x16,
        0x08, 0xae, 0xd2, 0xa6,
        0x2b, 0xf7, 0x15, 0x88,
        0x09, 0xcf, 0x4f, 0x3c};
    
    AES_CTR_PRNG prng;
    
    aes_ctr_prng_init(&prng, key_X, nonce);

    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < N; j++) 
        {
        	matrix[i][j] = (int) aes_ctr_random_int(&prng) % Q;
		
		if(matrix[i][j] < 0)
            		matrix[i][j] = matrix[i][j] + Q;
        }
    }
}

void generate_Random_matrix_X_id(int matrix[N][N], uint8_t nonce[8]) 
{
    int i, j;
    
    uint8_t key_I[16] = {
        0x2b, 0x7e, 0x15, 0x16,
        0x28, 0xae, 0xd2, 0xa6,
        0xab, 0xf7, 0x15, 0x88,
        0x09, 0xcf, 0x4f, 0x3c};
    
    AES_CTR_PRNG prng;
    
    aes_ctr_prng_init(&prng, key_I, nonce);

    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < N; j++) 
        {
        	matrix[i][j] = (int) aes_ctr_random_int(&prng) % Q;
	
		if(matrix[i][j] < 0)
            		matrix[i][j] = matrix[i][j] + Q;
        }
    }
}

void generate_random(uint8_t output[16],uint8_t nonce[8]) 
{
    uint8_t key[16] = {
        0x00, 0x01, 0x02, 0x03,
        0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b,
        0x0c, 0x0d, 0x0e, 0x0f
    };
    
    AES_CTR_PRNG prng;
    
    aes_ctr_prng_init(&prng, key, nonce);
    aes_ctr_prng_get_128bit(&prng, output);
}

void generate_Random_matrix(int matrix[N][M - W]) 
{
    int i, j;
    
    uint8_t key[16] = {
        0x00, 0x01, 0x02, 0x03,
        0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b,
        0x0c, 0x0d, 0x0e, 0x0f
    };
    
    uint8_t nonce[8] = {
        0xaa, 0xbb, 0xcc, 0xdd,
        0xee, 0xff, 0x11, 0x22
    };
    
    AES_CTR_PRNG prng;
    
    aes_ctr_prng_init(&prng, key, nonce);
    
    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < M - W; j++) 
        {
        	matrix[i][j] = (int) aes_ctr_random_int(&prng) % Q;
        	
		if(matrix[i][j] < 0)
	            	matrix[i][j] = matrix[i][j] + Q;
        }
    }
}

void generate_small_matrix(int matrix[M - W][W],uint8_t key[16])
{
    int i, j;
    
    uint8_t nonce[8] = {
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00
    };
    
    AES_CTR_PRNG prng;
    
    aes_ctr_prng_init(&prng, key, nonce);
    
    for(i = 0;i < M - W;i++)
    {
	for(j = 0;j< W;j++)
		matrix[i][j] = SampleZ(3.0,&prng);
    }
}

/*void generate_gadget_matrix(int G[N][W]) 
{ 
    int i, j;
    
    for (i = 0; i < N; i++) 
    {
        for(j = 0; j < K; j++)
        	G[i][(i * K) + j] = pow(2, j);
            
        for(j = K; j < W; j++) 
            	G[i][((i * K) + j)% W] = 0;
    }
}
*/
void Gen_Trapdoor(int A[N][M], int A_prime[N][M - W], int R[M - W][W],uint8_t key[16])
{
    int i, j, k, sum;
    //int G[N][W];
    
    generate_small_matrix(R,key);
    
    //generate_gadget_matrix(G);

    for (i = 0; i < N; i++) 
    {
        for (j = 0; j < M - W; j++) 
            A[i][j] = A_prime[i][j];
    }

    for (i = 0; i < N; i++) 
    {
        for (j = M - W; j < M ; j++) 
        {
        	sum = 0;
        	
            	for (k = 0; k < M - W; k++)
                	sum += A_prime[i][k] * R[k][j - M + W];

            	A[i][j] = ( G[i][j - M + W] - sum) % Q;
            
            	if(A[i][j] < 0)
                	A[i][j] = A[i][j] + Q;
        }
    }
}

int SampleZ(double s, AES_CTR_PRNG *prng) 
{
    int x;
    
    double u_2, accept_prob;

    do 
    {
        x = aes_ctr_random_int_range(prng,s);
            
        u_2 = aes_ctr_random_real(prng) ;
        
        accept_prob = exp(-((M_PI * x * x)/(s * s)));    
        
    } while(u_2 > accept_prob);

    return x;
}

void SampleD(double input[ROWS][COLS],int v[ROWS], AES_CTR_PRNG *prng)
{
    int i, j, z;
    
    double s_prime, norm;
    double output[ROWS][COLS] = {0};
    
    compute_GSO(input, output);
    
    for (i = COLS - 1; i >= 0; i--)
    {
        norm = vector_norm(output[i]);
         
        s_prime = S/norm;
        
        z = SampleZ(s_prime, prng);
        
        for(j = 0; j < ROWS; j++)   
            v[j] = v[j] + (z * input[j][i]);
    }

}

void SampleD_SL(int u[N], int result[N * K], AES_CTR_PRNG *prng)
{
    int i, j, k, X;
    int u_prime[N], output[K], count[P] = {0};
    
    for(i = 0; i < N; i++)
    {   
	for(j = 0; j < (K/L); j++)
	{
            X = 0;
            
            u_prime[i] = (u[i] % P);
            
            if(u_prime[i] < 0)
                u_prime[i] = u_prime[i] + P;
            
            for(k = 0; k < L; k++)
            {
		output[(j * L) + k] = x1[count[u_prime[i]]][u_prime[i]][k];
		
                X = X + (output[(j * L) + k] * (1 << k));
            }
            
            count[u_prime[i]]++;
            
            u[i] = (u[i] - X)/P;
	}

        for(j = 0; j < K; j++)
            result[(i * K) + j] = output[j];	    
    }
}

void SampleD_MP(int A_prime[N][M - W], int R[M - W][W], int u[N], int x[M]) 
{
    int i, j, k, sum; 
    
    int y[M - W] = {0}, W_prime[N] = {0}, w[N] = {0}, v[N] = {0};
    int z[N * K], result[M], /*G[N][W],*/ ext_R[M][W];
   
    for(i = 0; i < M - W; i++)
    {
        sum = 0;
        
        for(j = 0; j < W; j++)
        	sum += R[i][j] * p1[j + M - W];
      
        y[i] = p1[i] - sum;
    }

    for(i = 0; i < N; i++)
    {
        sum = 0;
        
        for(j = 0; j < M - W; j++)
            	sum += A_prime[i][j] * y[j];

        W_prime[i] = sum % Q;
        
        if(W_prime[i] < 0)
            	W_prime[i] = W_prime[i] + Q;
    }

    //generate_gadget_matrix(G);

    for(i = 0; i < N; i++)
    {
        sum = 0;
        
        for(j = 0; j < W; j++)
        	sum += G[i][j] * p1[j + M - W];
        
        w[i] = sum % Q;
        
        if(w[i] < 0)
        	w[i] = w[i] + Q;
    }

     for(i = 0; i < N; i++)
     { 
         v[i] = (u[i] - W_prime[i] - w[i]) % Q;
       
         if(v[i] < 0)
	 	v[i] = v[i] + Q;
    }
	
    SampleD_SL(v,z,&prng1);

    for(i = 0; i < M - W; i++)
    {
	for(j = 0; j < W; j++)		
		ext_R[i][j] = R[i][j];
    }

    for(i = M - W; i < M; i++)
    {
	for(j = 0; j < W; j++)
	{
		k = j + M - W;
		if(i == k)
			ext_R[i][j] = 1;
		else
			ext_R[i][j] = 0;
	}		
    }

    for(i = 0; i < M; i++)
    {
        sum = 0;
        
        for(j = 0; j < W; j++)
        	sum = sum + (ext_R[i][j] * z[j]);

        result[i] = sum;
    }
  
    for(i = 0; i < M; i++)
        x[i] = p1[i] + result[i];
}

#endif

