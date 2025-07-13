
#include <sys/time.h>

#define M 110
#define W 100
#define N 5
#define Q (1<<20)
#define O (1<<17)
#define P (1<<2)
#define S 2.0
#define ROWS 2
#define COLS 2
#define K 20
#define L 2
#define c 0.0
#define sizew 2
#define sizeid 50000
#define sizewid 100000
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
	
	struct timeval GetTimeStamp()
	{
		struct timeval tv;
		gettimeofday(&tv,NULL);
    	return tv;
	}

    typedef struct {
        AES_KEY aes_key;                 // AES key schedule
        uint8_t counter[AES_BLOCK_SIZE]; // 8-byte nonce + 8-byte counter
        uint8_t stream_block[AES_BLOCK_SIZE]; // Output of AES encryption
        int stream_index;                // Position in stream block
    } AES_CTR_PRNG;
    
    // Initialize the AES-CTR PRNG with a random key and nonce
    void aes_ctr_prng_init(AES_CTR_PRNG *prng, uint8_t key[16], uint8_t nonce[8]) {
    
        memcpy(prng->counter, nonce, 8); // Set the nonce part always to 0
        memset(prng->counter + 8, 0, 8);  // Set the counter part to 0
    
        AES_set_encrypt_key(key, 128, &prng->aes_key); // Set AES encryption key
        prng->stream_index = AES_BLOCK_SIZE;  // Start at the end to trigger a new block
    }
    
    // Generate 16 bytes of pseudorandom output by encrypting the counter
    void aes_ctr_generate_block(AES_CTR_PRNG *prng) {
        AES_encrypt(prng->counter, prng->stream_block, &prng->aes_key);
        
        // Increment counter (only last 8 bytes)
        for (int i = 15; i >= 8; i--) {
            if (++prng->counter[i]) break;  // Stop if no carry
        }
        
        prng->stream_index = 0;  // Reset stream index
    }
    
    // Generate a random 32-bit integer
    uint32_t aes_ctr_random_int(AES_CTR_PRNG *prng) {
        if (prng->stream_index + sizeof(uint32_t) > AES_BLOCK_SIZE) {
            aes_ctr_generate_block(prng);
        }
        
        uint32_t rand_value;
        memcpy(&rand_value, prng->stream_block + prng->stream_index, sizeof(uint32_t));
        prng->stream_index += sizeof(uint32_t);
        
        return rand_value;
    }
    
    
    // Get exactly 16 bytes (128 bits) of pseudorandom output
    void aes_ctr_prng_get_128bit(AES_CTR_PRNG *prng, uint8_t out[16]) {
        if (prng->stream_index == AES_BLOCK_SIZE) {
            aes_ctr_generate_block(prng);
        }
    
        int remaining = AES_BLOCK_SIZE - prng->stream_index;
        if (remaining >= 16) {
            memcpy(out, prng->stream_block + prng->stream_index, 16);
            prng->stream_index += 16;
        } else {
            memcpy(out, prng->stream_block + prng->stream_index, remaining);
            aes_ctr_generate_block(prng);
            memcpy(out + remaining, prng->stream_block, 16 - remaining);
            prng->stream_index = 16 - remaining;
        }
    }
// Generate a random real number in [0,1]
double aes_ctr_random_real(AES_CTR_PRNG *prng) {
    return ((double)aes_ctr_random_int(prng) / (double)UINT32_MAX);
}

int aes_ctr_random_int_range(AES_CTR_PRNG *prng, int s) {
    int min = c - s;
    int max = c + s;
    uint32_t range = (uint32_t)(max - min + 1);

    // Generate a random number in [0, range - 1], then shift to [min, max]
    uint32_t rand_val = aes_ctr_random_int(prng) % range;
    return min + rand_val;
}


int SampleZ(double s, AES_CTR_PRNG *prng) {
    int x;
    double u_2, accept_prob;

    do {
        x = aes_ctr_random_int_range(prng,s);    
        u_2 = aes_ctr_random_real(prng) ;
        accept_prob = exp(-((M_PI * x * x)/(s * s)));    
    } while(u_2 > accept_prob);

    return x;
}

int rounding(double x){
    int i;
    int y;
    i = floor(x);
    if((x - i) <= 0.5)
        y = i;
    else 
        y = i + 1;
    return y;
}

double dot_product(double v1[ROWS], double v2[ROWS]) {
    double sum = 0.0;
    for (int i = 0; i < ROWS; i++) {
        sum += v1[i] * v2[i];
    }
    return sum;
}


int init_dot_product(int v1[ROWS], int v2[ROWS]) {
    int sum = 0;
    for (int i = 0; i < ROWS; i++) {
        sum += v1[i] * v2[i];
    }
    return sum;
}


double vector_norm(double vector[ROWS]) {
    return sqrt(dot_product(vector, vector));
}

int product(int v1[M], int v2[M]) {
    int sum = 0;
    for (int i = 0; i < M; i++) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

double norm(int vector[M]) {
    return sqrt(product(vector, vector));
}


void compute_GSO(double inp_basis[ROWS][COLS], double GSO_basis[ROWS][COLS])
{
    double inner_prod, norm;
	double mu_mat[ROWS][COLS] = {0};
    
    for(int i= 0; i < ROWS; i++)                            // Setting b_0^* = b_0
        GSO_basis[i][0] = inp_basis[i][0];
        
    for(int j = 1; j < COLS; j++)
    {
        for(int i = 0; i < ROWS; i++)                       // Setting b_j^* = b_j
            GSO_basis[i][j] = inp_basis[i][j];
            
        for(int k = j - 1; k > -1; k--)                    // Update b_j^* = b_j^* - (sum(k = j - 1 down to 0) mu(j,k) . b_k^*)  
        {
            INNER_PROD(ROWS, inp_basis, GSO_basis, j, k, inner_prod);
            INNER_PROD(ROWS, GSO_basis, GSO_basis, k, k, norm);
            
            mu_mat[j][k] = inner_prod / norm;
            
            ADD_COL(ROWS, GSO_basis, j, k, mu_mat[j][k]);
           
        }
    }
}




void SampleD(double input[ROWS][COLS],int v[ROWS], AES_CTR_PRNG *prng){
    double s_prime;
    int i;
    int j;
    double norm;
    int z;
    double output[ROWS][COLS] = {0};
    compute_GSO(input, output);
    

    for (i = COLS - 1; i >= 0; i--){
        norm = vector_norm(output[i]); 
        s_prime = S/norm;
        z = SampleZ(s_prime, prng);
        for(j = 0; j < ROWS; j++){
           
            v[j] = v[j] + (z * input[j][i]);
        }
        
    }

}



void SampleD_SL(int u[N], int result[N * K], AES_CTR_PRNG *prng){
    int i, j, k, X;
    int v[L];
    int u_prime[N];
    int output[K];
    int p;
	int count[P] = {0};
    int cnt[P] = {0};
    int x[100][P][L] = {0};
    double input[ROWS][COLS];
    int g[ROWS];


    for(i = 0;i < ROWS;i++){
		for(j = 0;j < COLS;j++){
			if (i==j)
				input[i][j] = 1.0;
			else
				input[i][j] = 0.0;
		}
	}

    for(i = 0;i < L;i++){
        g[i] = pow(2,i);
    }
    

    for(i = 0;i <110; i++){
        for(p = 0;p < L;p++){
            v[p] = 0;
        }
        SampleD(input,v,prng);
    
        j = (init_dot_product(g, v) % P);
        if(j < 0){
            j = j + P;
        }
   
       for(p = 0;p < L;p++){
            x[cnt[j]][j][p] = v[p];
            
        }
        cnt[j]++;
        
    }
    


	for(i = 0;i < N;i++){
       
		for(j = 0;j < (K/L);j++){
            X = 0;
            u_prime[i] = (u[i] % P);
            if(u_prime[i] < 0)
                u_prime[i] = u_prime[i] + P;
            
            for(k = 0;k < L;k++){
		        output[(j * L) + k] = x[count[u_prime[i]]][u_prime[i]][k];
                X = X + (output[(j * L) + k] * (1 << k));
                
            }
            count[u_prime[i]]++;
            u[i] = (u[i] - X)/P;
            
            
		}

        for(j = 0;j < K;j++){
            result[(i * K) + j] = output[j];
        }
		
		    
	}
	
  
}



void generate_Random_matrix_X_w(int matrix[N][N], uint8_t nonce[8]) {

	AES_CTR_PRNG prng;
	uint8_t key_X[16] = {
        0x1b, 0x7f, 0x15, 0x16,
        0x08, 0xae, 0xd2, 0xa6,
        0x2b, 0xf7, 0x15, 0x88,
        0x09, 0xcf, 0x4f, 0x3c};
	
	aes_ctr_prng_init(&prng, key_X, nonce);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = (int) aes_ctr_random_int(&prng) % Q;
			if(matrix[i][j] < 0)
            	matrix[i][j] = matrix[i][j] + Q;
        }
    }
}


void generate_Random_matrix_X_id(int matrix[N][N], uint8_t nonce[8]) {

	AES_CTR_PRNG prng;
	uint8_t key_I[16] = {
        0x2b, 0x7e, 0x15, 0x16,
        0x28, 0xae, 0xd2, 0xa6,
        0xab, 0xf7, 0x15, 0x88,
        0x09, 0xcf, 0x4f, 0x3c};
	
	aes_ctr_prng_init(&prng, key_I, nonce);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = (int) aes_ctr_random_int(&prng) % Q;
			if(matrix[i][j] < 0)
            	matrix[i][j] = matrix[i][j] + Q;
        }
    }
}

void generate_random(uint8_t output[16],uint8_t nonce[8] ) {
    AES_CTR_PRNG prng;
    uint8_t key[16] = {
        0x00, 0x01, 0x02, 0x03,
        0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b,
        0x0c, 0x0d, 0x0e, 0x0f
    };
    

    aes_ctr_prng_init(&prng, key, nonce);
    aes_ctr_prng_get_128bit(&prng, output);

}


// Generate Random matrix A_bar
void generate_Random_matrix(int matrix[N][M - W]) {

    AES_CTR_PRNG prng;
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
    

    aes_ctr_prng_init(&prng, key, nonce);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M - W; j++) {
            matrix[i][j] = (int) aes_ctr_random_int(&prng) % Q;
			if(matrix[i][j] < 0)
            	matrix[i][j] = matrix[i][j] + Q;
        }
    }
}


// Generate small matrix R
void generate_small_matrix(int matrix[M - W][W],uint8_t key[16] ){
	int i, j;
    AES_CTR_PRNG prng;
    uint8_t nonce[8] = {
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00
    };

    aes_ctr_prng_init(&prng, key, nonce);
	for(i = 0;i < M - W;i++){
		for(j = 0;j< W;j++){
			
			matrix[i][j] = SampleZ(3.0,&prng);
           
		}

	}
}

// Generate gadget matrix G
void generate_gadget_matrix(int G[N][W]) { 
    int i;
    int j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < K; j++) {
            G[i][(i * K) + j] = pow(2, j);
        }
            
        for (j = K; j < W; j++) {
            G[i][((i * K) + j)% W] = 0;
            
        }
       
    }

}

void Gen_Trapdoor(int A[N][M], int A_prime[N][M - W], int R[M - W][W],uint8_t key[16]){
    int G[N][W];

    //generate_Random_matrix(A_prime,prng);
    generate_small_matrix(R,key);
    generate_gadget_matrix(G);

    // Copy A' to the left half of A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M - W; j++) {
            A[i][j] = A_prime[i][j];
        }
    }
    

    // Compute A = [A' | G - A'R ] mod Q

    for (int i = 0; i < N; i++) {
        for (int j = M - W; j < M ; j++) {
            int sum = 0;
            for (int k = 0; k < M - W; k++) {
                sum += A_prime[i][k] * R[k][j - M + W];
            }
            A[i][j] = ( G[i][j - M + W] - sum) % Q;
            if(A[i][j] < 0)
                A[i][j] = A[i][j] + Q;
    
        }
        
    }
}



void SampleD_MP(int A_prime[N][M - W], int R[M - W][W], /*int A[N][M],*/ int u[N], int x[M]) {
   
    int G[N][W]; 
    int p[M];
    int p_1[M - W] = {0};
    int p_2[W] = {0};
    int y[M - W] = {0};
    int W_prime[N] = {0};
    int w[N] = {0};
    int v[N] = {0};
    int i;
    int j;
    int k;
    //double input[ROWS][COLS];
    int z[N * K];
    int ext_R[M][W];
    int result[M];
    int sum;
    
    AES_CTR_PRNG prng;
	uint8_t key[16] = {
        0x1b, 0x7f, 0x15, 0x16,
        0x08, 0xae, 0xd2, 0xa6,
        0x2b, 0xf7, 0x15, 0x88,
        0x09, 0xcf, 0x4f, 0x3c};

    uint8_t nonce[8] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
	
	aes_ctr_prng_init(&prng, key, nonce);
    
    for(i = 0; i < M;i++){
        p[i] = SampleZ(6.0,&prng);
    }
   


    //Compute p_1
    for ( i = 0; i < M - W; i++) {
        p_1[i] = p[i];
    }
    

    //Compute p_2
    for ( i = 0; i < W; i++) {
        p_2[i] = p[i+ M - W];
    }
   

    //Compute y = p_1 - R * p_2
    for( i= 0; i< M - W ; i++){
         sum = 0;
        for( j = 0; j< W; j++)
            sum += R[i][j] * p_2[j];
        y[i] = p_1[i] - sum;
    }

    //Compute W' = A' * y
    for( i= 0; i< N ; i++){
         sum = 0;
        for( j = 0; j< M - W; j++){
            sum += A_prime[i][j] * y[j];
        }
        W_prime[i] = (sum % Q);
        if(W_prime[i] < 0)
            W_prime[i] = W_prime[i] + Q;
    }

    
    generate_gadget_matrix(G);

    //Compute W = G * p_2
    for( i= 0; i< N ; i++){
         sum = 0;
        for( j = 0; j< W; j++){
            sum += G[i][j] * p_2[j];
        }
        w[i] = (sum % Q);
        if(w[i] < 0)
            w[i] = w[i] + Q;
    }

    //Compute v
     for( i= 0; i< N ; i++){
        
       v[i] = (u[i] - W_prime[i] - w[i]) % Q;
	   if(v[i] < 0)
	   		v[i] = v[i] + Q;
    }
	
	SampleD_SL(v,z,&prng);

    for( i = 0; i < M - W; i++){
		for( j = 0; j < W; j++)		
			ext_R[i][j] = R[i][j];
	}

	for( i = M - W; i < M; i++){
		for( j = 0; j < W; j++){
			 k = j + M - W;
			if(i == k)
				ext_R[i][j] = 1;
			else
				ext_R[i][j] = 0;
		}		
	}

    //compute [R I]^t * z
    

    for(i = 0; i < M; i++){
        sum = 0;
        for(j = 0; j < W; j++){
            sum = sum + (ext_R[i][j] * z[j]);

        }
        result[i] = sum;
    }
    
    //compute final output
    for(i = 0; i < M;i++){
        x[i] = p[i] + result[i];
    }
      
}


void dectohex(uint8_t arr[8], int i) {

    for (int j = 0; j < 8; j++) {
        arr[7 - j] = (i >> (j * 8)) & 0xFF; // Big-endian
    }

}

int areMatricesEqualold(int mat1[N][N], int mat2[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (mat1[i][j] != mat2[i][j]) {
                return 0; 
            }
        }
    }
    return 1; 
}


int areMatricesEqual(int **mat1, int **mat2) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (mat1[i][j] != mat2[i][j]) {
                return 0; 
            }
        }
    }
    return 1; 
}
