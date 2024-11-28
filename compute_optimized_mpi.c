#include <omp.h>
#include <x86intrin.h>

#include "compute.h"

// Computes the convolution of two matrices
int convolve(matrix_t *a_matrix, matrix_t *b_matrix, matrix_t **output_matrix) {
  // TODO: convolve matrix a and matrix b, and store the resulting matrix in
  uint32_t rows_b = b_matrix->rows;
  uint32_t cols_b = b_matrix->cols;
  uint32_t rows_a = a_matrix->rows;
  uint32_t cols_a = a_matrix->cols;

  if(rows_b > rows_a || cols_b > cols_a){
      return -1;
  }

  uint32_t total_b = rows_b * cols_b;
  uint32_t out_rows = rows_a - rows_b + 1; 
  uint32_t out_cols = cols_a - cols_b + 1;
  matrix_t* temp_out = (matrix_t*) malloc(sizeof(matrix_t));
  temp_out->rows = out_rows;
  temp_out->cols = out_cols;

  int32_t *data = (int32_t*) calloc((out_rows * out_cols), sizeof(int32_t));
  if (data == NULL){
      printf("%s","malloc failed.");
      }
  temp_out->data = data;
    

  for (int i = 0; i < total_b/2; i++){
          int32_t curr = b_matrix->data[i];
          int32_t opp = b_matrix->data[total_b - i - 1];
          b_matrix->data[i] = opp;
          b_matrix->data[total_b - i - 1] = curr;
        
      }

// begin convolve
   
   int col_diff = cols_a - cols_b;
   int end_index = (rows_a - rows_b) * cols_a;
   int data_ind = 0;


    if (cols_b >= 8){
        #pragma omp parallel for
        for (int i = 0; i <= end_index; i+=cols_a){
            for(int z = 0; z <= col_diff; z++){
                int data_ind = (i/cols_a)*(col_diff + 1) + z;
                __m256i add_rows = _mm256_setzero_si256();
                int32_t tail = 0;
                for (int j = 0; j < rows_b; j++){
                    int start_index = i + z + (j*cols_a);
                    int x = j * cols_b;
                    int b_row_end = cols_b * (j + 1);
                    for (; x + 8 <= (b_row_end/8) * 8; x+=8){
                            __m256i a_row = _mm256_loadu_si256((__m256i*) (a_matrix->data + start_index));
                            __m256i b_row = _mm256_loadu_si256((__m256i*) (b_matrix->data + x));
                            __m256i product = _mm256_mullo_epi32(a_row, b_row);
                            add_rows = _mm256_add_epi32(add_rows, product);
                            start_index += 8;
                        }
                        int32_t dot = 0;
                        for (; x < b_row_end; x++){
                            dot += a_matrix->data[start_index] * b_matrix->data[x];
                            start_index ++;
                        }
                        tail += dot; 
                }
                int tmp[8];
                _mm256_storeu_si256((__m256i*) tmp, add_rows);
                data[data_ind] = tail + tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] +tmp[5] + tmp[6] + tmp[7];
                data_ind ++;
            }
        }
    } else{
        for (int i = 0; i <= end_index; i+=cols_a){
            int start_index = i;
            for (int z = 0; z <= col_diff; z ++){
                int matrix_pt = start_index;
                int32_t dot = 0;
                for (int j = 0; j < total_b; j++){
                    dot += a_matrix->data[matrix_pt] * b_matrix->data[j];
                    if ((j+ 1) % cols_b == 0){
                        matrix_pt = matrix_pt + cols_a - cols_b + 1;
                    } else{
                        matrix_pt ++;
                    }
                }
                start_index ++;
                data[data_ind] = dot;
                data_ind++;
            }
        }
    }
  *output_matrix = temp_out;

  return 0;
}

// Executes a task
int execute_task(task_t *task) {
  matrix_t *a_matrix, *b_matrix, *output_matrix;

  char *a_matrix_path = get_a_matrix_path(task);
  if (read_matrix(a_matrix_path, &a_matrix)) {
    printf("Error reading matrix from %s\n", a_matrix_path);
    return -1;
  }
  free(a_matrix_path);

  char *b_matrix_path = get_b_matrix_path(task);
  if (read_matrix(b_matrix_path, &b_matrix)) {
    printf("Error reading matrix from %s\n", b_matrix_path);
    return -1;
  }
  free(b_matrix_path);

  if (convolve(a_matrix, b_matrix, &output_matrix)) {
    printf("convolve returned a non-zero integer\n");
    return -1;
  }

  char *output_matrix_path = get_output_matrix_path(task);
  if (write_matrix(output_matrix_path, output_matrix)) {
    printf("Error writing matrix to %s\n", output_matrix_path);
    return -1;
  }
  free(output_matrix_path);

  free(a_matrix->data);
  free(b_matrix->data);
  free(output_matrix->data);
  free(a_matrix);
  free(b_matrix);
  free(output_matrix);
  return 0;
}
