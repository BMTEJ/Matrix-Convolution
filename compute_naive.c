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

  int32_t *data = (int32_t*) malloc((out_rows * out_cols) * sizeof(int32_t));
  if (data == NULL){
      printf("%s","malloc failed.");
      }
  temp_out->data = data;

  /*vertical flip on b_matrix*/
  int rev_i = total_b - 1;
  int32_t curr_int = 0;
  int32_t rev_int = 0;
  for (int i = 0; i < total_b; i++){
      if (rev_i <= i) break;
        curr_int = b_matrix->data[i];
        rev_int = b_matrix->data[rev_i];
        b_matrix->data[rev_i] = curr_int;
        b_matrix->data[i] = rev_int;
        rev_i --;
      }

/* begin convolve*/
   
   int col_diff = cols_a - cols_b;
   int end_index = (rows_a - rows_b) * cols_a;
   int start_index = 0;
   int data_ind = 0;
   int matrix_pt = 0;
   int32_t dot = 0;

   /*Outer loop hits start of matrix_a rows until end_index*/
   for (int i = 0; i <= end_index; i+=cols_a){
       start_index = i;
        
       /*Second loop goes down each row of a*/
    for (int z = 0; z <= col_diff; z ++){
          matrix_pt = start_index;
          dot = 0;
          
          /*Navigate through b*/
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
