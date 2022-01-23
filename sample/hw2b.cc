#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <sched.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <omp.h>
#include <mpi.h>


struct timespec diff(struct timespec start, struct timespec end) {
  struct timespec temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp;
}

void write_png(const char* filename, int iters, int width, int height, const int* buffer) {
    FILE* fp = fopen(filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_filter(png_ptr, 0, PNG_NO_FILTERS);
    png_write_info(png_ptr, info_ptr);
    png_set_compression_level(png_ptr, 1);
    size_t row_size = 3 * width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);
    for (int y = 0; y < height; ++y) {
        memset(row, 0, row_size);
        for (int x = 0; x < width; ++x) {
            int p = buffer[(height - 1 - y) * width + x];
            png_bytep color = row + x * 3;
            if (p != iters) {
                if (p & 16) {
                    color[0] = 240;
                    color[1] = color[2] = p % 16 * 16;
                } else {
                    color[0] = p % 16 * 16;
                }
            }
        }
        png_write_row(png_ptr, row);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

int main(int argc, char** argv) {
    struct timespec start_t, end_t;
    clock_gettime(CLOCK_MONOTONIC, &start_t);
    double CPU_TIME, COMM_TIME;
    CPU_TIME = COMM_TIME = 0;

    /* MPI Initialize */
    int rank, size;
	MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // printf("From process %d out of %d\n", rank, size);

    /* Get # of CPU */
    cpu_set_t cpu_set;
    omp_set_nested(true);
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    int num_threads = CPU_COUNT(&cpu_set);
    // printf("%d cpus available\n", CPU_COUNT(&cpu_set));

    /* Argument parsing */
    assert(argc == 9);
    const char* filename = argv[1];
    int iters = strtol(argv[2], 0, 10);
    double left = strtod(argv[3], 0);
    double right = strtod(argv[4], 0);
    double lower = strtod(argv[5], 0);
    double upper = strtod(argv[6], 0);
    int width = strtol(argv[7], 0, 10);
    int height = strtol(argv[8], 0, 10);

    /* Variables */
    int i, root, row, col;
    double x0,y0;
    double dx = (right - left) / (double)width;
    double dy = (upper - lower) / (double)height; 

    // how many height divide into all process and should one node(process) get
    int height_per_process = (double)height / size;

    // remain how many height that are looking for allocate evenly to each process
    int height_remain = height % size;

    // total height for each process should get
    int n_height = height_per_process + (rank<height_remain?1:0);
    
    // each process get how many data
    int start_row, end_row, prev_start_row, prev_end_row;
    if(rank==0){
        start_row = 0; // 0*height_per_process
        end_row = start_row + height_per_process + (rank<height_remain?1:0);
    }
    else{
        start_row = rank*height_per_process + ((rank-1)<height_remain?rank:height_remain);
        end_row = start_row + height_per_process + (rank<height_remain?1:0);
    }
    // printf("Start_Row:%d | End_Row:%d\n", start_row, end_row);

    /* allocate memory for image */
    int* root_image = (int*)malloc(width * height * sizeof(int));
    assert(root_image);
    int* image = (int*)malloc(width * n_height * sizeof(int));
    assert(image);
    int* displs = (int*)malloc(size * sizeof(int));
    assert(displs);
    int* rcounts = (int*)malloc(size * sizeof(int));
    assert(rcounts);

    // double new_lower;
    // if(rank==0)
    //     new_lower = lower;
    // else{
    //     new_lower = start_row*dy + lower;
    // }
    
    #pragma omp parallel num_threads(num_threads)
    {   
        #pragma omp for schedule(dynamic)
            for (row=0; row<n_height; ++row) {
                //y0 = row*dy + lower; <-- wrong
                double y0 = (rank + size*row)*dy + lower;
                for (col=0; col<width; ++col) {
                    double x0 = col*dx + left;
                    int repeats = 0;
                    double x = 0;
                    double y = 0;
                    double length_squared = 0;
                    while (repeats < iters && length_squared < 4) {
                        double temp = x * x - y * y + x0;
                        y = 2 * x * y + y0;
                        x = temp;
                        length_squared = x * x + y * y;
                        ++repeats;
                    }
                    image[row * width + col] = repeats;
                }
            }
    }

    // Prepare displs, rcounts for each stride of MPI_Gatherv, further more details, pls google by yourself
    displs[0] = 0; //initial memory location before stride
    for(i=0;i<size;i++){
        if(i<height_remain)
            rcounts[i] = (height_per_process+1)*width;
        else
            rcounts[i] = height_per_process*width;

        if(i+1<size)
            displs[i+1] = displs[i] + rcounts[i];
    }
    
    clock_gettime(CLOCK_MONOTONIC, &end_t);
    struct timespec temp_t = diff(start_t, end_t);
    CPU_TIME += temp_t.tv_sec + (double)temp_t.tv_nsec/1000000000.0;
    

    clock_gettime(CLOCK_MONOTONIC, &start_t);
    MPI_Gatherv(image, n_height*width, MPI_INT, root_image, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    clock_gettime(CLOCK_MONOTONIC, &end_t);
    temp_t = diff(start_t, end_t);
    COMM_TIME += temp_t.tv_sec + (double)temp_t.tv_nsec/1000000000.0;

    
    if(rank==0){
        clock_gettime(CLOCK_MONOTONIC, &start_t);
        int* final_image = (int*)malloc(width * height * sizeof(int));
        int index = 0;
        for(i=0;i<n_height;i++){
            int stride = n_height;
            int check_out_of_remainder = 0;
            for(row=i; row<height && index<height*width; row+=stride){
                for(col=0; col<width; col++){
                    final_image[index] = root_image[row*width + col];
                    index++;
                }
                if(check_out_of_remainder >= height_remain){
                    stride = height_per_process;
                }
                check_out_of_remainder++;
            }
        }
        clock_gettime(CLOCK_MONOTONIC, &end_t);
        temp_t = diff(start_t, end_t);
        CPU_TIME += temp_t.tv_sec + (double)temp_t.tv_nsec/1000000000.0;
        write_png(filename, iters, width, height, final_image);
        free(final_image);
    }
    
    printf("RANK:%d | TOTAL_TIME:%lf | CPU_TIME:%lf | COMM_TIME:%lf\n", rank, CPU_TIME+COMM_TIME, CPU_TIME, COMM_TIME);
    free(displs);
    free(rcounts);
    free(image);
    free(root_image);
    MPI_Finalize();
}



