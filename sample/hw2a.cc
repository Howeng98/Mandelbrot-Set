//##############################################################################################//
/* Don't seperate evenly the whole height to all the pthread, it will cause some thread may become idle if the height computation is same enough and create a lot of idle time!!! */
//##############################################################################################//
// if(temp_args->thread_id != 0){
//     start_row = temp_args->thread_id * n_per_thread 
//     + (temp_args->thread_id - 1 < remainder ? (temp_args->thread_id) : remainder);
//     end_row = start_row + n_per_thread + (temp_args->thread_id < remainder ? 1 : 0);
// }
// else{
//     start_row = temp_args->thread_id * n_per_thread;
//     end_row = start_row + n_per_thread + (temp_args->thread_id < remainder ? 1 : 0);
// }
//##############################################################################################//

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
#include <time.h>

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER; //creates a mutex variable
int current_row;

typedef struct{
    double left, right, upper, lower;
    unsigned int width;
    unsigned int height;
    int thread_id;
    int num_threads;
    int num_iter;
    int* image;
} Worker_Args;

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

void* worker(void* args){
    struct timespec worker_start_t, worker_end_t;
    double CPU_TIME=0;
    clock_gettime(CLOCK_MONOTONIC, &worker_start_t);
    Worker_Args* temp_args = (Worker_Args*)args;
    // printf("This is worker %d\n", temp_args->thread_id);

    int start_row=0, end_row;
    int n_per_thread = temp_args->height/temp_args->num_threads;
    int remainder = temp_args->height%temp_args->num_threads;
    int row, col, i, j;
    int repeats;
    double dx = (temp_args->right - temp_args->left) / temp_args->width;
    double dy = (temp_args->upper - temp_args->lower) / temp_args->height;
    double x, y, x0, y0, length_squared, t; 

    
    while(true){
        // Each pthread keep taking height to their working stage, make sure utilizing
        // mutex_lock section on height
        pthread_mutex_lock(&mutex);
        if(current_row < temp_args->height){
            start_row = current_row;
            current_row++;
        }
        else{
            start_row = current_row;
        }
        pthread_mutex_unlock(&mutex);
        
        if(start_row < temp_args->height){
            y0 = start_row * dy + temp_args->lower;
        
            for(col=0;col<temp_args->width;++col){
                x0 = col * dx + temp_args->left;
                repeats = 0;
                x = 0;
                y = 0;
                length_squared = 0;

                while (repeats < temp_args->num_iter && length_squared < 4) {
                    t = x * x - y * y + x0;
                    y = 2 * x * y + y0;
                    x = t;
                    length_squared = x * x + y * y;
                    ++repeats;
                }
                temp_args->image[start_row * temp_args->width + col] = repeats;
            }
        }
        
        if(start_row == temp_args->height){
            break;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &worker_end_t);
    struct timespec temp_t = diff(worker_start_t, worker_end_t);
    CPU_TIME = temp_t.tv_sec + (double)temp_t.tv_nsec/1000000000.0;
    printf("Thread_ID:%d | CPU_TIME:%lf\n", temp_args->thread_id,CPU_TIME);
    pthread_exit(NULL);
}

int main(int argc, char** argv) {
    struct timespec start_t, end_t;
    clock_gettime(CLOCK_MONOTONIC, &start_t);
    
    double CPU_TIME, IO_TIME;
    
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    // printf("%d cpus available\n", CPU_COUNT(&cpu_set));
    int num_threads = CPU_COUNT(&cpu_set);
    // printf("%d num threads\n", num_threads);

    /* argument parsing */
    assert(argc == 9);
    const char* filename = argv[1];
    int iters = strtol(argv[2], 0, 10);
    double left = strtod(argv[3], 0);
    double right = strtod(argv[4], 0);
    double lower = strtod(argv[5], 0);
    double upper = strtod(argv[6], 0);
    int width = strtol(argv[7], 0, 10);
    int height = strtol(argv[8], 0, 10);
    int i;

    /* allocate memory for image */
    int* image = (int*)malloc(width * height * sizeof(int));
    assert(image);
    current_row = 0;
    pthread_mutex_init(&mutex, NULL);
    pthread_t threads[num_threads];
    Worker_Args args[num_threads];

    /* assign argument */
    for(i=0;i<num_threads;i++){
        args[i].left = left;
        args[i].right = right;
        args[i].lower = lower;
        args[i].upper = upper;
        args[i].width = width;
        args[i].height = height;
        args[i].num_threads = num_threads;
        args[i].num_iter = iters;
        args[i].thread_id = i;
        args[i].image = image;
    }

    for(i=0;i<num_threads;i++){
        pthread_create(&threads[i], NULL, worker, &args[i]);
    }

    for(i=0;i<num_threads;i++){
        pthread_join(threads[i], NULL);
    }

    clock_gettime(CLOCK_MONOTONIC, &end_t);

    // calculate CPU_TIME
    struct timespec temp_t = diff(start_t, end_t);
    CPU_TIME = temp_t.tv_sec + (double)temp_t.tv_nsec/1000000000.0;
    // printf("CPU_TIME:%lf\n", CPU_TIME);

    /* draw and cleanup */
    write_png(filename, iters, width, height, image);
    free(image);
}
