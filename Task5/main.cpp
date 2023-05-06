#include <iostream>
#include <mpi.h>
#include <queue>
#include <random>
#include <pthread.h>


enum constants {
    BOUNDS_QUEUE = 10,
    BOUNDS_SIZE_TASK = 10000
};

int iterationCounter = 0;
double globalRes = 0;

int get_repeat_number(const int random_value, const int rank, const int countProcess) {
    return abs(rank - (random_value % countProcess));
}

int get_random_value(const int bounds, const int rank) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<int> dist(1, (rank + 1) * bounds);
    return dist(gen);
}

std::queue<int> initialize_queue_task(const int size_queue, const int rank, const int countProcess) {
    std::queue<int> queue;

    for(int i = 0; i < size_queue; ++i) {
        queue.push(get_repeat_number(get_random_value(BOUNDS_SIZE_TASK, rank),  rank, countProcess));
    }
}

void* execute_task(void* queue) {
    pthread_mutex_t mutex;
    auto queue_tasks = (std::queue<int>*)queue;

    pthread_mutex_lock(&mutex);
    int task_extent;
    while(true) {
        pthread_mutex_lock(&mutex);

        task_extent = queue_tasks->front();
        queue_tasks->pop();

        for (int i = 0; i < task_extent; ++i) {
            globalRes += sqrt(i);
        }

        pthread_mutex_unlock(&mutex);
    }

}


//void* run_worker(void* queue) {
//
//    auto queue_tasks = (std::queue<int>*)queue;
//    int task_extent;
//
//    while(true) {
//        int size_queue = queue_tasks->size();
//        for (int i = 0; i < size_queue; ++i) {
//            task_extent = queue_tasks->front();
//            queue_tasks->pop();
//            execute_task(task_extent);
//        }
//
//    }
//}

void run_pthread(const int rank, const int countProcess) {
    int size_queue = get_random_value(BOUNDS_QUEUE, rank);
    std::queue<int> queue_tasks = initialize_queue_task(size_queue, rank, countProcess);

    pthread_mutex_t mutex;
    pthread_mutexattr_t mutex_attrs;

    pthread_mutexattr_init(&mutex_attrs);
    pthread_mutex_init(&mutex, &mutex_attrs);
    pthread_mutexattr_destroy(&mutex_attrs);

    pthread_attr_t attrs;

    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

    pthread_t thread_execute;
    pthread_t thread_request_send;  // отправляет сообщение
    pthread_t thread_request_recv; // принимает запрашивает

    pthread_create(&thread_execute, &attrs, execute_task, &queue_tasks);

    while(true) {

    }

    //pthread_create(&thread_send, &attrs, vec_partial_product, &queue_tasks);
    //pthread_create(&thread_receive, &attrs, vec_partial_product, &context[i]);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, countThread;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &countThread);

    run_pthread(rank, countThread);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
