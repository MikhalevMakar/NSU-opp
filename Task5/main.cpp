#include <iostream>
#include <mpi.h>
#include <queue>
#include <random>
#include <utility>
#include <pthread.h>

enum constants {
    BOUNDS_QUEUE = 10,
    BOUNDS_SIZE_TASK = 10000,
    ROOT = 0
};

typedef struct {
    const int CountThread;
    pthread_mutex_t *Lock;
    std::queue<int> Queue;
    const int Rank;
} Context;

int iterationCounter = 0;
double globalRes = 0;

int get_repeat_number(const int random_value, const int rank, const int countProcess) {
    return abs(rank - (random_value % countProcess));
}

int get_random_value(const int bounds, const int rank, const int count_process) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<int> dist(count_process, (rank + 1) * bounds + count_process);
    return dist(gen);
}

std::queue<int> initialize_queue_task(const int size_queue, const int rank, const int count_process) {
    std::queue<int> queue;

    for (int i = 0; i < size_queue; ++i) {

        queue.push(
                   get_repeat_number(get_random_value(BOUNDS_SIZE_TASK, rank, count_process),
                                     rank,
                                     count_process)
                   );
    }
}

void* task_execute(void* _context) {
    auto context = (Context *) _context;

    int task_extent;
    for(int i = 0; i < context->Queue.size(); ++i) {
        while (true) {
            pthread_mutex_lock(context->Lock);

            task_extent = context->Queue.front();
            context->Queue.pop();

            for (int j = 0; j < task_extent; ++j) {
                globalRes += sqrt(j);
            }

            pthread_mutex_unlock(context->Lock);

        }
        ++iterationCounter;
    }
}

void* task_wait(void* _context) { // если очередь пуста, то ищем задачки
    auto context = (Context *) _context;

    MPI_Status status;
    int recvTask;

    pthread_mutex_lock(context->Lock);

    while(true) {
        //pthread_cond_wait(&finishTasks, &sendThreadMutex);

        for (int i = 0; context->CountThread; ++i) {
            if (i == context->Rank) continue;

            MPI_Send(&context->Rank, 1, MPI_INT, context->Rank, i, MPI_COMM_WORLD);
            MPI_Recv(&recvTask, 1, MPI_INT, i, context->Rank, MPI_COMM_WORLD, &status);
        }

    }

    pthread_mutex_unlock(context->Lock);


    return NULL;
}

void* task_send(void* _context) {
    auto context = (Context *) _context;

    MPI_Status status;

    while (true) {


    }

    return NULL;
}

Context fill_context(const int count_process, const int rank,
                     pthread_mutex_t* mutex, std::queue<int> queue) {
    Context context = {
            .CountThread = count_process,
            .Lock = mutex,
            .Queue = std::move(queue),
            .Rank = rank

    };
    return context;
}

void run_pthread(const int rank, const int count_process) {
    int size_queue = get_random_value(BOUNDS_QUEUE, rank, count_process);
    std::queue<int> queue_tasks = initialize_queue_task(size_queue, rank, count_process);

    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    pthread_attr_t attrs;

    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

    pthread_t thread_execute;
    pthread_t thread_request_send;
    pthread_t thread_request_recv;

    Context context = fill_context(count_process, rank,
                                   &mutex, queue_tasks);

    pthread_create(&thread_execute, &attrs, task_execute, &context);
    pthread_create(&thread_request_send, &attrs, task_wait, &context);
    pthread_create(&thread_request_recv, &attrs, task_send, &context);

    pthread_attr_destroy(&attrs);

    pthread_join(thread_execute, nullptr);
    pthread_join(thread_request_send, nullptr);
    pthread_join(thread_request_recv, nullptr);

    pthread_mutex_destroy(&mutex);
}

int main(int argc, char **argv) {
    int provider = MPI_THREAD_MULTIPLE;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provider);


    int rank, countThread;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &countThread);

    double start_time = MPI_Wtime();

    run_pthread(rank, countThread);

    double end_time = MPI_Wtime();

    if (rank == ROOT)
        std::cout << std::endl << "TIME: " << end_time - start_time << " seconds" << std::endl;

    MPI_Finalize();
    return EXIT_SUCCESS;
}
