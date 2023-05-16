#include <iostream>
#include <mpi.h>
#include <queue>
#include <random>
#include <utility>
#include <pthread.h>

enum constants {
    BOUNDS_QUEUE = 10,
    BOUNDS_SIZE_TASK = 10000,
    ROOT = 0,
    PROCESS_FULFILLED_TASK = -1
};

typedef struct {
    const int CountThread;
    pthread_mutex_t *Lock;
    std::queue<ssize_t> Queue;
    const int Rank;
    pthread_cond_t Cond;
    bool StatusRun;
} Context;

int iterationCounter = 0;
double globalRes = 0;

int get_repeat_number(const int random_value, const int rank, const int countProcess) {
    return abs(rank - (random_value % countProcess));
}

ssize_t get_random_value(const int bounds, const int rank, const int count_process) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<ssize_t> dist(count_process, (rank + 1) * bounds + count_process);
    return dist(gen);
}

std::queue<ssize_t> initialize_queue_task(const ssize_t size_queue, const int rank, const int count_process) {
    std::queue<ssize_t> queue;

    for (int i = 0; i < size_queue; ++i) {

        queue.push(
                   get_repeat_number(get_random_value(BOUNDS_SIZE_TASK, rank, count_process),
                                     rank,
                                     count_process)
                   );
    }
    return queue;
}

void task_execute(Context* context) {

}

void* worker(void* _context) {
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
    ssize_t recvTask;

    pthread_mutex_lock(context->Lock);

    printf("Run");
    while(context->StatusRun) {
        pthread_cond_wait(&context->Cond, context->Lock);
        int number_proc_completed = 0;

        for (int i = 0; context->CountThread; ++i) {
            if (i == context->Rank) continue;

            MPI_Send(&context->Rank, 1, MPI_INT, i, MPI_ANY_SOURCE, MPI_COMM_WORLD);
            MPI_Recv(&recvTask, 1, MPI_INT, i, MPI_ANY_SOURCE, MPI_COMM_WORLD, &status);

            if(recvTask != PROCESS_FULFILLED_TASK) {
                context->Queue.push(recvTask);
            } else { ++number_proc_completed; }
        }

        if(number_proc_completed == context->CountThread)
            context->StatusRun = false;
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
                     pthread_mutex_t* mutex, std::queue<ssize_t> queue) {
    Context context = {
            .CountThread = count_process,
            .Lock = mutex,
            .Queue = std::move(queue),
            .Rank = rank,
            .Cond = PTHREAD_COND_INITIALIZER,
            .StatusRun = true
    };
    return context;
    }

void run_pthread(const int rank, const int count_process) {
    ssize_t size_queue = get_random_value(BOUNDS_QUEUE, rank, count_process);
    std::queue<ssize_t> queue_tasks = initialize_queue_task(size_queue, rank, count_process);

    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_attr_t attrs;

    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

    pthread_t thread_worker;
    pthread_t thread_task_wait;
    pthread_t thread_task_send;

    Context context = fill_context(count_process, rank,
                                   &mutex, queue_tasks);

    //pthread_create(&thread_worker, &attrs, worker, &context);
    pthread_create(&thread_task_wait, &attrs, task_wait, &context);
    //pthread_create(&thread_task_send, &attrs, task_send, &context);

    pthread_attr_destroy(&attrs);

    //pthread_join(thread_worker, nullptr);
    pthread_join(thread_task_wait, nullptr);
    //pthread_join(thread_task_send, nullptr);

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
