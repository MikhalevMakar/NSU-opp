#include <iostream>
#include <mpi.h>
#include <queue>
#include <random>
#include <utility>
#include <pthread.h>
#include "ConcurrentQueue.h"

enum constants {
    BOUNDS_QUEUE = 10,
    BOUNDS_SIZE_TASK = 10000,
    ROOT = 0,
    PROCESS_FULFILLED_TASK = -1,
    TAG_REQUEST_TASK = 123,
    TAG_SEND_TASK = 321
};

typedef struct {
    const int CountThread;
    pthread_mutex_t* Lock;
    ConcurrentQueue<ssize_t>* Queue;
    const int Rank;
    pthread_cond_t Cond;
    bool StatusRun;
    int Count_task_execute;
} Context;

double globalRes = 0;

void* task_send(void* _context) {
    auto context = (Context*) _context;
    MPI_Status status;
    int rankSender;
    ssize_t task;

    while(context->StatusRun) {
        MPI_Recv(&rankSender, 1, MPI_INT, MPI_ANY_SOURCE, TAG_REQUEST_TASK, MPI_COMM_WORLD, &status);

        pthread_mutex_lock(context->Lock);

        if (!context->Queue->empty()) {
            task = context->Queue->pop();
        } else {
            task = PROCESS_FULFILLED_TASK;
        }

        pthread_mutex_unlock(context->Lock);

        MPI_Send(&task, 1, MPI_INT, TAG_SEND_TASK, MPI_ANY_TAG, MPI_COMM_WORLD);
    }
    return nullptr;
}

void* task_wait(void* _context) {
    auto context = (Context *) _context;

    MPI_Status status;
    ssize_t recvTask;

    pthread_mutex_lock(context->Lock);

    while(context->StatusRun) {

        pthread_cond_wait(&context->Cond, context->Lock);

        int number_proc_completed = 0;

        std::cout << "Get up thread task wait\n";

        for (int i = 0; context->CountThread; ++i) {
            if (i == context->Rank) continue;

            MPI_Send(&context->Rank, 1, MPI_INT, i, TAG_REQUEST_TASK, MPI_COMM_WORLD);
            MPI_Recv(&recvTask, 1, MPI_INT, i, TAG_SEND_TASK, MPI_COMM_WORLD, &status);

            if(recvTask != PROCESS_FULFILLED_TASK) {
                context->Queue->push(recvTask);
            } else { ++number_proc_completed; }
        }

        if(number_proc_completed == context->CountThread)
            context->StatusRun = false;
    }

    pthread_mutex_unlock(context->Lock);
    return nullptr;
}

void task_execute(ssize_t task_extent) {
    for (int j = 0; j < task_extent; ++j) {
        globalRes += sqrt(j);
    }
}

void* worker(void* _context) {
    auto context = (Context *) _context;

    while (context->StatusRun) {
        std::cout << "SIZE QUEUE : \n" << context->Queue.size();
        pthread_mutex_lock(context->Lock);

        unsigned size_queue = context->Queue->size();
        for(int i = 0; i < size_queue; ++i) {
            task_execute(context->Queue->pop());
        }

        pthread_cond_signal(&context->Cond);

        pthread_mutex_unlock(context->Lock);
    }

    return nullptr;
}


Context fill_context(const int count_process, const int rank,
                     pthread_mutex_t* mutex, pthread_cond_t* cond) {

    Context context = {
            .CountThread = count_process,
            .Lock = mutex,
            .Queue = new ConcurrentQueue<ssize_t>(mutex, cond),
            .Rank = rank,
            .Cond = PTHREAD_COND_INITIALIZER,
            .StatusRun = true,
            .Count_task_execute = 0
    };
    
    return context;
}

void run_pthread(const int rank, const int count_process) {
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_attr_t attrs;

    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

    pthread_t thread_worker;
    pthread_t thread_task_wait;
    pthread_t thread_task_send;

    Context context = fill_context(count_process, rank, &mutex, cond);
    context.Queue->filling(rank*BOUNDS_QUEUE);

    pthread_create(&thread_worker, &attrs, worker, &context);
    pthread_create(&thread_task_wait, &attrs, task_wait, &context);
    pthread_create(&thread_task_send, &attrs, task_send, &context);

    pthread_attr_destroy(&attrs);

    pthread_join(thread_worker, nullptr);
    pthread_join(thread_task_wait, nullptr);
    pthread_join(thread_task_send, nullptr);

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