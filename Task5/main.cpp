#include <iostream>
#include <mpi.h>
#include <pthread.h>
#include "ConcurrentQueue.h"

enum constants {
    BOUNDS_QUEUE = 10,
    BOUNDS_SIZE_TASK = 10000,
    ROOT = 0,
    PROCESS_FULFILLED_TASK = -1,
    TAG_REQUEST_TASK = 123,
    TAG_SEND_TASK = 321,
    INCREMENT = 1,
    STOP_WORK = -1
};

typedef struct {
    const int CountThread;
    pthread_mutex_t* Lock;
    std::unique_ptr<ConcurrentQueue<ssize_t>> Queue;
    const int Rank;
    pthread_cond_t* Cond;
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
        if(rankSender == STOP_WORK) continue;
        if (!context->Queue->empty()) {
            task = context->Queue->pop();
        } else {
            task = PROCESS_FULFILLED_TASK;
        }

        MPI_Send(&task, 1, MPI_INT, rankSender, TAG_SEND_TASK, MPI_COMM_WORLD);
    }

    printf("finish task_send, status: %d\n", context->StatusRun);
    return nullptr;
}

void* task_wait(void* _context) {
    auto context = (Context *) _context;

    MPI_Status status;
    ssize_t recvTask;

    pthread_mutex_lock(context->Lock);

    while(context->StatusRun) {

        while (!context->Cond) {
            pthread_cond_wait(context->Cond, context->Lock);
        }

        pthread_mutex_unlock(context->Lock);

        int number_proc_completed = 0;

        std::cout << "Get up thread task wait\n";

        for (int i = 0; i < context->CountThread; ++i) {
            if (i == context->Rank) continue;

            MPI_Send(&context->Rank, 1, MPI_INT, i, TAG_REQUEST_TASK, MPI_COMM_WORLD);
            MPI_Recv(&recvTask, 1, MPI_INT, i, TAG_SEND_TASK, MPI_COMM_WORLD, &status);

            if (recvTask != PROCESS_FULFILLED_TASK) {
                context->Queue->push(recvTask);
            } else { ++number_proc_completed; }
        }

        if (number_proc_completed == context->CountThread - INCREMENT) {
            recvTask = STOP_WORK;
            MPI_Send(&recvTask, 1, MPI_INT, context->Rank, TAG_REQUEST_TASK, MPI_COMM_WORLD);
            context->StatusRun = false;
        }
    }

    printf("finish task_wait, status: %d\n", context->StatusRun);
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
        unsigned size_queue = context->Queue->size();
        for(int i = 0; i < size_queue; ++i) {
            task_execute(context->Queue->pop());
            ++context->Count_task_execute;
        }
       // pthread_mutex_lock(context->Lock);
        pthread_cond_signal(context->Cond);

//        while(!context->Cond) {
  //         pthread_cond_wait(context->Cond, context->Lock);
//        }
       // pthread_mutex_lock(context->Lock);
    }

    printf("finish worker, count task_execute: %d\n", context->Count_task_execute);
    return nullptr;
}


Context fill_context(const int count_process, const int rank,
                     pthread_mutex_t* mutex, pthread_cond_t* cond) {

    Context context = {
            .CountThread = count_process,
            .Lock = mutex,
            .Queue = std::make_unique<ConcurrentQueue<ssize_t>>(mutex, cond),
            .Rank = rank,
            .Cond = cond,
            .StatusRun = true,
            .Count_task_execute = 0
    };

    return context;
}

void run_pthread(const int rank, const int count_process) {
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, nullptr);

    pthread_attr_t attrs;
    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

    pthread_t thread_worker;
    pthread_t thread_task_wait;
    pthread_t thread_task_send;

    pthread_cond_t cond;
    pthread_cond_init(&cond, nullptr);

    Context context = fill_context(count_process, rank, &mutex, &cond);
    context.Queue->filling((rank+INCREMENT)*BOUNDS_QUEUE, (rank+INCREMENT), BOUNDS_SIZE_TASK);

    printf("SIZE QUEUE : %zd\n", context.Queue->size());

    pthread_create(&thread_worker, &attrs, worker, &context);
    pthread_create(&thread_task_wait, &attrs, task_wait, &context);
    pthread_create(&thread_task_send, &attrs, task_send, &context);



    pthread_join(thread_worker, nullptr);
    pthread_join(thread_task_wait, nullptr);
    pthread_join(thread_task_send, nullptr);

    pthread_attr_destroy(&attrs);
    pthread_cond_destroy(&cond);
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