#include <iostream>
#include <mpi.h>
#include <pthread.h>
#include "ConcurrentQueue.h"

typedef struct {
    const int CountThread;
    pthread_mutex_t* Lock;
    std::unique_ptr<ConcurrentQueue<int>> Queue;
    const int Rank;
    pthread_cond_t* CondWait;
    pthread_cond_t* CondWork;
    bool StatusRun;
    int CountTaskExecute;
    double ResultSum;
} Context;

void* task_send(void* _context) {
    auto context = (Context*) _context;
    MPI_Status status;
    int rankSender;

    printf("RANK: %d, Start task_send\n", context->Rank);

    while(context->StatusRun) {
        MPI_Recv(&rankSender, 1, MPI_INT, MPI_ANY_SOURCE, TAG_REQUEST_TASK, MPI_COMM_WORLD, &status);

        if(rankSender == STOP_WORK) continue;

        int task = context->Queue->pop();
        MPI_Send(&task, 1, MPI_INT, rankSender, TAG_SEND_TASK, MPI_COMM_WORLD);
    }

    printf("RANK: %d, finish task_send, status: %d\n", context->Rank, context->StatusRun);
    return nullptr;
}

void* task_wait(void* _context) {
    auto context = (Context *) _context;

    MPI_Status status;
    int recvTask;

    printf("RANK: %d, start task_wait\n", context->Rank);

    while(context->StatusRun) {

        while (!context->Queue->empty()) {

             pthread_mutex_lock(context->Lock);
             pthread_cond_signal(context->CondWork);
             pthread_cond_wait(context->CondWait, context->Lock);
             pthread_mutex_unlock(context->Lock);
        }

        int number_proc_completed = 0;

        for (int i = 0; i < context->CountThread; ++i) {

            MPI_Send(&context->Rank, 1, MPI_INT, i, TAG_REQUEST_TASK, MPI_COMM_WORLD);
            MPI_Recv(&recvTask, 1, MPI_INT, i, TAG_SEND_TASK, MPI_COMM_WORLD, &status);
            if (recvTask != PROCESS_FULFILLED_TASK) {
                context->Queue->push(recvTask);
            } else { ++number_proc_completed; }
        }

        if(number_proc_completed == context->CountThread) {
            recvTask = STOP_WORK;
            MPI_Barrier(MPI_COMM_WORLD);
            context->StatusRun = false;
            MPI_Send(&recvTask, 1, MPI_INT, context->Rank, TAG_REQUEST_TASK, MPI_COMM_WORLD);

            pthread_mutex_lock(context->Lock);
            pthread_cond_signal(context->CondWork);
            pthread_mutex_unlock(context->Lock);
        }
    }

    printf("RANK: %d, finish task_wait, status: %d\n", context->Rank, context->StatusRun);
    return nullptr;
}

void task_execute(ssize_t task_extent, Context* context) {
    for (int j = 0; j < task_extent; ++j) {
        context->ResultSum += sqrt(j);
    }
}

void* worker(void* _context) {
    auto context = (Context *) _context;

    while (context->StatusRun) {
        while(true) {
            int task = context->Queue->pop();
            if (task == STOP_WORK) break;
            task_execute(task, context);
            ++context->CountTaskExecute;
        }

        pthread_mutex_lock(context->Lock);
        pthread_cond_signal(context->CondWait);
        pthread_mutex_unlock(context->Lock);

        while(context->Queue->empty() && context->StatusRun) {
          pthread_mutex_lock(context->Lock);
          pthread_cond_wait(context->CondWork, context->Lock);
          pthread_mutex_unlock(context->Lock);
        }
    }

    printf("RANK: %d, finish worker, count task_execute: %d\n", context->Rank, context->CountTaskExecute);
    return nullptr;
}

Context fill_context(const int count_process, const int rank, pthread_mutex_t* mutex,
                     pthread_cond_t* cond_wait, pthread_cond_t* cond_work) {

    Context context = {
            .CountThread = count_process,
            .Lock = mutex,
            .Queue = std::make_unique<ConcurrentQueue<int>>(mutex, cond_work),
            .Rank = rank,
            .CondWait = cond_wait,
            .CondWork = cond_work,
            .StatusRun = true,
            .CountTaskExecute = 0,
            .ResultSum = 0
    };

    return context;
}

void run_pthread(const int rank, const int count_process) {
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&mutex, nullptr);

    pthread_attr_t attrs;
    pthread_attr_init(&attrs);
    pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE);

    pthread_t thread_worker;
    pthread_t thread_task_wait;
    pthread_t thread_task_send;

    pthread_cond_t cond_wait = PTHREAD_COND_INITIALIZER,
                   cond_work = PTHREAD_COND_INITIALIZER;

    pthread_cond_init(&cond_wait, nullptr);
    pthread_cond_init(&cond_work, nullptr);

    Context context = fill_context(count_process, rank, &mutex, &cond_wait, &cond_work);
    if (context.Rank == ROOT) {
        context.Queue->filling((rank + INCREMENT) * BOUNDS_QUEUE, (rank + INCREMENT)*100, BOUNDS_SIZE_TASK);
    }

    printf("RANK: %d, SIZE QUEUE : %u\n", context.Rank, context.Queue->size());

    MPI_Barrier(MPI_COMM_WORLD);

    pthread_create(&thread_worker, &attrs, worker, &context);
    pthread_create(&thread_task_wait, &attrs, task_wait, &context);
    pthread_create(&thread_task_send, &attrs, task_send, &context);

    pthread_join(thread_worker, nullptr);
    pthread_join(thread_task_wait, nullptr);
    pthread_join(thread_task_send, nullptr);

    pthread_attr_destroy(&attrs);
    pthread_cond_destroy(&cond_wait);
}

int main(int argc, char **argv) {
    int provider;
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