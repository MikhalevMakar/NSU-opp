//
// Created by Макар Михалёв on 17.05.2023.
//

#ifndef TASK5_CONCURRENTQUEUE_H
#define TASK5_CONCURRENTQUEUE_H

#include <atomic>
#include <queue>
#include <random>
#include <utility>


template <typename T>
class ConcurrentQueue {
private:
    std::queue<T> queue;
    pthread_mutex_t* lock;
    pthread_cond_t* cond;

public:
    ConcurrentQueue(pthread_mutex_t* lock, pthread_cond_t* cond) {
        this->lock = lock;
        this->cond = cond;
    }

    void push(T task) {
        pthread_mutex_lock(lock);
        queue.add(task);
        pthread_cond_signal(cond);
        pthread_mutex_unlock(lock);
    }

    T pop() {
        pthread_mutex_lock(lock);
        T task = queue.front();
        queue.pop();
        pthread_mutex_unlock(lock);

        return task;
    }

    void filling(int size_queue, int initial_boundary, int final_boundary) {
         pthread_mutex_lock(lock);
         for (int i = 0; i < size_queue; ++i) {

         queue.push(
                    get_repeat_number(get_random_value(BOUNDS_SIZE_TASK, rank, count_process),
                                      rank,
                                      count_process)
                    );
         }

         pthread_mutex_unlock(lock);
    }

    bool empty() {
        pthread_mutex_lock(lock);
        bool isEmpty = queue.empty();
        pthread_mutex_unlock(lock);

        return isEmpty;
    }

    ssize_t size() {
         pthread_mutex_lock(lock);
         ssize_t size = queue.size();
         pthread_mutex_unlock(lock);

         return size;
    }

private:
    ssize_t get_repeat_number(const ssize_t random_value, const int rank, const int countProcess) {
        return abs(rank - (random_value % countProcess));
    }

    ssize_t get_random_value(const int bounds, const int rank, const int count_process) {
        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_int_distribution<ssize_t> dist(count_process, (rank + 1) * bounds + count_process);
        return dist(gen);
    }
};

#endif //TASK5_CONCURRENTQUEUE_H
