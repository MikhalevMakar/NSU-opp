//
// Created by Макар Михалёв on 17.05.2023.
//

#ifndef TASK5_CONCURRENTQUEUE_H
#define TASK5_CONCURRENTQUEUE_H

#include <atomic>
#include <queue>
#include <random>
#include <utility>

enum constants {
    BOUNDS_QUEUE = 1000,
    BOUNDS_SIZE_TASK = 10000000,
    ROOT = 0,
    PROCESS_FULFILLED_TASK = -1,
    TAG_REQUEST_TASK = 123,
    TAG_SEND_TASK = 321,
    INCREMENT = 1,
    STOP_WORK = -1
};

template <typename T>
class ConcurrentQueue {
private:
    std::queue<T> queue;
    pthread_mutex_t* lock{};
    pthread_cond_t* cond{};

public:
    ConcurrentQueue(pthread_mutex_t* lock, pthread_cond_t* cond) {
        this->lock = lock;
        this->cond = cond;
    }

    void push(const T& task) {
        pthread_mutex_lock(lock);
        queue.push(task);
        pthread_cond_signal(cond);
        pthread_mutex_unlock(lock);
    }

    T pop() {
        T task;
        pthread_mutex_lock(lock);
        if (!queue.empty()) {
             task = queue.front();
             queue.pop();
        } else {
            task = PROCESS_FULFILLED_TASK;
        }
        pthread_mutex_unlock(lock);
        return task;
    }

    bool empty() {
        pthread_mutex_lock(lock);
        bool isEmpty = queue.empty();
        pthread_mutex_unlock(lock);

        return isEmpty;
    }

    unsigned size() {
         pthread_mutex_lock(lock);
         unsigned size = queue.size();
         pthread_mutex_unlock(lock);

         return size;
    }
};

#endif //TASK5_CONCURRENTQUEUE_H
