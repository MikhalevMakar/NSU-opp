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
    BOUNDS_QUEUE = 100,
    BOUNDS_SIZE_TASK = 100000,
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

    void filling(const int size_queue, const T initial_boundary, const T final_boundary) {
         pthread_mutex_lock(lock);
         for (int i = 0; i < size_queue; ++i) {

         queue.push(
                    get_repeat_number(initial_boundary,
                                      final_boundary)
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

    unsigned size() {
         pthread_mutex_lock(lock);
         unsigned size = queue.size();
         pthread_mutex_unlock(lock);

         return size;
    }

private:
    T get_repeat_number(const T initial_boundary, const T final_boundary) {
        return abs(initial_boundary - (final_boundary / get_random_value(initial_boundary, final_boundary)));
    }

    ssize_t get_random_value(const T initial_boundary, const T final_boundary) {
        std::random_device rd;
        std::mt19937 gen(rd());

        std::uniform_int_distribution<ssize_t> dist(initial_boundary, final_boundary);
        return dist(gen);
    }
};

#endif //TASK5_CONCURRENTQUEUE_H
