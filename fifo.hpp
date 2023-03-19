#pragma once

#include <queue>
#include <mutex>
#include <condition_variable>

template <typename type>
class fifo: private std::queue<type> {
private:
    std::mutex mtx;
    std::condition_variable cv;
public:
    void write(type& elem) {
        std::lock_guard<std::mutex> lock(mtx);
        std::queue::push(elem);
        cv.notify_one();
    }

    type& read(void) { // blocking read
        std::unique_lock<std::mutex> lock(mtx);
        if (std::queue::empty()) { // main & callback threads, no spurious wakeup
            cv.wait(lock);
        }
        type& elem = std::queue::front();
        std::queue::pop();
        return elem;
    }

    bool peek(type& elem) { // non-blocking read
        std::lock_guard<std::mutex> lock(mtx);
        if (std::queue::empty()) {
            return false;
        }

        elem = std::queue::front();
        std::queue::pop();
        return true;
    }

    int size(void) {
        std::lock_guard<std::mutex> lock(mtx);
        return std::queue::size();
    }
};
