#pragma once

#include <queue>
#include <mutex>
#include <condition_variable>

template <typename type>
class fifo {
private:
    std::queue<type> q;
    std::mutex mtx;
    std::condition_variable cv;
public:
    // void write(type& elem) {
    //     std::lock_guard<std::mutex> lock(mtx);
    //     q.push(elem);
    //     cv.notify_one();
    // }

    void write(type elem) {
        std::lock_guard<std::mutex> lock(mtx);
        q.push(elem);
        cv.notify_one();
    }

    type read(void) { // blocking read
        std::unique_lock<std::mutex> lock(mtx);
        if (q.empty()) { // main & callback threads, no spurious wakeup
            cv.wait(lock);
        }
        type elem = q.front();
        q.pop();
        return elem;
    }

    bool peek(type& elem) { // non-blocking attempt to read
        std::lock_guard<std::mutex> lock(mtx);
        if (q.empty()) {
            return false;
        } else {
            elem = q.front();
            q.pop();
            return true;
        }
    }

    unsigned size(void) {
        std::lock_guard<std::mutex> lock(mtx);
        return q.size();
    }
};
