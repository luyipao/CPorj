#include <iostream>
#include <thread>
#include <vector>
#include <functional> // for std::ref

// 这是你要让每个线程执行的任务
void infinite_loop() {
    for (int i = 0; ; ++i) {
        i += 1;
        // 为了防止编译器优化掉这个看似无用的循环，可以加一点操作
        // 但对于大多数现代编译器，死循环本身就不会被优化掉
    }
}

int main() {
    // 获取你的电脑支持多少个并发线程（通常是核心数或逻辑核心数）
    const unsigned int num_threads = std::thread::hardware_concurrency();

    std::cout << "你的电脑支持 " << num_threads << " 个并发线程。" << std::endl;
    std::cout << "现在将创建 " << num_threads << " 个线程来跑满CPU..." << std::endl;

    // 创建一个线程容器
    std::vector<std::thread> threads;

    // 创建并启动线程
    for (unsigned int i = 0; i < num_threads; ++i) {
        threads.emplace_back(infinite_loop);
    }

    // 等待所有线程执行完毕（在这个例子里，它们永远不会结束）
    for (auto& t : threads) {
        t.join();
    }

    return 0;
}