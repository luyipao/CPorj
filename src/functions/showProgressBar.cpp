#include "../include/showProgressBar.h"
void showProgressBar(int total, int current) {
    const int barWidth = 70; // 进度条宽度
    float progress = static_cast<float>(current) / total;
    int pos = barWidth * progress;

    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r"; // \r 使光标回到行首
    std::cout.flush(); // 刷新输出
}