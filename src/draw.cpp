#include <fstream>
#include <cmath>

void plotToFile() {
    std::ofstream data("data.txt"), script("plot.plt");

    // 生成数据文件
    for (double x = 0; x <= 10; x += 0.1) {
        data << x << " " << sin(x) << "\n";
    }
    data.close();

    // 生成GNUplot脚本
    script << "set terminal png\n"
        << "set output 'plot.png'\n"
        << "plot 'data.txt' with lines title 'y=sin(x)'\n";
    script.close();

    // 执行绘图命令（需安装GNUplot）
    system("gnuplot plot.plt");
}

int main() {
    plotToFile();
    return 0;
}