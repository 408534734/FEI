# FEI
#### Formula and Equation Interperter
#### 表达式、方程解释器
## 功能简介
因为在数值计算编程中往往涉及到在打码中插入大量的公式，而在将公式键入代码时非常容易出现各种错误，所以开发了本解释器，可以将要键入代码的公式用另外的文档以**MarkDown**格式和**LaTeX数学公式**格式存于`FEI.md`中，在预览中实现非常人性化的可视化
##使用方法
- 添加头文件`FEI.h`和程序文件`FEI.cpp`
- 添加将编译程序`FEI.py`与`*.cpp`在统一目录下
- 在编译前前使用`python3 FEI.py 输入文件 输出`命令将`*.cpp`进行预编译的预编译
- 强烈建议直接使用与本项目一样的目录结构，直接相应的文件夹内添加代码文件，然后直接`make`编译即可
## 文件结构
- `bin` 存放编译后的可执行文件
    - `main` 示例程序的可执行文件
- `build` 存放编译出的可链接文件
    - `*.o` 可链接文件
- `doc`
    - `FEI.html` FEI.md的预览文件
    - `FEI_files` FEI.md的预览文件缓存文件夹
    - **`FEI.md` 存放数学公式的文件**
    - `希腊字母.xls` 可识别的希腊字母
- `include`
    - `FEI.h` FEI对象的头文件
- `src`
    - `main.cpp` 主程序
    - `FEI.cpp` FEI对象的函数体
    - `FEI.py` FEI系统的编译器
- `ans_output.m` 示例程序的输出结果
- `result.jpg` 示例程序输出结果的作图
- `makefile` 嗯，这是makefile
- `readme.md` 嗯，这就是我
