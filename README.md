# 程序说明

## 前言
欢迎来到我的程序:-)！当你来到这里时，说明你在核专业课程设计这门课中可能遇到了一些困难。我是在大四下的时候学习这门课的（没错，就是大四下，大家都除了毕设无所事事时，我们有大作业，还是1-2周要做完），当时正值毕业设计开题报告审查，于是不可避免地要忙碌。前人栽树，后人乘凉，我也是参考网络的资源，自己完成了此门课程的任务。但无奈网络资源鱼龙混杂，且其中的代码包含不少错误，有时甚至误人子弟。于是乎，本着资源互惠共享的精神，我这里有一份成型的报告和程序代码可供你参考。此程序为我个人独立完成，并且尽量给予了良好的注释，方便参考。希望能够给你帮助！

## 文件夹介绍
4. `NTH-caculation`：matlab工作文件夹，所有的程序均在其中。
3. `fig`：结果图以及功率归一化因子源数据图。
2. `report_files`：我的报告。
1. `refprop_install`：热物性查询程序文件，具体使用说明在文件夹里。

## 关于功率归一化因子
当时为了确定功率归一化因子，我费了一番功夫。最后确定采用图解源数据加高斯拟合的方法得到功率归一化因子的生成公式。坐标图数据读取的网页为[WebPlotDigitizer](https://apps.automeris.io/wpd/index.zh_CN.html)。提取教程[参考这篇](https://blog.csdn.net/YanLu99/article/details/114172184)。

## 关于堆芯热功率
相信细心的你一定发现了，code中的堆芯热功率写的是2800MW，与指导书上的3400MW相比有所出入。没错，经过测试，3400MW的热功率会使得堆芯冷却剂过沸，导致堆芯爆炸。指导书并不严谨，这令我失望。我没有向老师指出这个问题，老师也没有发现code的端倪，说明老师对此门课也没有深入的理解。如果你有兴趣，可以与老师探讨~
