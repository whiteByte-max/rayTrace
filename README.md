# rayTrace
how to implement a simplified ray tracing demo based on smallpt learning

包含了一些光线追踪学习及实践的项目

1、smallpt_custom 显而易见，参考smallpt为范例实现的一个简单的光线追踪渲染器，基本和原版无异，按面向对象结构调整了代码，主要用于我入门，做的过程中将 光线追踪的定义、发展、类别以及基本的算法实现过了一遍。
顺便重温了 积分、线代、概率相关的数学知识。此外，做的过程中对信号采样相关知识通过书籍做了一遍基本学习，和对其中一些基本点深入看了些文章。也为后面在这个领域深入打下最基本的基础。

每个子像素采样10次的效果：
![image](https://github.com/user-attachments/assets/25dbc701-e70d-4804-ab8e-722d04301c6a)
每个子像素采样50次的效果：
![image](https://github.com/user-attachments/assets/2658c45d-044c-4aa5-b84c-173376c3b983)
每个子像素采样200次的效果：
![image](https://github.com/user-attachments/assets/180c4de9-b658-43e1-905a-2e58d5cb5a85)


基本上可以看到，随着子像素采样次数的增加，画面整体亮度在增长，而且画面看起来更精致了些，细碎的噪声感觉也少了很多，基本符合预期，但就是效率不太行，实际上还是有很多无效采样，以及过程是可以加速的，这个随着我后面的学习继续改进吧。
