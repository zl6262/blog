| 数学分析           | 线性代数          | 概率论与数理统计 | 复变函数     | 大学物理 | 离散数学 |
| :----------------- | ----------------- | ---------------- | ------------ | -------- | -------- |
| *信号与系统*  8.11 | 数字信号处理 8.11 | 模拟电子电路     | 数字电子电路 |          |          |
| 通信原理 8.11      | 信息论 8.11       | 移动通信 8.11    |              |          |          |
| 计算机系统原理     | 数据结构          | C++              | 数据库原理   |          |          |
| 电磁场与电磁波     | 微波技术基础      |                  |              |          |          |

# 数学分析

#### 1.

> 

# 信号与系统

### 第一章 绪论

##### Question 1.信号的分类

> 模拟信号：幅值和时间都为连续的信号
>
> 数字信号：幅值和时间都为离散的信号
>
> 离散信号：

##### Question 2.信号的运算

> $$
> 微分：\frac{d}{dt}f(t)		积分：\int_{-\infty}^{t}f(\delta)d\delta
> $$
>
> 

##### Question 3.冲激函数的定义

> $$
> 矩形脉冲：\delta(t)=lim_{\tau\rightarrow0}\frac{1}{\tau}[u(t+\frac{\tau}{2})-u(t-\frac{\tau}{2})]
> $$
>
> $$
> Sa(t)信号： lim_{k\rightarrow\infty}[\frac{k}{\pi}Sa(kt)]
> $$

##### Question 4.

> 积分单元(延迟单元）、相加、倍乘

##### Question 5.

线性系统：齐次性：叠加性和均匀性

时不变系统：系统的参数不随时间变化; 系统响应与激励施加于系统的时刻无关

因果性

Question 6.

### 第二章 连续时间系统的时域分析

##### Question 1.时域经典法解微分方程

齐次解+特解（系数由初始状态确定）

Question 2.起始点的跳变

在没有冲击电流（或阶跃电压）强迫作用于电容的条件下，电容两端电压Vc(t)不发生跳变；在没有冲激电压（或阶跃电流）强迫作用于电感的条件下，流经电感的电流IL(t)不发生跳变。

Question 3.

零输入响应：齐次解；零状态响应：齐次解+特解。

自由响应：齐次解；强迫相应：特解。

稳态响应；暂态响应

Question 4.

冲激响应：t>0, e(t)=0；冲激响应引入能量存储作用，响应形式与零输入响应相同（齐次解）。确定系数：奇异函数系数匹配。

阶跃响应：应增加特解项； 
$$
齐次性：\sum y(t+k\Delta t)\Delta t
$$
Question 5.卷积：限于线性时不变系统；线性时变系统，冲激响应由冲激加入时间\tau和响应观测时间t决定

Question 6.反卷积

求逆系统。

Question 7.算子表示微分方程
$$
Ae^{at}u(t)=\frac{A}{p-a}\delta(t)
$$

### 第三章 傅里叶变换

Question 1.傅里叶级数存在条件(充分条件)

Dirichlet条件

* 在一周期内，间断点和极大值、极小值的个数为有限个
* 一周期内信号是绝对可积的

Question 2.傅里叶有限级数和最小方均误差
$$
误差函数\epsilon_N(t)=f(t)-S_N(t)
$$

$$
方均误差 E_N=\overline{{\epsilon_N(t)}^2}=\frac{1}{T_1}\int_{t_0}^{t_0+T_1}\epsilon_N^2(t)dt
$$

Question 3.典型周期信号的傅里叶级数
$$
周期矩形脉冲信号:
$$
Question 4.傅里叶变换

频谱密度函数 F(w)

充分条件：在无限区间内满足绝对可积条件
$$
E[u(t+\tau/2)-u(t-\tau/2)]  \leftrightarrow E\tau Sa(w\tau/2)
$$

$$
f(0)\tau =F(0) \qquad F(0)B=2\pi f(0)\qquad B=\frac{2\pi}{\tau}\qquad信号的等效脉冲宽度与占有的等效带宽成反比
$$

Question 5.周期信号的傅里叶变换

在允许冲击函数存在并认为它具有意义的前提下， 绝对可积条件就成为不必要的限制了。
$$
f(t)=\sum_{n=-\infty}^{\infty}F_ne^{jnw_1t}\qquad F_n=\frac{1}{T_1}\int_{-\frac{T_1}{2}}^{\frac{T_1}{2}}f(t)e^{-jnw_1t}dt
$$

$$
从周期性脉冲序列f(t)中截取一个周期，F_0(w)=\int_{-\frac{T_1}{2}}^{\frac{T_1}{2}}f(t)e^{-jwt}dt
$$

$$
F_n=\frac{1}{T_1}F_0(w)|_{w=nw_1}
$$

利用单脉冲的傅里叶变换式可以很方便地求出周期性脉冲序列的傅里叶系数。

Question 6.抽样
$$
抽样函数：Sa(t)=\frac{sin\space t}{t}
$$
自然抽样（矩形脉冲抽样）

理想抽样
$$
时域抽样\space T_s=\frac{1}{2f_m} \qquad
频域抽样 f_s=\frac{1}{2t_m}
$$

### 第四章 拉普拉斯变换、连续时间系统的s域分析

##### Question 1.拉式变换要点

将微分、积分运算转换为乘除运算；将积分、微分方程转换为代数方程。

扩大了傅里叶函数的种类

从零、极点分析系统特性。
$$
拉式变换存在条件\qquad lim_{t\rightarrow\infty}f(t)e^{-\sigma t}=0 \space(\sigma>\sigma_0)
$$

##### Question 2.拉式变换的性质

$$
初值定理 \space f(0_+)=lim_{s\rightarrow\infty}[sF(s)-ks]
$$

$$
终值定理\space lim_{t\rightarrow\infty}f(t)=lim_{s\rightarrow0}sF(s) \qquad终值定理只有在sF(s)在s平面的虚轴上及其右边都解析时(原点除外)才可使用。
$$

$$
s\rightarrow0相当于直流状态，所以得到电路的稳定的终值f(\infty)；而s\rightarrow\infty，相当于接入信号的突变(高频分量), 所以给出相应的初值f(0_+)
$$

Question 3. 拉普拉斯逆变换

##### Question 4.系统函数H(s)

回路电压方程 
$$
V=ZI\qquad  Z:各回路的s域互阻抗或自阻抗
\\第k个回路电流:I_k(s)=\frac{\Delta_{1k}}{\Delta}V_1(s)+\frac{\Delta_{2k}}{\Delta}V_2(s)+···\frac{\Delta_{lk}}{\Delta}V_l(s)
\\\Delta为Z方阵的行列式，\Delta_{jk}是\Delta的元素Z_{jk}的代数余子式
$$
H(p)可以用来说明零状态和零输入特性，而H(s)只能用来说明零状态特性。

##### Question 5.

系统与最小相移网络

> s
