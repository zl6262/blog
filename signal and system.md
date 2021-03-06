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

##### Question 5. 全通网络与最小相移网络

* 全通网络：系统的极点位于左半平面，零点位于右半平面，且零点与极点对于jw轴互为镜像。在传输系统中常用来进行相位校正，如相位均衡器或移相器。
* 最小相移网络：零点仅位于左半平面或jw轴的网络。 
* 非最小相移函数可以表示为最小相移函数与全通函数的乘积。(网络则对应级联)

##### Question 6. 线性系统的稳定性
* 一种是由极点类型进行判断。
* 另一种是BIBO(有限输入，有限输出)。

##### Question 7. 双边拉普拉斯变换
注意由收敛域进行反变换的形式。

##### Question 8. 拉普拉斯变换与傅里叶变换的关系
收敛边界位于虚轴时：
$$
F(s)=F_a(s)+\sum_{n=1}^N \frac{K_n}{s-jw_n}
$$ 
$$
\begin{aligned}
F[f(t)] &=F_a(jw)+F[\sum_{n=1}^N K_ne^{jw_nt}u(t) \\
      &=F_a(jw)+\sum_{n=1}^N K_n \{\delta(w-w_n)*[\pi\delta(w)+\frac{1}{jw}]\} \\
      &=F_a(jw)+\sum_{n=1}^N \frac{K_n}{j(w-w_n)} + \sum_{n=1}^N K_n\pi\delta(w-w_n)\\
      &=F(s)|_{s=jw} +\sum_{n=1}^N K_n \pi \delta(w-w_n)
\end{aligned}
$$

### 第五章 傅里叶变换应用于通信系统  ----滤波、调制与抽样

##### Question 1. 无失真传输

$$
r(t)=Ke(t-t_0)
$$

$$
R(jw)=KE(jw)e^{-jwt_0}
$$

要保证没有相位失真，必须使响应中各频率分量与激励中各对应分量滞后相同的时间。 

群时延: 系统相频特性对频率的导数并取负号。

用群时延间接表达相位特性的好处是便于实际测量，有助于理解调幅波传输过程的波形变化。

##### Question 2. 理想低通滤波器

$$
|H(jw)|=\left\{\begin{matrix} 
1 \qquad(-w_c<w<w_c)
\\0 \qquad(w为其他值)
\end{matrix} \right.
$$

$$
h(t)=\frac{w_c}{\pi}Sa[w_c(t-t_0)]
$$

理想低通的阶跃响应的上升时间与系统的截止频率（带宽）成反比。

##### Question 3. 系统的物理可实现性、佩利-维纳准则

* 时域：因果性，h(t)=0    t<0
* 频率：H(jw)满足平方可积条件  (功率有限信号)

佩利-维纳准则
$$
\int_{-\infty}^{\infty}\frac{|\space ln|H(jw)|\space|}{1+w^2}dw<\infty
$$

##### Question 4.利用希尔伯特变换研究系统函数的约束特性

因果信号
$$
h(t)=h(t)u(t)
$$

$$
R(w)=X(w)*\frac{1}{\pi w}
$$

$$
X(w)=R(w)*(-\frac{1}{\pi w})
$$

##### Question 5.从抽样信号恢复连续时间信号

1. 冲激抽样信号恢复连续时间信号
   $$
   f(t)=\sum f(nT)Sa[w_c(t-nT_s)]
   $$

2. 零阶保持抽样

   为复原F(w)频谱，接收端应用具有补偿特性的低通滤波器，存在相位滞后（时间上滞后T_s/2）

3. 一阶保持抽样：三角脉冲

   内插补偿低通滤波器，相移特性为零（非因果系统）。

##### Question 6. 数字通信

数字信号经中继器转发，噪声不会积累，且存在整形再生。

##### Question 7. 码速与带宽

矩形脉冲信号频谱函数的第一个零点1/T。 B=f=1/T。

B与矩形脉冲的宽度有关。

f与周期有关。

### 第六章 信号的矢量空间分析

##### Question 1. 范数

正定性、正齐性、三角形不等式

p阶范数
$$
||x||_p=\left\{
\begin{matrix}	
[\sum_{i=1}^{N}||x_i||^p]^{1/p}\space(对于1\leq p<\infty)
\\ max_{1\leq i\leq N} \space|x_i| \space(对于p\rightarrow\infty)
\end{matrix}\right.
$$
二阶范数的平方表示信号的能量

##### Question 2.能量谱与功率谱

维纳辛钦：相关函数的傅氏变换为功率谱密度。

##### Question 3.匹配滤波器

为了滤波器输出端的信号瞬时功率与噪声平均功率之比值为最大。

增强信号分量，减弱噪声分量。
$$
H(jw)=k[S(jw)e^{jwt_m}]^*
$$

$$
H(jw)=kS(-jw)e^{-jwt_m}
$$

$$
h(t)=ks(t_m-t)
$$

$$
t_m应尽可能小,使判决迅速,t_m=T比t_m>T合适。
$$

匹配滤波器是线性系统的最佳滤波器。

### 第七章 离散时间系统的时域分析

##### Question 1. 常系数线性差分方程的方法

* 迭代法
* 时域经典法----齐次解+特解 再用初始条件求系数
* 分别求零状态响应与零输入响应（初始条件）
* 变换域方法

### 第八章 z变换、离散时间系统的z域分析

##### Question 1. 围线积分

$$
X(z)=\sum x(n)z^{-n} \qquad \qquad x(n)=\frac{1}{2\pi j}\oint X(z)z^{n-1}dz
$$

$$
x(n)=\sum Res[X(z)z^{n-1}]_{z=z_m}
$$

$$
Res[X(z)z^{n-1}]_{z=z_m}=\frac{1}{(s-1)!}\{\frac{d^{s-1}}{dz^{s-1}}[(z-z_m)^sX(z)z^{n-1}]\}_{z=z_m}
$$

##### Question 2. 系统分析

离散时间系统稳定的充分必要条件是单位样值响应h(n)绝对可和。

对于稳定系统H(z)的收敛域应包含单位圆在内。

##### Question 3. 序列的傅里叶变换DTFT

x(n)逆变换积分限 -pi --- pi 是围线积分的变形。

x(n)绝对可和只是傅里叶变换存在的充分条件。

##### Question 4. 频率响应

$$
e^{sT}=e^{jwT} \longrightarrow z^{-n}=e^{-nsT}=e^{-jnwT}
$$

$$
重复频率 \qquad w_s*T=2\pi \qquad w_s=\frac{2\pi}{T}
\\ -\frac{w_s}{2} \leq w \leq \frac{w_s}{2}
$$

##### Question 5. 全通系统

全通系统零点与极点的模量互为倒数、幅角相等。

冲激响应不变法由于s - z平面多值映射，不能用于高通和带阻滤波器

### 第九章 离散傅里叶变换以及其他离散正交变换

##### Question 1.时域与频域对应关系

$$
时域采样间隔T_1 -- 频域周期 w_s
\\频域采样间隔f_1--时域周期T_s
\\X(k)=W^{nk}x(n)
\\x(n)=\frac{1}{N}W^{-nk}X(k)
$$

矩形脉冲序列的DFT仅在k=0样点取得N值，在其余(N-1)个点都为零。
$$
X(k)=N\delta(k)
$$

##### Question 2.离散傅里叶变换的性质

对序列的位移赋予新的解释：
$$
先将x(n)周期延拓成x_p(n)，然后移位得到x_p(n-m)，最后取主值，得到了x(n)的所谓圆周移位序列x_p(n-m)R_N(n)。
\\有限长序列x(n)的圆周移位序列一般写作x((n-m))_NR_N(n)
$$

##### Question 3. 卷积

圆周卷积只在[0,N-1]区间内进行。
$$
y(n)=\sum_{m=0}^{N-1}x(m)h((n-m))_NR_N(n)
$$

线性卷积：
$$
y(n)=\sum x(m)h(n-m)
\\ 0\leq m \leq N-1 \qquad 0\leq n-m \leq M-1
\\ 0 \leq n \leq N+M-2
$$
圆卷积的两序列长度相等。补零扩展，这样圆卷积时，向右移去的零值，从左边出现仍取零值。

补零扩展后的长度L满足
$$
L \geq N+M-1
$$
圆卷积可借助FFT技术。

##### Question 4. 内插

$$
x(n)=\frac{1}{N}\sum X(k)W^{-nk}
\\X(z)=\sum [\frac{1}{N}\sum X(k)W^{-nk}]z^{-n}
$$

##### Question 5. FFT

1. 计算圆周卷积，两次FFT，一次IFFT
2. 即位运算；自然顺序（需要额外寄存器)
3. 快速卷积：对于卷积序列之一特别长，无法全都补零到相同长度，则采取分段卷积。重叠保留法，重叠相加法。
4. 快速相关（功率谱计算）：两次FFT、一次IFFT计算相关函数（求信号的时延）。

##### Question 6. 对连续时间信号分析的逼近（信号频谱分析）

1. 时间有限信号：抽样导致频谱混叠（aliasing）。
2. 频率有限信号：时域截短（加窗），频域与Sa函数卷积，频谱泄露(leakage)。
3. 周期信号：若为频谱无限的周期信号，经抽样后也会产生频谱混叠。若频谱有限，则可避免。

##### Question 7. 电平数值的相对变动

从DFT计算结果转至连续信号的频谱时，需考虑电平数值的相对变动。

连续信号的频谱函数计算抽样信号的频谱时，需倍乘系数1/T_s。

1. 对一个时间有限连续信号进行傅里叶分析，将DFT计算结果乘以系数T_s，则可得到其近似频谱。
2. 由频谱合成波形。如果已知某信号的频谱在正、负频谱范围内共占据频带f_s，利用IDFT计算之结果乘以系数f_s即可获得其近似的时间波形。
3. 用DFT来求其一周期函数的傅里叶级数近似式。这时，因子1/N最好是放在正变换式中，而不是在逆变换式中，采取这一措施后，无论正、逆DFT变换都可以直接表示所需结果，无需再乘转换系数。

##### Question 8.OFDM

有效性：提高码元速率。

可靠性：噪声干扰和传输信道非理想特性引入的失真。

均衡：寻找与信道模型相应的逆系统，使信道与均衡器组合提供逼近理想传输特性的系统模型。

OFDM：多载波后，每个子载波传输的数据率只是原有单载波系统的1/N。每个子载波承载的码元周期是原单载波系统码元周期的N倍，使每个子载波传输码元之周期可以远远大于信道时延之扩展，提高抗码间干扰能力。

### 第十章 模拟与数字滤波器

##### Question 1. 模拟滤波器的逼近

工作参数法

1.给定频率特性模平方|H_a(j\Omega)|^2，求系统函数H_a(s)
$$
h_a(t)=L^{-1}[H_a(s)]为t的实函数
\\则H_a(j\Omega)具有共轭对称性，H_a(j\Omega)=H_a^*(-j\Omega)
\\|H_a(j\Omega)|^2|_{j\Omega=s}=H_a(j\Omega)H_a^*(j\Omega)=H_a(j\Omega)H_a(-j\Omega)=H_a(s)H_a(-s)
$$
2.若限定最小相位，则只能取所有左半平面的零极点。

##### Question 2. Butterworth滤波器

最大平坦性（0点处各阶导数为0）
$$
|H_a(j\Omega)|^2=\frac{1}{1+(\frac{\Omega}{\Omega_c})^{2N}}
$$
N为奇数，极点的幅角为2\pi /(2N)的整数倍。

N为偶数，极点的幅角为\pi /2N的奇数倍。

##### Question 3. 稳定系统

$$
H_a(s)=\frac{B(s)}{A(s)}=\frac{b_ms^m+b_{m-1}s^{m-1}+···+b_1s+b_0}{a_ns^n+a_{n-1}s^{n-1}+···+a_1s+a_0}
$$
稳定的必要条件：分母是霍尔维茨多项式（其零点再s平面左半平面内或j\Omega轴上的单阶）。

稳定系统的另一限制：
$$
H_a(s)|_{s\rightarrow\infty} \approx \frac{b_m}{a_n}s^{m-n}
$$

1. 若m>n，则\infty是H_a(s)的极点，其阶次为m-n，无穷极点可认为位于j\Omega轴上，因而按照稳定条件，m-n最多为1。
2. 同理s=0也在虚轴上，所以H_a(s)分母的最低阶次与分子的最低阶次之差也不能超过一次。

##### Question 4. IIR滤波器（无限冲激响应数字滤波器）

数字滤波器是一种对输入信号进行离散时间处理的系统。输入信号可以为连续信号、抽样序列或数字信号。

IIR的冲激响应h(n)为无限长。
$$
H(z)=\frac{\sum_{r=0}^Mb_rz^{-r}}{1+\sum_{k=1}^Na_kz^{-k}}
\\若a_k \neq0，则构成IIR滤波器，反而为FIR滤波器
$$

1. 频响曲线以2\pi为周期重复出现（一般选取T=1)
2. 完成s域与z域之间的映射转换。

##### Question 5. 冲激响应不变法

$$
h(n)=h_a(nT)=h_a(t)|_{t=nT}
\\H_a(s)=\sum \frac{A_k}{s-s_k}
\\h(n)=h_a(t)|_{t=nT}=\sum A_ke^{s_knT}u(n)
\\H(z)=\sum \frac{A_k}{1-e^{s_kT}z^{-1}}
$$

$$
\frac{1}{s-s_k} \longrightarrow \frac{1}{1-e^{s_kT}z^{-1}}
$$

极点s_k映射到z平面是位于z=e^{s_k T}处的极点。
$$
L[h_a(t)\sum\delta(t-nT)]=\sum h_a(nT)e^{-nsT}=\sum h(n)z^{-n}|_{z=e^{sT}}=H(z)|_{z=e^{sT}}
$$

$$
L[h_a(t)\sum \delta(t-nT)]=L[h_a(t)·\frac{1}{T}\sum e^{j\frac{2\pi}{T}kt}]=\frac{1}{T}\sum H_a(s+j\frac{2\pi}{T}k)
$$

$$
H(z)|_{z=e^{sT}}=\frac{1}{T}\sum H_a(s+j\frac{2\pi}{T}k)
$$

z平面与s平面为多值映射关系。

只适用于低通滤波器或限带的高通或带通场合。
$$
h(n)的设计值采用Th_a(nT),以保持转换后数字滤波器增益不变。
$$

##### Question 6. 双线性变换法

把s域压缩至0到+-\pi之间，形成单值映射。
$$
s平面整个变换到一个中介的s_1平面的一个水平窄带\Omega:-\frac{\pi}{T} \rightarrow \frac{\pi}{T}之中，再经过z=e^{s_1T}进行映射。
$$
过程:
$$
\Omega=tan(\frac{\Omega_1T}2)
\\通过欧拉公式
\\s=\frac{e^{j\frac{s_1T}{2}}-e^{-j\frac{s_1T}{2}}}{e^{j\frac{s_1T}2}+e^{-j\frac{s_1T}2}}=\frac{1-e^{-s_1T}}{1+e^{-s_1T}}
\\z=e^{s_1T}
\\s=\frac{1-z^{-1}}{1+z^{-1}}
\\为使AF与DF存在频率对应关系
\\s=\frac{2}{T} \frac{1-z^{-1}}{1+z^{-1}}
\\ z=\frac{\frac{T}2+s}{\frac{T}2-s}
\Omega
$$
进行预畸变校正。同时考虑引入的相位畸变是否被允许。

##### Question 7. FIR滤波器

$$
H(z)=\sum_{n=0}^{N-1}h(n)z^{-n}
$$

有N-1个零点，且为有限制，全部极点位于z平面的原点，因此系统是稳定的。

FIR可以做到严格的线性相移，而IIR滤波器的相移特性呈现非线性，这是FIR滤波器最突出的优点。

|            |         |                                  |
| ---------- | ------- | -------------------------------- |
| h(n)偶对称 | N为奇数 | 低通                             |
|            | N为偶数 | 无法实现高通、带阻特性           |
| h(n)奇对称 | N为奇数 | 无法实现低通、高通及带阻滤波特性 |
|            | N为偶数 | 无法实现低通、带阻滤波器         |

* h(n)只要满足偶对称或奇对称条件，他的相频特性就是现行的，且相移常数\alpha=(N-1)/2。在h(n)为奇对称情况，滤波器有固定的90°相移，在微分器、希尔伯特变换器及信号正交处理中特别有用。
* 

##### Question 8. 窗函数法

$$
h_d(n)=\frac{1}{2\pi}\int_{-\pi}^{\pi}H_d(jw)e^{jwn}dw
$$

$$
\epsilon^2=\frac{1}{2\pi}\int_{-\pi}^{\pi}|H_d(e^{jw})-H(e^{jw})|^2dw=min
$$

频域均方误差最小的意义下进行逼近。
$$
h(n)=h_d(n)R_N(n) \qquad R_N(n):矩形窗函数
$$

$$
H_d(e^{jw})= e^{jw\alpha},|w|\leq w_c
$$

$$
h_d(n)=\frac{sin[w_c(n-\alpha)]}{\pi(n-\alpha)} \qquad \alpha=\frac{N-1}{2}
$$

$$
h_d(n)无限长、非因果，h(n)长度为N,具有偶对称的线性相位FIR滤波器单位样值响应。
$$

$$
加窗截短后的H_g(w)特性：在w_c附近形成过渡带，过渡带两边出现正、负肩峰（肩峰的增量值约为9\%,Gibbs现象）
$$

所以出现各种窗函数。

##### Question 9. 比较IIR、FIR

* IIR可以用较少的阶数得到很好的选频特性，但是用相位特性的非线性作为代价的。
* IIR滤波器用递归结构，反馈支路的存在使系统对稳定性要求高，有限字长效应影响大，设计不当会引起震荡。
* FIR没有可控制的极点，因此要达到与IIR相当的选频特性，需要的阶数很高，导致信号通过系统延时增加。但是能做到严格的线性相移。
* FIR采用非递归结构，因此系统始终稳定，且有限字节影响小。
* 对相位要求不敏感的场合(如语音)，可选IIR滤波器；而对波形上携带信息的传输系统(如图像信号等)，对系统相位特性要求较高，可考虑用FIR滤波器。

### 第十一章 反馈系统

##### Question 1. 信号流图

梅森增益公式
$$
H=\frac{1}{\Delta}\sum_k g_k \Delta_k
\\ \Delta:流图的特征行列式
\\ \Delta=1-(所有不同环路的增益之和)+(每两个互不接触环路增益乘积之和)-(每三个互不接触环路增益乘积之和)+···
\\  =1-\sum_aL_a+\sum_{b,c}L_bL_c-\sum_{d,e,f}L_dL_eL_f
$$

$$
k:表示由源点到阱点之间第k条前向通路的标号
\\ g_k:表示由源点到阱点之间第k条前向通路的增益
\\ \Delta_k:称为对于第k条前向通路特征行列式的余因子。它是除去与第k条前向通路相接触的环路外，余下的特征行列式(或者说在\Delta式中只留下与该通路不接触者)。
$$

### 第十二章 系统的状态变量分析


