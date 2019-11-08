##————————————————————————————————缺失值处理——————————————————————————————————

###——————————————————————————背景知识————————————————————————————————————

####1.处理缺失值的步骤
##步骤：（1）识别缺失数据；（2）检查导致数据缺失的原因；（3）删除包含缺失值的实例或用合理的数值代替（插补）缺失值

##缺失值数据的分类：
#（1）完全随机缺失：若某变量的缺失数据与其他任何观测或未观测变量都不相关，则数据为完全随机缺失（MCAR）。
#（2）随机缺失：若某变量上的缺失数据与其他观测变量相关，与它自己的未观测值不相关，则数据为随机缺失（MAR）。
#（3）非随机缺失：若缺失数据不属于MCAR或MAR，则数据为非随机缺失（NIMAR）。

###理解缺失值数据的来由和影响
##识别缺失数据的数目、分布和模式有两个目的：
#（1）分析生成缺失数据的潜在机制；
#（2）评价缺失数据对回答实质性问题的影响。即：缺失数据的比例有多大？缺失数据是否集中在少数几个变量上，抑或广泛存在？缺失是随机产生的吗？缺失数据间的相关性或与可观测数据间的相关性，是否可以表明产生缺失值的机制呢？

#若缺失数据集中在几个相对不太重要的变量上，则可以删除这些变量，然后再进行正常的数据分析；
#若有一小部分数据随机分布在整个数据集中（MCAR），则可以分析数据完整的实例，这样仍可得到可靠有效的结果；
#若以假定数据是MCAR或MAR，则可以应用多重插补法来获得有铲的结论。
#若数据是NIMAR，则需要借助专门的方法，收集新数据，或加入一个相对更容易、更有收益的行业。


###2.识别缺失值
#NA：代表缺失值；
#NaN：代表不可能的值；
#Inf：代表正无穷；
#-Inf：代表负无穷。
#?is.na（）：识别缺失值；
#is.nan（）：识别不可能值；
#is.infinite（）：无穷值。


###（1）列表显示缺失值
#mice包中的md.pattern()函数可以生成一个以矩阵或数据框形式展示缺失值模式的表格
library(mice)
library(VIM)
data(sleep,package = "VIM");head(sleep);str(sleep) 
md.pattern(sleep)#如Missing data pattern所示，可知sleep一共有38个NA值。最左方，可见，完全没有缺失值的行数为42行，仅span有缺失值的行数有2个，如此类推。最右方，对应模式的行中，分别有几个缺失值，可见最高可达3个，其中都包括“Dream”、“NonD” 两列。最下方，显示出每列有几个缺失值。


###（2）图形探究缺失数据
#VIM包中提供大量能可视化数据集中缺失值模式的函数：aggr（）、matrixplot（）、scattMiss（）

library("VIM")  
dev.new();aggr(sleep,prop=FALSE,numbers=TRUE)  
dev.new();aggr(sleep,prop=TRUE,numbers=TRUE)#用比例代替了计数
dev.new();matrixplot(sleep) #matrixplot()函数可生成展示每个实例数据的图形

#marginplot（）函数可生成一幅散点图，在图形边界展示两个变量的缺失值信息。
library("VIM")  
dev.new();marginplot(sleep[c("Gest","Dream")],pch=c(20),col=c("darkgray","red","blue")) 


###（3）用相关性探索缺失值
#影子矩阵：用指示变量替代数据集中的数据（1表示缺失，0表示存在），这样生成的矩阵有时称作影子矩阵。
#求这些指示变量间和它们与初始（可观测）变量间的相关性，有且于观察哪些变量常一起缺失，以及分析变量“缺失”与其他变量间的关系。

x <- as.data.frame(abs(is.na(sleep)));x#sleep的影子矩阵
y <- x[which(sd(x)>0),]
cor(y)
cor(sleep,y,use = "pairwise.complete.obs") 
?cor

options(digits=1)
cor(na.omit(sleep)) 
cor(sleep,use="complete.obs")  
fit <- lm(Dream~Span+Gest,data=na.omit(sleep))
summary(fit) 


#pmm (Predictive Mean Matching)-For numeric variables
#logreg(Logistic Regression)—For Binary Variables( with 2 levels)
#polyreg(Bayesian polytomous regression)—For Factor Variables (>= 2 levels)
#Proportional odds model (ordered, >= 2 levels)