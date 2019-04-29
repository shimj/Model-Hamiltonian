## Model Py DOC

#### Update note

* 2019-3-22
  * 修改了pgroup下的query.py的match_info函数的注释，添加了`(op3x3 is "gamma_nu" in the note, and it should be noted that the value of time reversal operator is a minus identity rather than the identity.)`
  * 判断P和T是否同时存在，改为判断PT是否存在。并且不再判断名字，而是判断op3x3和is_au。

#### Code Structure

* model_hamiltonian
  * pgroup (package)
    * json (folder): 存储点群的json文件
    * get (module)
      * get_data函数: 以dict类型返回整个json数据的copy
    * query (module): json数据接口
      * get_the_first_irred_rep函数: 根据轨道和排布查询得到不可约表示，如果用特征标来查找，或者单单用轨道、单单用排布来查找，则返回第一个匹配的轨道和排布对应的不可约表示。如果使用了含有TR的json数据来查找，那么返回的矩阵list中最后一个是TR的表示矩阵。另外，如果用特征标来查找，那么除了生成元的特征标外，还要提供给出恒元的特征标（也即不可约表示的维数），恒元的特征标要放在list的最前面。
      * match_info函数: 处理多（或一个）能级，得到一个大的可约（或不可约）表示，然后再与其他信息（生成元的名字、3x3矩阵、是否反幺正）一起打包（a list of dict）后返回。核心部分具体来说，接收a list of “轨道和排布的组合”，或者a list of “特征标列”，对每个元素调用get_the_first_irred_rep函数（其实是分情况调用了multi_orbs_to_rep和multi_chars_to_rep，再间接调用了上面这个函数）得到多个不可约表示表示，然后对多个不可约表示做直和，即得所有带的总的表示。另外，这个函数还接收一个可选参数：基矢的过渡矩阵，如果输入了过渡矩阵，那么会把最后那个直和矩阵做相似变换后再返回。（注意一点，如果基矢进行了两次变换，最后的过渡矩阵是第一个过渡矩阵乘第二个过渡矩阵，别搞反了，乘完再作为参数传入函数。）
  * model_hamiltonian (module): 核心模块
    * model_detail函数: 核心函数，返回的信息很详细
    * simple_calc函数: 调用model_detail函数，把输出信息变得更简单实用。输入轨道信息或者特征标信息（数组类型），数组元素个数等于不可约表示数（不考虑偶然简并的话就是能级数）；返回一个数组（a），每个元素都是一个tuple（a[0]），元素个数由系数非零的Hermitian基矢个数决定。每个tuple包含三个元素，第一个元素是从1开始的序号：`a[0][0]=1,a[1][0]=2...`；第二个元素是系数（k的多项式）；第三个元素是Hermitian矩阵。注：最终的Hamiltonian不会输出，需要拿返回的结果自己算，就是$\sum_i a[i][1]*a[i][2]$。
    * energy函数: 传入一个包含多个Hermitian矩阵的list，得到由这些矩阵组合出来的矩阵的本征值的一般表达式，h1、h2、h3 ... 被假设为这些矩阵的系数。（如果无法得到本征值，计算会消耗很长时间，并给出无意义的分段函数结果，因此在调用此函数时可以考虑使用try，以便手动终止计算。）