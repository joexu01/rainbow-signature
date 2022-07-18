#define FINITE_FIELD_Q 257
#define LAYERS_U 5
#define VARIABLE_COUNT_N 33
#define VINEGAR_VARIABLE_V1 6

/* --------------！对一篇毕业设计论文的重新实现！------------- */

/*
	宏定义参数：
	FINITE_FIELD_Q      -- 有限域大小
	LAYERS_U            -- 层数，决定式(2.3)中的v_i的个数
	VARIABLE_COUNT_N    -- 变量总数n，也就是式(2.3)中的n
	VINEGAR_VARIABLE_V1 -- 宏定义重要的醋v_1
	详情见论文 3.2.1小节 油醋变量参数设计
*/

#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <vector>
#include <iostream>
#include <ctime>
#include <string>

using namespace std;
using namespace NTL;

// 如果在 main 函数里把 verbose 设置为 true，运行时会打印中间变量信息
bool verbose = false;

// 用来打印矩阵和向量
void print_matrix_with_info(string, int, mat_ZZ_p*);
void print_vector_with_info(string, int, vec_ZZ_p*);

class RainbowSignature {
public:
	vector<int> vinegar;        // 存储变量 v_i
	vector<int> oil;            // 存储变量 o_i

	/* -------------------------------------------------------------- */
	/*
		中心多项式 F 根据论文中 式(2.7) (2.8)生成，具体生成过程见 3.2.2 节
		需要注意的是，式(2.7)中指出，中心多项式中矩阵和向量的个数分别为 n-v_1
		个，但是在 3.2.2 节中，存储空间却开辟了 n 个，因此 mat_A 和 vec_B的
		前 v_1 个矩阵/向量是0矩阵/向量
	*/
	vector<mat_ZZ_p*> mat_A;    // 中心多项式 F 中的矩阵 A
	vector<vec_ZZ_p*> vec_B;    // 中心多项式 F 中的向量 B

	mat_ZZ_p* mat_M1;           // 可逆仿射变换 L_1 中的矩阵 M1 见式 (2.9)
	vec_ZZ_p* vec_C1;           // 可逆仿射变换 L_1 中的向量 C1 

	mat_ZZ_p* mat_M2;           // 可逆仿射变换 L_2 中的矩阵 M2 见式 (2.10)
	vec_ZZ_p* vec_C2;           // 可逆仿射变换 L_2 中的向量 C2
	/*
		以上便是私钥组成部分 L_1、L_2 以及 F ，见 3.2.3 节 私钥生成模块设计
	*/

	/* -------------------------------------------------------------- */
	vector<mat_ZZ_p*> mat_A_pp; // A''  --  对应式 (3.10) 中的变量命名
	vector<vec_ZZ_p*> vec_B_pp; // B''
	vector<ZZ_p> num_C_p;       // C'
	/*
		以上便是公钥式 (2.11) 表达成式 (3.10) 后的表示，公钥还包括有限域 Q
	*/

	RainbowSignature();
	void set_variables(vector<int>);
	void generate_coefficient_mat_vec();
	void generate_reversible_affine_transformation();
	void generate_public_key();
	vec_ZZ_p sign(vec_ZZ_p);

	void check(vec_ZZ_p, vec_ZZ_p);
};

// 构造函数，初始化一些变量的长度
RainbowSignature::RainbowSignature() {
	vinegar.resize(LAYERS_U);
	oil.resize(LAYERS_U - 1);

	vinegar[0] = VINEGAR_VARIABLE_V1;
	*(vinegar.rbegin()) = VARIABLE_COUNT_N;

	mat_A.resize(VARIABLE_COUNT_N);
	vec_B.resize(VARIABLE_COUNT_N);

	mat_M1 = NULL;
	mat_M2 = NULL;
	vec_C1 = NULL;
	vec_C2 = NULL;

	mat_A_pp.resize(VARIABLE_COUNT_N - VINEGAR_VARIABLE_V1);
	vec_B_pp.resize(VARIABLE_COUNT_N - VINEGAR_VARIABLE_V1);
	num_C_p.resize(VARIABLE_COUNT_N - VINEGAR_VARIABLE_V1);
}

// set_variables 方法设置 v_i 以及 o_i
void RainbowSignature::set_variables(vector<int> rest_v_i) {
	// 先把 v_1 和 v_n 分别设置为宏定义的 v_1 和 变量数 n
	vinegar[0] = VINEGAR_VARIABLE_V1;
	vinegar[LAYERS_U - 1] = VARIABLE_COUNT_N;

	// 然后再设置 v_2 到 v_(n-1) 
	int rest_idx = 0;
	for (int i = 1; i < LAYERS_U - 1; i++) {
		vinegar[i] = rest_v_i[rest_idx];
		rest_idx++;
	}

	// 根据 vinegar 生成 oil
	for (int i = 0; i < LAYERS_U - 1; i++) {
		oil[i] = vinegar[i + 1] - vinegar[i];
	}
}

// generate_coefficient_mat_vec 生成中心多项式矩阵的系数 A[n] 和向量 B[n]
void RainbowSignature::generate_coefficient_mat_vec() {
	for (auto v_val = vinegar.begin(); v_val < vinegar.end() - 1; v_val++) {
		int v_l = *v_val;
		int v_l_next = *(v_val + 1);

		// 前 v_1 个矩阵和向量置 0 不用
		for (int i = 0; i < VINEGAR_VARIABLE_V1; i++) {
			mat_A[i] = new mat_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N, VARIABLE_COUNT_N);
			vec_B[i] = new vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);
		}

		// 这个 for 循环的索引 k 是根据式 (2.7) 中的 k 设定的
		for (int k = v_l + 1; k <= v_l_next; k++) {
			mat_ZZ_p* temp_mat_ptr = new mat_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N, VARIABLE_COUNT_N);

			for (int i = 0; i < v_l; i++) {
				for (int j = i; j < v_l_next; j++) {
					(*temp_mat_ptr)[i][j] = rand() % FINITE_FIELD_Q;
				}
			}

			vec_ZZ_p* temp_vec_ptr = new vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);
			for (int i = 0; i < v_l_next; i++) {
				(*temp_vec_ptr)[i] = rand() % FINITE_FIELD_Q;
			}

			mat_A[k - 1] = temp_mat_ptr;
			vec_B[k - 1] = temp_vec_ptr;
		}
	}

	// 打印这个阶段生成的变量信息
	if (verbose) {
		for (int i = 0; i < VARIABLE_COUNT_N; i++) {
			print_matrix_with_info("matrix A", i, mat_A[i]);
			print_vector_with_info("vector B", i, vec_B[i]);
		}
	}
}

// generate_reversible_affine_transformation 生成可逆仿射变换 L_1 L_2 的矩阵和向量
// 详情见 3.2.3 节 私钥生成模块设计
void RainbowSignature::generate_reversible_affine_transformation() {
	int M1_len = VARIABLE_COUNT_N - vinegar[0];
	ZZ_p d = ZZ_p(0);

	// 初始化一个 M1_len × M1_len 大小的矩阵，M1_len 的大小是 n - v_1
	mat_M1 = new mat_ZZ_p(INIT_SIZE, M1_len, M1_len);
	while (true) {
		for (int i = 0; i < M1_len; i++) {
			for (int j = 0; j < M1_len; j++) {
				(*mat_M1)[i][j] = rand() % FINITE_FIELD_Q;
			}
		}

		// 计算随机生成的矩阵的行列式值，不为0就成功，为0就重新生成
		ZZ_p d_generate = determinant(*mat_M1);
		if (d_generate != d) {
			cout << "M1 determinant: " << d_generate << endl;
			break;
		}
		clear(*mat_M1);
	}

	// 初始化一个 M1_len 长度的向量
	vec_C1 = new vec_ZZ_p(INIT_SIZE, M1_len);
	for (int i = 0; i < M1_len; i++) {
		(*vec_C1)[i] = rand() % FINITE_FIELD_Q;
	}

	// 同理，但是 M2 的行和列都是 n 
	mat_M2 = new mat_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N, VARIABLE_COUNT_N);
	while (true) {
		for (int i = 0; i < VARIABLE_COUNT_N; i++) {
			for (int j = 0; j < VARIABLE_COUNT_N; j++) {
				(*mat_M2)[i][j] = rand() % FINITE_FIELD_Q;
			}
		}
		ZZ_p d_generate = determinant((*mat_M2));
		if (d_generate != d) {
			cout << "M2 determinant: " << d_generate << endl;
			break;
		}
		clear(*mat_M2);
	}

	vec_C2 = new vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);
	for (int i = 0; i < VARIABLE_COUNT_N; i++) {
		(*vec_C2)[i] = rand() % FINITE_FIELD_Q;
	}

	// 打印这个阶段生成的变量信息
	if (verbose) {
		print_matrix_with_info("matrix M1", 0, mat_M1);
		print_vector_with_info("vector C1", 0, vec_C1);
		print_matrix_with_info("matrix M2", 0, mat_M2);
		print_vector_with_info("vector C2", 0, vec_C2);
	}
}

// generate_public_key 生成公钥 A''(k) B''(k) C'(k)，见 3.2.4 节
void RainbowSignature::generate_public_key() {

	/*
		值得一提的是，此处使用了矩阵、向量相乘的 API -- mul(result, A, B)
		表示 result = A × B 
		还使用到了两个向量的内积 API -- InnerProduct(result, a, b)
		表示 result = a.b
	*/

	// 首先生成 A'(k)，根据式 (3.4) 和式 (3.5)
	vector<mat_ZZ_p*> mat_A_p(VARIABLE_COUNT_N - vinegar[0]);
	mat_ZZ_p temp_mat = mat_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N, VARIABLE_COUNT_N);
	
	mat_ZZ_p mat_M2_T = mat_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N, VARIABLE_COUNT_N);
	transpose(mat_M2_T, *mat_M2); // 对 M2 矩阵转置
	for (int k = vinegar[0]; k < VARIABLE_COUNT_N; k++) {
		mat_ZZ_p* result = new mat_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N, VARIABLE_COUNT_N);
		
		mul(temp_mat, mat_M2_T, *(mat_A[k]));  // M2_T * A[k]
		mul(*result, temp_mat, *mat_M2);       // 再乘以 M2
		mat_A_p[k - vinegar[0]] = result;
		clear(temp_mat);
	}
	
	// 根据式 (3.4) 计算 B'(k)
	vector<vec_ZZ_p*> vec_B_p(VARIABLE_COUNT_N - vinegar[0]);
	vec_ZZ_p temp_vec = vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);

	for (int k = vinegar[0]; k < VARIABLE_COUNT_N; k++) {
		vec_ZZ_p* result = new vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);
		mat_ZZ_p mat_A_k_T = mat_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N, VARIABLE_COUNT_N);
		transpose(mat_A_k_T, *(mat_A[k]));

		mul(*result, *vec_C2, mat_A_k_T);
		mul(*result, *result, *mat_M2);

		// vec_ZZ_p temp = vec_ZZ_p(INIT _SIZE, VARIABLE_COUNT_N);
		mul(temp_vec, *vec_C2, *(mat_A[k]));
		mul(temp_vec, temp_vec, *mat_M2);
		add(*result, *result, temp_vec);

		clear(temp_vec);
		mul(temp_vec, *(vec_B[k]), *mat_M2);
		add(*result, *result, temp_vec);

		vec_B_p[k - vinegar[0]] = result;
		clear(temp_vec);
	}
	
	// 根据式 (3.4) 计算 C(k)
	vector<ZZ_p> num_C(VARIABLE_COUNT_N - vinegar[0]);
	ZZ_p temp_num = ZZ_p(0);

	for (int k = vinegar[0]; k < VARIABLE_COUNT_N; k++) {
		ZZ_p sum = ZZ_p(0);

		mul(temp_vec, *vec_C2, *(mat_A[k]));
		InnerProduct(sum, temp_vec, *vec_C2);  // InnerProduct 代表向量内积

		InnerProduct(temp_num, *(vec_B[k]), *vec_C2);
		sum += temp_num;
		num_C[k - vinegar[0]] = sum;
		sum = ZZ_p(0);
		temp_num = ZZ_p(0);
	}

	// 打印这个阶段生成的变量信息
	if (verbose) {
		for (int i = 0; i < VARIABLE_COUNT_N - vinegar[0]; i++) {
			print_matrix_with_info("MAT A\'", i, mat_A_p[i]);
			print_vector_with_info("VEC B\'", i, vec_B_p[i]);
			cout << "C_num: " << num_C[i] << endl << endl;
		}
	}

	clear(temp_mat);
	// 根据式 (3.10) 生成 A''(k)
	for (int k = 0; k < VARIABLE_COUNT_N - vinegar[0]; k++) {
		mat_ZZ_p* sum = new mat_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N, VARIABLE_COUNT_N);
		for (int i = 0; i < VARIABLE_COUNT_N - vinegar[0]; i++) {
			mul(temp_mat, (*mat_M1)[k][i], *(mat_A_p[i]));
			(*sum) += temp_mat;
			clear(temp_mat);
		}
		mat_A_pp[k] = sum;
	}

	clear(temp_vec);
	// 根据式 (3.10) 生成 B''(k)
	for (int k = 0; k < VARIABLE_COUNT_N - vinegar[0]; k++) {
		vec_ZZ_p* sum = new vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);
		for (int i = 0; i < VARIABLE_COUNT_N - vinegar[0]; i++) {
			mul(temp_vec, (*mat_M1)[k][i], *(vec_B_p[i]));
			(*sum) += temp_vec;
			clear(temp_vec);
		}
		vec_B_pp[k] = sum;
	}

	temp_num = ZZ_p(0);
	// 根据式 (3.10) 生成 C'(k)
	for (int k = 0; k < VARIABLE_COUNT_N - vinegar[0]; k++) {
		ZZ_p sum = (*vec_C1)[k];
		for (int i = 0; i < VARIABLE_COUNT_N - vinegar[0]; i++) {
			temp_num = (*mat_M1)[k][i] * num_C[i];
			sum += temp_num;
			temp_num = ZZ_p(0);
		}
		num_C_p[k] = sum;
	}

	// 打印这个阶段生成的变量信息
	if (verbose) {
		for (int i = 0; i < VARIABLE_COUNT_N - vinegar[0]; i++) {
			print_matrix_with_info("MAT A\'\'", i, mat_A_pp[i]);
			print_vector_with_info("VEC B\'\'", i, vec_B_pp[i]);
			cout << "NUM C\': " << "index: " << i << " Value: " << num_C_p[i] << endl << endl;
		}
	}
}

// sign 对消息签名，返回值是消息的签名结果；见 3.2.5 节 消息签名模块设计
vec_ZZ_p RainbowSignature::sign(vec_ZZ_p msg) {
	long msg_length = msg.length();
	// 核对消息的长度是否为 n - v_1
	if (msg_length != long(VARIABLE_COUNT_N - vinegar[0])) {
		cout << "Error: length of message unmatch" << endl;
		exit(1);
	}

	// 根据式 (3.11) 求 Y'
	vec_ZZ_p vec_y_minus_C1 = msg - *vec_C1;

	if (verbose) {
		cout << "y - C1: " << vec_y_minus_C1 << endl << endl;
	}

	mat_ZZ_p mat_M1_T = transpose(*mat_M1);
	vec_ZZ_p vec_y_p = vec_ZZ_p(INIT_SIZE, msg_length);
	ZZ_p d = ZZ_p(0);

	/* 
		这里为了调用 solve API 解方程组，给等式两边同时取了转置
		也就是说原始的式子是 M1 * Y' = Y - C1
		现在我们计算的是    Y' * M1_T  = Y - C1
		因为 solve 这个 API 计算的是 x * A = b
	*/

	solve(d, vec_y_p, mat_M1_T, vec_y_minus_C1);
	// d 是系数矩阵的行列式，如果 d = 0 说明式子无解
	if (d == ZZ_p(0)) {
		cout << "no solution" << endl;
		exit(1);
	}

	if (verbose) {
		cout << "vec_y_p: " << vec_y_p << endl << endl;
	}

	// 接下来对中心多项式 F 求逆，也是求出𝑋 = (𝑥_1,… , 𝑥_n) 的过程
	// 详情见论文正文第15页（共24页）

	vec_ZZ_p vec_x_value = vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);

restart:
	
	{
		// 先用rand()函数随机设置𝑥1, …, 𝑥𝑣1的值
		int x_idx = 0;
		for (int i = 0; i < VINEGAR_VARIABLE_V1; i++) {
			//vec_x_value[i] = rand() % FINITE_FIELD_Q;
			vec_x_value[i] = 1;
			x_idx++;
		}

		// 以层数为计数器
		int k_begin = vinegar[0];  // 索引 k 是为了取 A[k] B[k] 的值参与运算
		int k_end = vinegar[0];
		for (int l = 0; l < LAYERS_U - 1; l++) {
			k_end += oil[l];

			// 先计算式 (3.13) 的各项系数，并定义一个大小为 𝑜_𝑙 × 𝑜_𝑙 的 mat_ZZ_p 矩阵和一个长度为 𝑜_𝑙 的 vec_ZZ_p 的向量

			mat_ZZ_p mat_left = mat_ZZ_p(INIT_SIZE, oil[l], oil[l]);
			vec_ZZ_p vec_right = vec_ZZ_p(INIT_SIZE, oil[l]);

			int mat_left_idx_row = 0;
			int mat_left_idx_col = 0;

			for (int k = k_begin; k < k_end; k++) {
				// 下面这个 for 循环计算式 (3.13) 等号左边式子的系数
				for (int j = vinegar[l]; j < vinegar[l + 1]; j++) {
					ZZ_p mat_sum = (*(vec_B[k]))[j];

					for (int i = 0; i < vinegar[l]; i++) {
						mat_sum += vec_x_value[i] * (*(mat_A[k]))[i][j];
					}

					mat_left[mat_left_idx_row][mat_left_idx_col] = mat_sum;
					mat_left_idx_col++;
				}

				// 下面计算式 (3.13) 等号右边的常数向量
				ZZ_p vec_sum = vec_y_p[k - vinegar[0]];
				for (int i = 0; i < vinegar[l]; i++) {
					for (int j = i; j < vinegar[l]; j++) {
						vec_sum -= ((*(mat_A[k]))[i][j]) * vec_x_value[i] * vec_x_value[j];
					}
				}
				for (int i = 0; i < vinegar[l]; i++) {
					vec_sum -= (*(vec_B[k]))[i] * vec_x_value[i];
				}
				vec_right[mat_left_idx_row] = vec_sum;

				mat_left_idx_row++;
				mat_left_idx_col = 0;
			}

			if (verbose) {
				print_matrix_with_info("mat left - layer index: ", l, &mat_left);
				print_vector_with_info("vec right - layer index", l, &vec_right);
			}

			// solve 函数计算对应油变量𝑥_𝑣𝑖+1,… , 𝑥_𝑣(𝑖+1) 的值
			ZZ_p d = ZZ_p(0);
			vec_ZZ_p result = vec_ZZ_p(INIT_SIZE, oil[l]);
			mat_ZZ_p mat_left_T = transpose(mat_left);

			solve(d, result, mat_left_T, vec_right);
			if (d == ZZ_p(0)) {
				// 若无解，该过程就重新开始
				cout << "there's no solution" << endl;
				goto restart;
			}
			
			for (int i = 0; i < result.length(); i++) {
				vec_x_value[x_idx] = result[i];
				x_idx++;
			}
			k_begin += oil[l];
		}
	}

	if (verbose) {
		cout << "解出的 x: " << vec_x_value << endl;
	}

	// 接下来就是对 L_2 求逆
	vec_ZZ_p vec_x_minus_C2 = vec_x_value - *vec_C2;
	mat_ZZ_p mat_M2_T = transpose(*mat_M2);
	vec_ZZ_p vec_x_p = vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);
	d = ZZ_p(0);

	solve(d, vec_x_p, mat_M2_T, vec_x_minus_C2);
	if (d == ZZ_p(0)) {
		cout << "no solution" << endl;
		exit(1);
	}

	// 解 X' 就是消息 Y 的一个签名
	cout << "Signature(X\'): " << vec_x_p << endl;

	return vec_x_p;
}

// check 方法验证签名的有效性，将消息𝑌的签名𝑋′代入二次多变量方程组𝑃(x)即可
// 见 3.2.6 节 签名验证模块设计
void RainbowSignature::check(vec_ZZ_p signature, vec_ZZ_p msg) {
	vec_ZZ_p result = vec_ZZ_p(INIT_SIZE, mat_A_pp.size());
	int result_idx = 0;
	for (int i = 0; i < mat_A_pp.size(); i++) {
		ZZ_p sum = ZZ_p(0);

		vec_ZZ_p temp_vec = vec_ZZ_p(INIT_SIZE, VARIABLE_COUNT_N);
		mul(temp_vec, signature, *(mat_A_pp[i]));
		ZZ_p temp_num = ZZ_p(0);
		InnerProduct(temp_num, temp_vec, signature);
		sum += temp_num;

		temp_num = ZZ_p(0);
		InnerProduct(temp_num, *(vec_B_pp[i]), signature);
		sum += temp_num;

		sum += num_C_p[i];

		result[result_idx] = sum;
		result_idx++;
	}

	cout << endl << "Original Message:    " << msg << endl;
	cout << endl << "Verification Result: " << result << endl;
}

void print_preset_variables() {
	cout << "Finite Field Q: " << FINITE_FIELD_Q << endl;
	cout << "Layers U: " << LAYERS_U << endl;
	cout << "Variables Count N: " << VARIABLE_COUNT_N << endl;
	cout << endl;
}

void print_matrix_with_info(string name, int index, mat_ZZ_p* mat) {
	cout << name << " index: " << index << " size: " << mat->NumRows() << "×" << mat->NumCols() << endl << *mat << endl << endl;
}

void print_vector_with_info(string name, int index, vec_ZZ_p* vec) {
	cout << name << " index: " << index << " size: " << vec->length() << endl << *vec << endl << endl;
}

int main()
{
	print_preset_variables();
	srand(time(NULL));
	// verbose = true;

	// 在所有含有 ZZ_p 成分的模块中实现有限域上的运算
	ZZ field_number = ZZ(FINITE_FIELD_Q);
	ZZ_p::init(field_number);

	// !---------------测试数据以4.3节为例---------------！

	// 此处填写除 v_1 v_n 之外的变量
	vector<int> rest_v_i = vector<int>{ 12,17,22 };

	RainbowSignature sign = RainbowSignature();
	sign.set_variables(rest_v_i);
	
	for (int i = 0; i < LAYERS_U; i++) {
		cout << "vinegar: " << (sign.vinegar[i]) << endl;
	}
	for (int i = 0; i < LAYERS_U - 1; i++) {
		cout << "oil: " << sign.oil[i] << endl;
	}

	sign.generate_coefficient_mat_vec();
	sign.generate_reversible_affine_transformation();
	sign.generate_public_key();

	// 此处填写要签名的消息
	vector<int> message = vector<int>{
		1,0,5,9,23, 100,201,56,203,88, 178,222,123,234,56, 0,0,0,0,0, 155,86,85,137,22, 0,0 };

	// 把 int 转换成 ZZ_p
	vec_ZZ_p msg = vec_ZZ_p(INIT_SIZE, 27);
	for (int i = 0; i < message.size(); i++) {
		msg[i] = ZZ_p(message[i]);
	}
	vec_ZZ_p sig = sign.sign(msg);
	sign.check(sig, msg);

	return 0;
}
