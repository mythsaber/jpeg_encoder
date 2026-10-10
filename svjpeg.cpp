#include<iostream>
#include<string>
#include<cstdio>
#include<vector>
#include<list>
#include<map>
#include<algorithm>
#include<numeric>
#include<tuple>
#include<utility>
#include<functional>  //std::greater
#include<fstream>
#include<assert.h>
#include<cmath>
#include<iterator>
using std::vector;
using std::cout;
using std::endl;


//输入的frequency不要求有序
vector<int> get_depth_in_huffman_tree(const vector<int>& frequency)
{
	vector<int> res(frequency.size(), 1);
	if (frequency.size() <= 2) return res;

	for (auto& i : res) i = 0;

	struct node
	{
		int index; //index in frequency, -1表示不是叶节点
		int frequency;
		std::list<int> leaves_idx; //该结点下的叶节点在原vector中的索引
		node(int index_, int frequency_,std::list<int> leaves_idx_)
			:index(index_), frequency(frequency_),leaves_idx(leaves_idx_){}
		bool operator <(const node& rhs)
		{
			return frequency < rhs.frequency;
		}
	};
	std::list<node> medium_list;
	for (int i = 0; i < (int)frequency.size(); i++)
	{
		node temp(i, frequency[i],std::list<int>());
		medium_list.push_back(temp);
	}
	medium_list.sort();  //std::less<node>, left<right，因此可提供成员函数，也可非成员函数

	while (medium_list.size()!=1)
	{
		node x = *medium_list.begin();
		medium_list.pop_front();
		node y = *medium_list.begin();
		medium_list.pop_front();
		std::list<int> all_leaves = x.leaves_idx.empty() ? std::list<int>{x.index} : std::move(x.leaves_idx);
		if (y.leaves_idx.empty()) 
		{
			all_leaves.push_back(y.index);
		}
		else
		{
			all_leaves.splice(all_leaves.end(),y.leaves_idx);
		}
		node father(-1,x.frequency+y.frequency,all_leaves);//new node

		//在有序链表中插入一个结点
		auto insert_pos=std::lower_bound(medium_list.begin(), medium_list.end(), father); //std::less<node>, left<right
		medium_list.insert(insert_pos,father);

		for (auto i : father.leaves_idx)
		{
			res[i]++;
		}
	}
	return res;
}

void  test_get_depth_in_huffman_tree()
{
	vector<int> frequency{22,10,11,12,13,23};
	vector<int> depth = get_depth_in_huffman_tree(frequency);
	for (auto i : depth) cout << i << " "; cout << endl;
}

//输出vetcor<pair<编码长度n，编码（16位整数的低n位）>>
vector<std::pair<int,uint16_t>> canonical_huffman_encode(const vector<int>& depth)
{
	if (depth.empty()) return vector<std::pair<int, uint16_t>>();

	vector<int> len(depth); //码长可能被调整，故拷贝一份

	//std::tuple<在初始frequency中的索引、编码长度、编码（16位整数的低n位)
	vector<std::tuple<int, int, uint16_t>> idx_len_code(len.size());

	for (int attempt = 0; ; attempt++)
	{
		for (unsigned i = 0; i < len.size(); i++)
		{
			std::get<0>(idx_len_code[i]) = i; //std::get<>()返回引用
			std::get<1>(idx_len_code[i]) = len[i];
			std::get<2>(idx_len_code[i]) = 0;
		}

		//按编码长度由小到大排列(stable保证等长时的相对次序确定)
		std::stable_sort(idx_len_code.begin(), idx_len_code.end(),
			[](const std::tuple<int, int, uint16_t>& lhs, const std::tuple<int, int, uint16_t>& rhs)->bool
		{
			return std::get<1>(lhs) < std::get<1>(rhs);
		});
		std::get<2>(idx_len_code[0]) = 0;
		for (unsigned i= 1; i <  idx_len_code.size(); i++)
		{
			int steps = std::get<1>(idx_len_code[i]) - std::get<1>(idx_len_code[i-1]);
			std::get<2>(idx_len_code[i]) = (std::get<2>(idx_len_code[i - 1]) + 1) << steps;
		}

		/*
		T.81规定不能使用全1的码字：解码器读到全1只能当作末尾的比特填充，
		否则最后那个字节的填充会被误认成码字。规范编码的码字单调递增，
		所以只有最后一个（也是最大的）码字有可能全为1，把它加长1位后
		重新分配即可(加长只会让Kraft不等式更容易满足)。
		*/
		int last_len = std::get<1>(idx_len_code.back());
		bool all_ones = last_len < 16 &&
			std::get<2>(idx_len_code.back()) == (uint16_t)((1u << last_len) - 1);
		if (!all_ones || attempt >= 16) break;

		len[std::get<0>(idx_len_code.back())]++;
	}

	//返回与初始输入的frequency对应的vetcor<pair<编码长度n，编码（16位整数的低n位>>
	vector<std::pair<int, uint16_t>> res(idx_len_code.size());
	for (unsigned i = 0; i < idx_len_code.size(); i++)
	{
		int index = std::get<0>(idx_len_code[i]);
		res[index] = std::make_pair(std::get<1>(idx_len_code[i]), std::get<2>(idx_len_code[i]));
	}
	return res;
}

void test_canonical_huffman_encode()
{
	vector<int> frequency{ 22,10,11,12,13,23 };
	vector<int> depth = get_depth_in_huffman_tree(frequency);
	vector<std::pair<int, uint16_t>> len_code = canonical_huffman_encode(depth);
	for (auto i : frequency) cout << i << " "; cout << endl;
	for (unsigned i=0; i < frequency.size(); i++)
	{
		int len = std::get<0>(len_code[i]);
		uint16_t code = std::get<1>(len_code[i]);

		char str[17];
		_itoa_s((int)code, str, 17,2);   //转换成字符串，进制基数为2
		printf("len=%d,code=%s\n", len, str);
	}

}

//由frequency[256]得到出现次数不为0的符号的频数
//vector<int>& symbol为符号的值，uint16_t
//vector<int>& symbol_frequency为对应符号的频数
void get_symbol_frequency(int* frequency, int len, vector<int>& symbol,
	vector<int>& symbol_frequency)
{
	assert(len == 256);
	assert(symbol.empty());
	assert(symbol_frequency.empty());

	vector<int> res;
	for (int i = 0; i < len; i++)
	{
		if (frequency[i] != 0)
		{
			symbol.push_back(i);
			symbol_frequency.push_back(frequency[i]);
		}
	}
}

//对[0,255]中的每一个数huffman编码
//出现次数为0的pair<huffman编码长度n为0，huffman编码（16位整数的低n位）为0>
vector<std::pair<int, uint16_t>> get_each_uint8_t_huffman(int* frequency, int len)
{
	assert(len == 256);

	vector<int> symbol;
	vector<int> symbol_frequency;
	get_symbol_frequency(frequency, len, symbol, symbol_frequency);

	vector<int> depth = get_depth_in_huffman_tree(symbol_frequency);
	vector<std::pair<int, uint16_t>> symbol_huffman_code;
	symbol_huffman_code = canonical_huffman_encode(depth);

	vector<std::pair<int, uint16_t>> huffman_len_code(len, std::make_pair(0, 0));
	for (unsigned int i = 0; i < symbol.size(); i++)
	{
		huffman_len_code[symbol[i]] = symbol_huffman_code[i];
	}
	return huffman_len_code;
}

void test_get_each_uint8_t_frequency_write_huffman_table()
{
	vector<int> frequency;
	frequency.insert(frequency.end(), 10, 0);
	frequency.insert(frequency.end(), { 22,10,11,12,13,23 });
	frequency.resize(256);

	vector<std::pair<int, uint16_t>> huffman_table =
		get_each_uint8_t_huffman(&frequency[0], 256);

	for (unsigned i = 0; i < huffman_table.size(); i++)
	{
		int len = huffman_table[i].first;
		uint16_t code = huffman_table[i].second;

		char str[17];
		_itoa_s((int)code, str, 17, 2);   //转换成字符串，进制基数为2
		printf("i=%d: len=%d,code=%s\n", i,len, str);
	}

	//要写入jpeg文件的huffman表相关内容
	auto write_huffman = [](const vector<std::pair<int, uint16_t>>& table, std::ofstream& ofs)
	{
		assert(table.size() == 256);

		//长度为n的码字有几个
		uint8_t nums[16] = {0};
		for (const auto& i : table) nums[i.first]++;

		vector<uint8_t> not_zerolencode_idx;
		vector<std::pair<int, uint16_t>> not_zerolencode_table;
		for (unsigned int i = 0; i < table.size(); i++)
		{
				if (table[i].first != 0)
				{
					not_zerolencode_idx.push_back((uint8_t)i);
					not_zerolencode_table.push_back(table[i]);
				}
		}

		//获取长度不为0的码字由小到大排序后的索引
		vector<uint8_t> idx(not_zerolencode_table.size());
		for (unsigned int i = 0; i < idx.size(); i++) idx[i] = (uint8_t)i;
		std::sort(idx.begin(), idx.end(),
			[&not_zerolencode_table](int lhs, int rhs) {
			return not_zerolencode_table[lhs].second <not_zerolencode_table[rhs].second;
		});

		vector<uint8_t> sorted_not_zerolencode_idx(not_zerolencode_idx.size());
		for (unsigned int i = 0; i < sorted_not_zerolencode_idx.size(); i++)
		{
			sorted_not_zerolencode_idx[i] = not_zerolencode_idx[idx[i]];
		}

		ofs.write((char*)nums, 16);
		ofs.write((char*)&sorted_not_zerolencode_idx[0], sorted_not_zerolencode_idx.size());
	};
	
	std::string path = "C:/Users/myth/Desktop/";
	std::ofstream ofs(path+"test_get_each_uint8_t_frequency_write_huffman_table",
		std::ios::binary);
	if (!ofs)
	{
		cout << "Couldn't open txt" << endl;
		return;
	}
	cout << "0-----------------------------------------------0" << endl;
	write_huffman(huffman_table, ofs);
}

//将RGB24格式像素数据转换为YUV420P格式像素数据
//对应的采样因子为：
//Y：水平采用因子=2； 垂直采样因子=2；
//V：水平采用因子  = 1；垂直采样因子 = 1
//U：水平采用因子 = 1；垂直采样因子 = 1
void rgb24_to_yuv420P(const unsigned char *rgbbuf, int w, int h, unsigned char *yuvbuf)
{
	unsigned char* ptrY = yuvbuf;
	int iy = 0;
	unsigned char* ptrU = yuvbuf + w * h; //ptrV
	int iu = 0;
	unsigned char* ptrV = yuvbuf + w * h * 5 / 4; //ptrU
	int iv = 0;
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			unsigned char r = rgbbuf[3 * (i*w + j)];
			unsigned char g = rgbbuf[3 * (i*w + j) + 1];
			unsigned char b = rgbbuf[3 * (i*w + j) + 2];
			unsigned char y = (unsigned char)((66 * r + 129 * g + 25 * b + 128) >> 8) + 16;
			unsigned char u = (unsigned char)((-38 * r - 74 * g + 112 * b + 128) >> 8) + 128;
			unsigned char v = (unsigned char)((112 * r - 94 * g - 18 * b + 128) >> 8) + 128;
			//转换后Y的范围是[16,235]，UV 的范围是[16,239]，跟[0,255]很接近

			ptrY[iy++] = y;
			if ((i & 1) == 0 && (j & 1) == 0)  //必须加()，注意优先级
			{
				ptrU[iu++] = u;
			}
			else if ((i & 1) != 0 && (j & 1) == 0)
			{
				ptrV[iv++] = v;
			}
		}
	}
}

//将RGB24格式像素数据转换为YUV411P格式像素数据,w*h*1.5=w*h+w*h/4+w*h/4
void rgb24_to_yuv411P(const unsigned char *rgbbuf, int w, int h,unsigned char *yuvbuf)
{
	unsigned char* ptrY = yuvbuf;
	int iy = 0;
	unsigned char* ptrU = yuvbuf + w*h;
	int iu = 0;
	unsigned char* ptrV = yuvbuf + w*h * 5 / 4;
	int iv = 0;
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++)
		{
			unsigned char r = rgbbuf[3 * (i*w + j)];
			unsigned char g = rgbbuf[3 * (i*w + j) + 1];
			unsigned char b = rgbbuf[3 * (i*w + j) + 2];
			unsigned char y = (unsigned char)((66 * r + 129 * g + 25 * b + 128) >> 8) + 16;
			unsigned char u = (unsigned char)((-38 * r - 74 * g + 112 * b + 128) >> 8) + 128;
			unsigned char v = (unsigned char)((112 * r - 94 * g - 18 * b + 128) >> 8) + 128;
			//转换后Y的范围是[16,235]，UV 的范围是[16,239]，跟[0,255]很接近

			ptrY[iy++] = y;
			if (j % 4 == 0)
			{
				ptrU[iu++] = u;
				ptrV[iv++] = v;
			}
		}
	}
}

void test_rgb24_to_yuv411P()
{
	const char* src = "colorbar_640x360.rgb";
	vector<unsigned char> rgbbuf(640 * 360 * 3);
	std::ifstream ifs(src, std::ios::binary);
	if (!ifs) { cout << "not open " << src << endl;  return; }
	ifs.read((char*)&rgbbuf[0], 640 * 360 * 3);
	vector<unsigned char>  yuvbuf(640 * 360 * 3 / 2);
	rgb24_to_yuv411P(&rgbbuf[0], 640, 360, &yuvbuf[0]);
	
	const char* dest = "colorbar_yuv411P_640x360.yuv";
	std::ofstream ofs(dest, std::ios::binary);
	if (!ifs) { cout << "not open " << dest << endl;  return; }
	ofs.write((char*)&yuvbuf[0], 640 * 360 * 3 / 2);
}

//将图像的长和宽都补为n的倍数
void padto(vector<uint8_t>& data, int& width, int& height, int n)
{
	assert((n == 8) || (n == 16) | (n == 32));
	int addw = width%n==0?0:n - width%n;
	int addh = height%n == 0 ? 0 : n - height%n;
	if (addw == 0 && addh == 0) return;
	int neww = width + addw;
	int newh = height + addh;
	vector<uint8_t> res_img(neww*newh * 3);
	for (int i = 0; i < height; i++)
	{
		memcpy(&res_img[i*neww*3],&data[i*width*3],width*3);
		for (int k = 0; k < addw; k++)
		{
			memcpy(&res_img[i*neww * 3 + width * 3+3*k], &res_img[i*neww * 3 + width * 3 - 3], 3);
		}
	}
	for (int k = 0; k < addh; k++)
	{
		memcpy(&res_img[(height+k)*neww * 3],&res_img[(height-1)*neww*3],neww*3);
	}

	data = std::move(res_img);
	width = neww;
	height = newh;
}

void test_padto()
{
	vector<uint8_t> data{ 1,1,1, 2,2,2, 3,3,3, 4,4,4,
										   5,5,5, 6,6,6, 7,7,8, 9,10,11,
										  12,13,14, 25,26,27, 38,37,36, 40,41,42 };
	int width = 4;
	int height = 3;
	const int n =8;
	padto(data, width, height, n);
	for (unsigned int i = 0; i < data.size(); i++)
	{
		cout <<(int)data[i] << " ";
		if (i%(n*3) == (n*3-1)) cout << endl;
	}
	cout << endl;

	padto(data, width, height, 5);
}

//DCT 变化需要定义域对称，在做 DCT 变化之前要先将矩阵中的减小
//左移 128，使得值域落在 -128~127
vector<int16_t> less_n_each_ele(const vector<uint8_t> val,int n)
{
	vector<int16_t> res(val.size());
	for (unsigned int i = 0; i < val.size(); i++)
	{
		res[i] = val[i] - n;
	}
	return res;
}

//对一个8*8的数组DCT变换
vector<int16_t> dct(const vector<int16_t>& input)
{
	const int N = (int)std::sqrt(input.size());
	assert( N== 8);
	vector<double> res_dbl(input.size());

	vector<vector<double>> dkj(N);
	for (int i = 0; i < N; i++)
	{
		dkj[i] = vector<double>(N,1);
	}
	dkj[0][0] = 0.5;
	for (int i = 1; i < N; i++) dkj[0][i] = 1 / std::sqrt(2);
	for(int i = 1; i < N; i++) dkj[i][0]= 1 / std::sqrt(2);

	constexpr double PI = 3.141592653;
	for (int k = 0; k < N; k++)
	{
		for (int j = 0; j < N; j++)
		{
			res_dbl[k*N + j] = dkj[k][j] * 2.0 / N;
			double sum = 0;
			for (int m = 0; m < N; m++)
			{
				for (int n = 0; n < N; n++)
				{
					sum+=input[m*N + n] *
						std::cos(PI*(2 * m + 1)*k / 2 / N)*
						std::cos(PI*(2 * n + 1)*j / 2 / N);
				}
			}
			res_dbl[k*N + j] *= sum;
		}
	}

	vector<int16_t> res(input.size());
	for (unsigned int i = 0; i < res.size(); i++)
	{
		res[i] = res_dbl[i] > 0 ? (int16_t)(res_dbl[i] + 0.5) : (int16_t)(res_dbl[i] - 0.5);
	}
	return res;
}

void test_dct()
{
	vector<int16_t> input_int{ 0, 1, 5, 6,14,15,27,28,
		2, 4, 7,13,16,26,29,42,
		3, 8,12,17,25,30,41,43,
		9,11,18,24,31,40,44,53,
		10,19,23,32,39,45,52,54,
		20,22,33,38,46,51,55,60,
		21,34,37,47,50,56,59,61,
		35,36,48,49,57,58,62,63 };
	vector<double> input_dbl(input_int.size());
	for (unsigned int i = 0; i < input_dbl.size(); i++)
	{
		input_dbl[i] = input_int[i];
	}

	auto dct_another=[](const double* input, std::int16_t* output) {
		constexpr double PI = 3.1415926;
		auto f = [](int x) { return x ? 1 : 1 / std::sqrt(2); };
		for (int v = 0; v < 8; ++v) for (int u = 0; u < 8; ++u) {
			double sum = 0;
			for (int y = 0; y < 8; ++y) for (int x = 0; x < 8; ++x)
				sum += input[y * 8 + x] * std::cos((2 * x + 1) * u * PI * 0.0625) * std::cos((2 * y + 1) * v * PI * 0.0625);
			output[v * 8 + u] = (uint16_t)std::round(sum * 0.25 * f(u) * f(v));
		}
	};
	int16_t output_another[64];
	dct_another(&input_dbl[0],&output_another[0]);

	vector<int16_t> output_mine = dct(input_int);

	for (int i = 0; i < 64; i++)
	{
		cout << output_another[i] << " ";
		if (i % 8 == 7) cout << endl;
	}
	cout <<"---------------------------------------"<< endl;
	for (int i = 0; i < 64; i++)
	{
		cout << output_mine[i] << " ";
		if (i % 8 == 7) cout << endl;
	}

}

//zigzag 排列，将相似频率组在一起
void zigzag_sort(vector<int16_t>&  x)
{
	assert(x.size() == 64);
	int src_table[64]={
		0, 1, 5, 6,14,15,27,28,
		2, 4, 7,13,16,26,29,42,
		3, 8,12,17,25,30,41,43,
		9,11,18,24,31,40,44,53,
		10,19,23,32,39,45,52,54,
		20,22,33,38,46,51,55,60,
		21,34,37,47,50,56,59,61,
		35,36,48,49,57,58,62,63
	};
	int dst_table[64];
	for (int i = 0; i < 64; i++)
	{
		dst_table[src_table[i]] = i;
	}

	vector<int16_t>  y(x.size());
	for (unsigned int i = 0; i < y.size(); i++)
	{
		y[i] = x[dst_table[i]];
	}

	x = y;
}

void test_zigzag_sort()
{
	int src_table[64] = {
		0, 1, 5, 6,14,15,27,28,
		2, 4, 7,13,16,26,29,42,
		3, 8,12,17,25,30,41,43,
		9,11,18,24,31,40,44,53,
		10,19,23,32,39,45,52,54,
		20,22,33,38,46,51,55,60,
		21,34,37,47,50,56,59,61,
		35,36,48,49,57,58,62,63
	};
	int dst_table[64];
	for (int i = 0; i < 64; i++)
	{
		dst_table[src_table[i]] = i;
	}
	for (int i = 0; i < 64; i++)
	{
		cout << dst_table[i] << " ";
		if (i % 8 == 7) cout << endl;
	}
}

//量化，除以量化表里对应值并四舍五入
void quantify(vector<int16_t>&  val, vector<uint8_t>&  table)
{
	assert(val.size() == table.size());

	//量化后的数据有可能超过 2 byte 有符号整数的处理范围.
	auto clip = [](int32_t x) ->int32_t
	{
		if (x > std::numeric_limits<int16_t>::max())
			return std::numeric_limits<int16_t>::max();
		else if (x < std::numeric_limits<int16_t>::min())
			return std::numeric_limits<int16_t>::min();
		else
			return x;
	};

	for (unsigned  i = 0; i < val.size(); i++)
	{
		double quotient = ((double)val[i]) / table[i];
		int32_t temp = quotient > 0 ? (int32_t)(quotient + 0.5) : (int32_t)(quotient - 0.5);
		val[i] = (int16_t)clip(temp);
	}
}

void test_quantify()
{
	cout << std::numeric_limits<int16_t>::max() << endl; //32767
	cout << std::numeric_limits<int16_t>::min() << endl; // -32768

	double quotient = 3.6;
	int32_t temp = quotient > 0 ? (int32_t)(quotient + 0.5) : (int32_t)(quotient - 0.5);
	cout << temp << endl;;

	quotient = -3.6;
	temp = quotient > 0 ? (int32_t)(quotient + 0.5) : (int32_t)(quotient - 0.5);
	cout << temp << endl;;
}

//对AC分量RLE编码,每组数字前一个表示 0 的数量的必须是 4 bit, 即 0~15
//vector<std::pair<前续0数，当前值>>
vector<std::pair<uint8_t,int16_t>> rle_encode(const vector<int16_t> input,unsigned first_bitnums)
{
	if (input.empty()) return vector<std::pair<uint8_t, int16_t>>();

	assert(first_bitnums >= 1 && first_bitnums <= 8);
	const int first_max = (int)std::pow(2, first_bitnums) - 1;
	vector<std::pair<uint8_t, int16_t>> res;

	unsigned all_zero_begin_idx= input.size();//由该idx开始，全为0
	for (int i = (int)input.size()-1; i >= 0; i--)
	{
		if (input[i] == 0) all_zero_begin_idx = i;
		else
		{
			break;
		}
	}
	uint8_t consecutive = 0;
	for (unsigned i = 0; i < all_zero_begin_idx; i++)
	{
		if (input[i] != 0)
		{
			res.push_back(std::make_pair(consecutive,input[i]));
			consecutive = 0;
		}
		else
		{
			consecutive++;
			if (consecutive == first_max+1)
			{
				res.push_back(std::make_pair(first_max, 0));
				consecutive = 0;
			}
		}
	}
	//EOB 是一个结束标记, 表示后面都是 0 了. 实际上, 我们用 (0,0) 表示 EOB
	//但是, 如果这组数字不以 0 结束, 那么就不需要 EOB
	if (all_zero_begin_idx != input.size())
	{
		res.push_back(std::make_pair(0, 0));
	}
	return res;
}
void test_rle_encode()
{
	vector<int16_t> a = { 57,45,0,0,0,0,23,0,-30,-16,0,0,1,0,0,0, 0 , 0 ,0 , 0,0,0,0 };
	vector<std::pair<uint8_t, int16_t>> r1 = rle_encode(a,4);
	for (auto& i : r1) { cout << "(" << (int)i.first << "," << i.second << ")  "; } cout << endl;

	vector<int16_t> b;
	b.push_back(57);
	b.insert(b.end(),18,0);
	b.push_back(3);
	b.insert(b.end(),4,0);
	b.push_back(2);
	b.insert(b.end(),33,0);
	b.push_back(895);
	b.insert(b.end(), 15, 0);
	b.push_back(777);
	b.insert(b.end(),100,0);
	vector<std::pair<uint8_t, int16_t>> r2 = rle_encode(b, 4);
	for (auto& i : r2) { cout << "(" << (int)i.first << "," << i.second << ")  "; } cout << endl;
}

//对一副图像的Y直流分量差分编码DPCM
vector<int16_t> y_dpcm_encode(const vector<vector<int16_t>>& blocks)
{
	/*
	3个颜色分量的直流变量是分开进行差分编码的
	Diff为差分校正变量，也就是直接解码出来的直流系数。
	但如果当前颜色分量单元是第一个单元，则解码出来的
	直流数值就是真正的直流变量
	*/
	assert(!blocks.empty());
	vector<int16_t> code(blocks.size());
	code[0] = blocks[0][0];
	for (unsigned int i = 1; i < blocks.size(); i++)
	{
		code[i] = blocks[i][0] - blocks[i - 1][0];
	}
	return code;
}

//对一副图像的UV直流分量差分编码DPCM
vector<int16_t> uv_dpcm_encode(const vector<vector<int16_t>>& blocks)
{
	/*
	3个颜色分量的直流变量是分开进行差分编码的
	Diff为差分校正变量，也就是直接解码出来的直流系数。
	但如果当前颜色分量单元是第一个单元，则解码出来的
	直流数值就是真正的直流变量
	*/
	assert(!blocks.empty());
	vector<int16_t> code(blocks.size());
	//u
	code[0] = blocks[0][0];
	for (unsigned int i = 1; i < blocks.size()/2; i++)
	{
		code[i] = blocks[i][0] - blocks[i - 1][0];
	}
	//v
	code[blocks.size() / 2] = blocks[blocks.size() / 2][0];
	for (unsigned int i = blocks.size() / 2+1; i < blocks.size(); i++)
	{
		code[i] = blocks[i][0] - blocks[i - 1][0];
	}
	return code;
}

//为了提高储存效率, JPEG 里并不直接保存数值, 而是将数值按位数分成 16 组
//std::pair<长度bit数,按分组编码后的值>
std::pair<uint8_t,uint16_t> group_encode(int16_t src) 
{
	std::pair<uint8_t, uint16_t> res;
	if (src == 0)
	{
		res.first = 0;
		res.second = 0;
		return res;
	}
	res.first = (uint8_t)log2(abs(src)) + 1;
	if (src > 0) res.second = src;
	if (src < 0)
	{
		int16_t abx_small = (int16_t)pow(2, res.first)-1;
		res.second = abx_small + src;
	}
	return res;
}
void test_group_encode()
{
	cout << "("<< (int)group_encode(0).first<<","<< group_encode(0).second<<")" << endl;
	for (int layer = 1; layer < 5; layer++)
	{
		int16_t abx_small = (int16_t)pow(2, layer - 1);
		int16_t abx_big = (int16_t)pow(2, layer )-1;
		for (int16_t i = -1 * abx_big; i <= -1 * abx_small;i++)
		{
			cout << "(" << (int)group_encode(i).first << "," << group_encode(i).second << ") " ;
		}
		cout << "  |  ";
		for (int16_t i =abx_small; i <= abx_big;i++)
		{
			cout << "(" << (int)group_encode(i).first << "," << group_encode(i).second << ")";
		}
		cout << endl;
	}
}

//两个uint8_t合为一个uint，高 4 位是前续 0 的个数, 低 4 位描述了后面数字的位数
uint8_t merge_to_symbol(int8_t zeronum,int8_t len)
{
	return (zeronum << 4)+len;
}

//vector<std::pair<uint8_t符号，对值分组编码后的bit数，按分组编码后的值>>
vector<std::tuple<uint8_t, uint8_t, uint16_t>> get_tuple(
	const vector<std::pair<uint8_t, int16_t>>& rle_code)
{
	vector<std::tuple<uint8_t, uint8_t, uint16_t>> tu_vec(rle_code.size());
	for (unsigned i = 0; i < rle_code.size(); i++)
	{
		uint8_t zeronum = rle_code[i].first;
		std::pair<uint8_t, uint16_t> group_code=group_encode(rle_code[i].second);
		uint8_t symbol = merge_to_symbol(zeronum,group_code.first);
		std::get<0>(tu_vec[i]) = symbol;//前续0数及长度组成的uint8_t
		std::get<1>(tu_vec[i]) = group_code.first;//对值分组编码后的bit数
		std::get<2>(tu_vec[i]) = group_code.second;//按分组编码后的值
	}
	return tu_vec;
}

//frequency[256],符号i出现的次数为frequency[i]
void get_each_uint8_t_frequency(int* frequency, int len,
	const vector<std::tuple<uint8_t, uint8_t, int16_t>>& tu_vec)
{
	assert(len == 256);
	for (const auto& i : tu_vec)
	{
		frequency[std::get<0>(i)]++;
	}
}

vector<vector<uint8_t>> get_y_blocks(const vector<uint8_t>& y_vals,int width,int height)
{
	int row = height / 8;//row*col个8*8的块
	int col = width / 8;
	vector<vector<uint8_t>> y_blocks(row*col);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			vector<uint8_t> block(8 * 8);
			for (int k = 0; k < 8; k++)
			{
				const uint8_t* src = &y_vals[(i * 8)*width + j * 8];
				memcpy(&block[k * 8], src+k*width, 8);
			}
			y_blocks[i*col + j]=block;
		}
	}
	return y_blocks;
}
vector<vector<uint8_t>> get_uv_blocks(const vector<uint8_t>& uv_vals, int width, int height)
{
	assert(uv_vals.size() == width*height / 2);
	vector<uint8_t> u_vals(uv_vals.begin(),uv_vals.begin()+height*width/4);
	vector<uint8_t> v_vals(uv_vals.begin() + height*width / 4, uv_vals.end());

	constexpr int sh = 2; //水平方向每sh个像素点采样一个u或v值
	constexpr int sv = 2; //垂直方向每sv个像素点采样一个u或v值
	const int row = height / sh / 8; //有row*col个block
	const int col = width / sv / 8;

	auto func = [row,col](const vector<uint8_t>& u_or_v_vals, int width, int height)
		->vector<vector<uint8_t>>
	{
		vector<vector<uint8_t>> u_or_v_blocks(row*col);
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				vector<uint8_t> block(8 * 8);
				//当前block首元素在容器u_or_v_vals中的索引
				for (int k = 0; k < 8; k++)
				{
					const uint8_t* src = &u_or_v_vals[i*col * 64 + j * 8];
					memcpy(&block[k * 8], src+k * (col*8), 8);
				}
				u_or_v_blocks[i*col + j] = block;
			}
		}
		return u_or_v_blocks;
	};

	vector<vector<uint8_t>> u_blocks = func(u_vals,width,height);
	vector<vector<uint8_t>> v_blocks = func(v_vals, width, height);
	vector<vector<uint8_t>> uv_blocks(u_blocks.size()+v_blocks.size());
	for (unsigned int i = 0; i < (uv_blocks.size() / 2);i++)
	{
		uv_blocks[i] = u_blocks[i];
	}
	for (unsigned int i = uv_blocks.size() / 2; i < uv_blocks.size(); i++)
	{
		uv_blocks[i] = v_blocks[i- uv_blocks.size() / 2];
	}
	return uv_blocks;
}

void test_get_y_blocks_get_uv_blocks()
{
	int width = 64;
	int height = 24;
	
	vector<uint8_t> y_vals(width*height);
	for (unsigned int i = 0; i < y_vals.size(); i++) y_vals[i] = (3 * i) % 128;
	vector<uint8_t> uv_vals(width*height/2);
	for (unsigned int i = 0; i < uv_vals.size() / 2; i++) uv_vals[i] = (i * 12 + 1) % 128;
	for (unsigned int i = uv_vals.size() / 2; i < uv_vals.size(); i++)
	{
		uv_vals[i] = ((i- uv_vals.size() / 2) * 12 + 2) % 128;
	}

	vector<vector<uint8_t>> y_blocks = get_y_blocks(y_vals,width,height);
	vector<vector<uint8_t>> uv_blocks = get_uv_blocks(uv_vals,width,height);

	std::ofstream ofs_src("C:/Users/myth/Desktop/src_temp.txt");
	cout << "output file:" << "C:/Users/myth/Desktop/src_temp.txt" << endl;
	for (int i = 0; i < width*height*3; i++)
	{
		ofs_src << i%128 << " ";
		if (i % (width * 3) == (width * 3 - 1)) ofs_src << endl;
	}

	std::ofstream  ofs_dst("C:/Users/myth/Desktop/dst_temp.txt");
	cout << "output file:" << "C:/Users/myth/Desktop/dst_temp.txt" << endl;
	ofs_dst << "------------------------------y--------------------------------" << endl;
	int row = height / 8;
	int col = width / 8;
	for (unsigned int i = 0; i < y_blocks.size(); i++)
	{
		assert(y_blocks[i].size() == 64);
		for (auto i : y_blocks[i]) ofs_dst << (int)i << " "; ofs_dst << endl;
	}
	ofs_dst << endl;
	ofs_dst << "------------------------------u--------------------------------" << endl;
	for (unsigned int i = 0; i < uv_blocks.size()/2; i++)
	{
		assert(uv_blocks[i].size() == 64);
		for (auto i : uv_blocks[i]) ofs_dst << (int)i << " "; ofs_dst << endl;
	}
	ofs_dst << endl;
	ofs_dst << "------------------------------v--------------------------------" << endl;
	for (unsigned int i = uv_blocks.size() / 2; i < uv_blocks.size(); i++)
	{
		assert(uv_blocks[i].size() == 64);
		for (auto i : uv_blocks[i]) ofs_dst << (int)i << " "; ofs_dst << endl;
	}
}

void avjpeg(vector<uint8_t>& rgbdata, int width, int height, const char* dst_path)
{
	padto(rgbdata, width, height, 16);
	vector<uint8_t> yvubuf(width*height*3/2);  //先存放Y再存放V再存放U
	rgb24_to_yuv420P(&rgbdata[0], width, height, &yvubuf[0]);

	if (false)
	{
		std::string yuv_path = "C:/Users/myth/Desktop/yuv420p_" + std::to_string(width) + "x" +
		std::to_string(height) + ".yuv";
		std::ofstream ofs_yuv(yuv_path, std::iostream::binary);
		if (!ofs_yuv) { cout << "Couldn't open" << yuv_path << endl; return; }
		ofs_yuv.write((char*)&yvubuf[0],yvubuf.size());
	}

	vector<uint8_t> y_vals(yvubuf.begin(), yvubuf.begin()+ width*height);
	//vector vu_vals的前一半元素为v值，后一半元素为u值
	vector<uint8_t> vu_vals(yvubuf.begin() + width*height, yvubuf.end());
	assert(y_vals.size() == width*height);
	assert(vu_vals.size() == width*height / 2);
	vector<vector<uint8_t>> y_blocks = get_y_blocks(y_vals,width,height);
	vector<vector<uint8_t>> uv_blocks = get_uv_blocks(vu_vals, width, height);

	vector<vector<int16_t>> y_blocks_less_128(y_blocks.size());
	for (unsigned int i = 0; i < y_blocks_less_128.size(); i++)
	{
		y_blocks_less_128[i] = less_n_each_ele(y_blocks[i],128);
	}
	vector<vector<int16_t>> uv_blocks_less_128(uv_blocks.size());
	for (unsigned  i = 0; i < uv_blocks_less_128.size(); i++)
	{
		uv_blocks_less_128[i] = less_n_each_ele(uv_blocks[i], 128);
	}

	if (true)
	{
		cout << endl << endl;
		auto show_block = [](const vector<int16_t>& block)
		{
			for (int i = 0; i < 64; i++)
			{
				printf("%5d ", (int)block[i]);
				if (i % 8 == 7) cout << endl;
			}
		};
		for (const auto& i : y_blocks_less_128) 
		{
			show_block(i); 
			cout << endl;
		}
		cout << "---------------------------------------u、v:------------------------------" << endl;
		for (const auto& i : uv_blocks_less_128)
		{
			show_block(i);
			cout << endl;
		}
	}

	vector<vector<int16_t>> y_blocks_dct(y_blocks_less_128.size());
	for (unsigned int i = 0; i <  y_blocks_dct.size(); i++)
	{
		y_blocks_dct[i] = dct(y_blocks_less_128[i]);
	}
	vector<vector<int16_t>> uv_blocks_dct(uv_blocks_less_128.size());
	for (unsigned int i = 0; i < uv_blocks_dct.size(); i++)
	{
		uv_blocks_dct[i] = dct(uv_blocks_less_128[i]);
	}

	if (true)
	{
		cout << endl << endl;
		auto show_block = [](const vector<int16_t>& block)
		{
			for (int i = 0; i < 64; i++)
			{
				printf("%5d ", (int)block[i]);
				if (i % 8 == 7) cout << endl;
			}
		};
		cout << "---------------------------------------y blocks_dct:------------------------------" << endl;
		for (const auto& i : y_blocks_dct)
		{
			show_block(i);
			cout << endl;
		}
		cout << "---------------------------------------u、v blocks_dct:------------------------------" << endl;
		for (const auto& i : uv_blocks_dct)
		{
			show_block(i);
			cout << endl;
		}
	}

	//zigzag
	for (unsigned int i = 0; i < y_blocks_dct.size(); i++)
	{
		zigzag_sort(y_blocks_dct[i]);
	}
	for (unsigned int i = 0; i < uv_blocks_dct.size(); i++)
	{
		zigzag_sort(uv_blocks_dct[i]);
	}
	if (true)
	{
		cout << endl << endl;
		auto show_block = [](const vector<int16_t>& block)
		{
			for (int i = 0; i < 64; i++)
			{
				printf("%5d ", (int)block[i]);
				if (i % 8 == 7) cout << endl;
			}
		};
		cout << "---------------------------------------y blocks_dct after zigzag:----------------------" << endl;
		for (const auto& i : y_blocks_dct)
		{
			show_block(i);
			cout << endl;
		}
		cout << "---------------------------------------u、v blocks_dct after zigzag:-------------------------" << endl;
		for (const auto& i : uv_blocks_dct)
		{
			show_block(i);
			cout << endl;
		}
	}

	//量化

#if 0
	vector<uint8_t> q_table_y{
		0x02, 0x01, 0x01, 0x01, 0x01, 0x01, 0x02, 0x01,
		0x01, 0x01, 0x02, 0x02, 0x02, 0x02, 0x02, 0x04,
		0x03, 0x02, 0x02, 0x02, 0x02, 0x05, 0x04, 0x04,
		0x03, 0x04, 0x06, 0x05, 0x06, 0x06, 0x06, 0x05,
		0x06, 0x06, 0x06, 0x07, 0x09, 0x08, 0x06, 0x07,
		0x09, 0x07, 0x06, 0x06, 0x08, 0x0B, 0x08, 0x09,
		0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x06, 0x08, 0x0B,
		0x0C, 0x0B, 0x0A, 0x0C, 0x09, 0x0A, 0x0A, 0x0A
	};
	vector<uint8_t> q_table_uv{
		0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x05, 0x03,
		0x03, 0x05, 0x0A, 0x07, 0x06, 0x07, 0x0A, 0x0A,
		0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A,
		0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A,
		0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A,
		0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A,
		0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A,
		0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A, 0x0A
	};
#endif
	vector<uint8_t> q_table_y{ 
		16,11,10,16,24,40,51,61,
		12,12,14,19,26,58,60,55,
		14,13,16,24,40,57,69,56,
		14,17,22,29,51,87,80,62,
		18,22,37,56,68,109,103,77,
		24,35,55,64,81,104,113,92,
		49,64,78,87,103,121,120,101,
		72,92,95,98,112,100,103,99 };
	vector<uint8_t> q_table_uv{
		17,18,24,47,99,99,99,99,
		18,21,26,66,99,99,99,99,
		24,26,56,99,99,99,99,99,
		47,66,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99,
		99,99,99,99,99,99,99,99 };
	for (unsigned int i = 0; i < y_blocks_dct.size(); i++)
	{
		quantify(y_blocks_dct[i], q_table_y);
	}
	for (unsigned int i = 0; i < uv_blocks_dct.size(); i++)
	{
		quantify(uv_blocks_dct[i], q_table_uv);
	}
	if (true)
	{
		cout << endl << endl;
		auto show_block = [](const vector<int16_t>& block)
		{
			for (int i = 0; i < 64; i++)
			{
				printf("%5d ", (int)block[i]);
				if (i % 8 == 7) cout << endl;
			}
		};
		cout << "---------------------------------------y blocks_dct after quan:----------------------" << endl;
		for (const auto& i : y_blocks_dct)
		{
			show_block(i);
			cout << endl;
		}
		cout << "---------------------------------------u、v blocks_dct after quan:-------------------------" << endl;
		for (const auto& i : uv_blocks_dct)
		{
			show_block(i);
			cout << endl;
		}
	}

	//对交流分量RLE(游程编码)
	vector<vector<std::pair<uint8_t, int16_t>>> y_ac_rle_code(y_blocks_dct.size());
	for (unsigned int i = 0; i < y_ac_rle_code.size(); i++)
	{
		y_ac_rle_code[i] = rle_encode(
			std::vector<int16_t>(y_blocks_dct[i].begin() + 1, y_blocks_dct[i].end()), 4);
	}
	vector<vector<std::pair<uint8_t, int16_t>>> uv_ac_rle_code(uv_blocks_dct.size());
	for (unsigned int i = 0; i < uv_ac_rle_code.size(); i++)
	{
		uv_ac_rle_code[i] = rle_encode(
			std::vector<int16_t>(uv_blocks_dct[i].begin() + 1, uv_blocks_dct[i].end()), 4);
	}
	if (true)
	{
		cout << endl << endl;
		auto show_block_rle = [](const vector<std::pair<uint8_t, int16_t>> block_rle)
		{
			for (const auto& p : block_rle)
			{
				printf("(%5d, %5d)   ", (int)p.first, (int)p.second);
			}
			cout << endl;
		};
		cout << "---------------------------------------y 交流 游程编码rle:----------------------" << endl;
		for (const auto& i : y_ac_rle_code)
		{
			show_block_rle(i);
			cout << endl;
		}
		cout << "---------------------------------------uv  交流 游程编码rle:--------------------" << endl;
		for (const auto& i : uv_ac_rle_code)
		{
			show_block_rle(i);
			cout << endl;
		}
	}

	//对直流分量差分编码
	vector<int16_t> y_dc_dpcm_code = y_dpcm_encode(y_blocks_dct);
	vector<int16_t> uv_dc_dpcm_code = uv_dpcm_encode(uv_blocks_dct);
	if (true)
	{
		cout << endl << endl;
		auto show_block_dpcm = [](const vector<int16_t> block_dpcm)
		{
			for (const auto& p : block_dpcm)
			{
				printf("%5d,   ", (int)p);
			}
			cout << endl;
		};
		cout << "---------------------------------------y 直流 差分编码:----------------------" << endl;
		show_block_dpcm(y_dc_dpcm_code);
		cout << endl;
		cout << "--------------------------------------uv  直流 差分编码:--------------------" << endl;
		show_block_dpcm(uv_dc_dpcm_code);
		cout << endl << endl;
	}

	//对直流分量分组，std::pair<长度bit数,按分组编码后的值>
	vector<std::pair<uint8_t, uint16_t>> y_dc_group_code(y_dc_dpcm_code.size());
	for (unsigned int i = 0; i < y_dc_group_code.size(); i++)
	{
		y_dc_group_code[i] = group_encode(y_dc_dpcm_code[i]);
	}
	vector<std::pair<uint8_t, uint16_t>> uv_dc_group_code(uv_dc_dpcm_code.size());
	for (unsigned int i = 0; i < uv_dc_group_code.size(); i++)
	{
		uv_dc_group_code[i] = group_encode(uv_dc_dpcm_code[i]);
	}

	//一张图片Y dc分量分组编码后码字长度的频数
	const int potential_symbol_num = 256;
	int y_dc_frequency[potential_symbol_num] = { 0 };//所有可能的码字长度为0、1、2……15
	for (const auto& i : y_dc_group_code)
	{
		y_dc_frequency[i.first]++;
	}
	//Y dc分量的符号数
	const int y_dc_symbol_num = std::accumulate(y_dc_frequency, 
		y_dc_frequency + potential_symbol_num, 0,
		[](int own, int cur) {return own + (cur == 0 ? 0 : 1); });

	//一张图片UV dc分量分组编码后码字长度的频数
	int uv_dc_frequency[potential_symbol_num] = { 0 };
	for (const auto& i : uv_dc_group_code)
	{
		uv_dc_frequency[i.first]++;
	}
	//UV dc分量的符号数
	const int uv_dc_symbol_num = std::accumulate(uv_dc_frequency, uv_dc_frequency + potential_symbol_num, 0,
		[](int own, int cur) {return own + (cur == 0 ? 0 : 1); });

	//为一张图片生成Y分量、UV分量直流huffman表，用于由uint8_t查看其huffman码字
	//对[0,255]中的每一个数huffman编码
	//出现次数为0的pair<huffman编码长度n为0，huffman编码（16位整数的低n位）为0>
	vector<std::pair<int, uint16_t>> y_dc_each_uint8_t_huffman_code =
		get_each_uint8_t_huffman(y_dc_frequency, potential_symbol_num);
	if (true)
	{
		for (int i = 0; i<potential_symbol_num; i++)
		{
			if (y_dc_frequency[i] != 0) printf("y_dc_frequency[%d] =%d  ", i, y_dc_frequency[i]);
		}
		cout << endl;

		for (const auto& i : y_dc_each_uint8_t_huffman_code)
		{
			if (i.first != 0) printf("(编码长度，码字)：(%d,%d)\n", i.first, i.second);
		}
		cout << endl;
	}
	vector<std::pair<int, uint16_t>> uv_dc_each_uint8_t_huffman_code =
		get_each_uint8_t_huffman(uv_dc_frequency, potential_symbol_num);


	//由每个block生成一个vector<std::tuple<uint8_t符号，对值分组编码后的bit数，按分组编码后的值>>
	vector<vector<std::tuple<uint8_t, uint8_t, uint16_t>>> y_ac_symbol_bitnums_groupcode;
	y_ac_symbol_bitnums_groupcode.resize(y_ac_rle_code.size());
	for (unsigned int i = 0; i < y_ac_symbol_bitnums_groupcode.size(); i++)
	{
		y_ac_symbol_bitnums_groupcode[i] = get_tuple(y_ac_rle_code[i]);
	}
	vector<vector<std::tuple<uint8_t, uint8_t, uint16_t>>> uv_ac_symbol_bitnums_groupcode;
	uv_ac_symbol_bitnums_groupcode.resize(uv_blocks_dct.size());
	for (unsigned int i = 0; i < uv_ac_symbol_bitnums_groupcode.size(); i++)
	{
		uv_ac_symbol_bitnums_groupcode[i] = get_tuple(uv_ac_rle_code[i]);
	}
	if (true)
	{
		cout << endl << endl;
		auto show_block_symbol_bitnums_groupcode = [](const vector<std::tuple<uint8_t, uint8_t, uint16_t>> block)
		{
			for (const auto& tu : block)
			{
				printf("(sumbol=%d, bitnum=%d,groupcode=%d)   ", 
					(int)std::get<0>(tu), (int)std::get<1>(tu), (int)std::get<2>(tu));
			}
			cout << endl;
		};
		cout << "---------------------------------------y 交流游程编码+分组编码再合并:----------------------" << endl;
		for (const auto& bk : y_ac_symbol_bitnums_groupcode)
		{
			show_block_symbol_bitnums_groupcode(bk);
		}
		cout << endl;
		cout << "--------------------------------------uv 交流游程编码+分组编码再合并:--------------------" << endl;
		for (const auto& bk : uv_ac_symbol_bitnums_groupcode)
		{
			show_block_symbol_bitnums_groupcode(bk);
		}
		cout << endl << endl;
	}

	//一张图片Y ac分量的频率
	int y_ac_frequency[potential_symbol_num] = {0};
	for (const auto& i: y_ac_symbol_bitnums_groupcode)
	{
		for (const auto& j : i)
		{
			y_ac_frequency[std::get<0>(j)]++;
		}
	}
	//一张图片UV ac分量的频率
	int uv_ac_frequency[potential_symbol_num] = { 0 };
	for (const auto& i : uv_ac_symbol_bitnums_groupcode)
	{
		for (const auto& j : i)
		{
			uv_ac_frequency[std::get<0>(j)]++;
		}
	}
	
	//Y ac分量的符号数
	const int y_ac_symbol_num = std::accumulate(y_ac_frequency, y_ac_frequency + potential_symbol_num, 0,
		[](int own, int cur) {return own + (cur == 0 ? 0 : 1); });
	//UV ac分量的符号数
	const int uv_ac_symbol_num = std::accumulate(uv_ac_frequency, uv_ac_frequency + potential_symbol_num, 0,
		[](int own, int cur) {return own + (cur == 0 ? 0 : 1); });

	//为一张图片生成Y分量、UV分量交流huffman表，用于由uint8_t查看其huffman码字
	//对[0,255]中的每一个数huffman编码
   //出现次数为0的pair<huffman编码长度n为0，huffman编码（16位整数的低n位）为0>
	vector<std::pair<int, uint16_t>> y_ac_each_uint8_t_huffman_code =
		get_each_uint8_t_huffman(y_ac_frequency, potential_symbol_num);
	if (true)
	{
		for (int i =0;i<potential_symbol_num;i++)
		{
			if (y_ac_frequency[i] != 0) printf("y_ac_frequency[%d] =%d  ", i, y_ac_frequency[i]);
		}
		cout << endl;

		for (const auto& i : y_ac_each_uint8_t_huffman_code)
		{
			if (i.first != 0) printf("(编码长度，码字)：(%d,%d)\n",i.first,i.second);
		}
		cout << endl;
	}
	vector<std::pair<int, uint16_t>> uv_ac_each_uint8_t_huffman_code =
		get_each_uint8_t_huffman(uv_ac_frequency, potential_symbol_num);

	auto write_huffman_table = [](const vector<std::pair<int, uint16_t>>& table, std::ofstream& ofs)
	{
		assert(table.size() == 256);

		//长度为n的码字有几个
		uint8_t nums[16] = { 0 };
		for (const auto& i : table)
		{
			if (i.first != 0)
			{
				nums[i.first-1]++;
			}
		}
		if (true)
		{
			cout << "长为i的码字有:" << endl;
			for (auto i : nums)
			{
				cout << (int)i << " ";
			}
			cout << endl;
		}
		

		vector<uint8_t> not_zerolencode_idx;
		vector<std::pair<int, uint16_t>> not_zerolencode_table;
		for (unsigned int i = 0; i < table.size(); i++)
		{
			if (table[i].first != 0)
			{
				not_zerolencode_idx.push_back((uint8_t)i);
				not_zerolencode_table.push_back(table[i]);
			}
		}

		//获取长度不为0的码字由小到大排序后的索引
		vector<uint8_t> idx(not_zerolencode_table.size());
		for (unsigned int i = 0; i < idx.size(); i++) idx[i] = (uint8_t)i;
		std::sort(idx.begin(), idx.end(),
			[&not_zerolencode_table](int lhs, int rhs) {
			return not_zerolencode_table[lhs].second <not_zerolencode_table[rhs].second;
		});

		vector<uint8_t> sorted_not_zerolencode_idx(not_zerolencode_idx.size());
		for (unsigned int i = 0; i < sorted_not_zerolencode_idx.size(); i++)
		{
			sorted_not_zerolencode_idx[i] = not_zerolencode_idx[idx[i]];
		}

		ofs.write((char*)nums, 16);
		ofs.write((char*)&sorted_not_zerolencode_idx[0], sorted_not_zerolencode_idx.size());
	};
	
	auto write = [](std::initializer_list<uint8_t> d, std::ofstream& ofs) 
	{ 
		ofs.write((char*)d.begin(),d.size());
	};
	auto write16=[](uint16_t val, std::ofstream& ofs)
	{
		//端模式：高字节在前, 低字节在后
		uint8_t high = (uint8_t)(val >> 8);
		uint8_t low = (uint8_t)(val);
		ofs << high << low;
	};
	auto write_marker = [write](std::uint8_t code,  std::ofstream& ofs) 
	{ 
		write({ 0xFF, code },ofs); 
	};

	uint8_t buf = 0, off = 0;
	auto write_bit = [write, &buf, &off](bool d, std::ofstream& ofs)
	{
		//将buf字节的off位置为1（d为true时）或0，该位之前为0
		buf |= ((uint8_t)d) << (7 - off);
		if (++off == 8)
		{
			/*
			在处理 JPG 文件时, 如果你碰到一个 0xFF, 而它后面的字节不是 0, 并且这个字节
			没有意义. 那么你遇到的 0xFF 字节必须被忽略. (一些 JPG 里, 常用用 0xFF 做某
			些填充用途) 如果你在做 huffman 编码时碰巧产生了一个 0xFF, 那么就用 0xFF 
			0x00 代替. 就是说在 jpeg 图形解码时碰到 FF00 就把它当作 FF 处理
			*/
			if (buf == 0xFF)
			{
				write({ 0xFF, 0x00 }, ofs);
			}			
			else
			{
				write({ buf }, ofs);
			}
			buf = off = 0;
		}
	};
	auto write_bits = [write_bit](std::uint16_t d, int s, std::ofstream& ofs) 
	{ 
		//将uint16_t d的低s位bit写入ofs
		for (int i = 0; i < s; ++i) 
			write_bit((d & (1 << (s - i - 1)))!=0,ofs); 
	};

	auto write_a_block = [write_bits](
		std::pair<uint8_t, uint16_t> dc_group_code, //直流分量
		vector<std::pair<int, uint16_t>>& dc_each_uint8_t_huffman_code, //直流huffman表
		const vector<std::tuple<uint8_t, uint8_t, uint16_t>>& ac_symbol_bitnums_groupcode,
		const vector<std::pair<int, uint16_t>>& ac_each_uint8_t_huffman_code,  //交流huffman表
		std::ofstream& ofs)
	{
		cout << "write_a_block():" << endl;

		//对直流分量的表示为：std::pair<长度bit数,按分组编码后的值>
		std::pair<int, uint16_t > hufCode = dc_each_uint8_t_huffman_code[dc_group_code.first];
		write_bits(hufCode.second,hufCode.first,ofs);  //RLE长度分量的huffman编码
		write_bits(dc_group_code.second,dc_group_code.first,ofs);  //RLE的码值分量
		printf("直流分量的分组编码为(长度=%d，码值=%d),  再对长度分量huffman编码(长度=%d, 码值=%d) \n",
			(int)dc_group_code.first, (int)dc_group_code.second, hufCode.first, (int)hufCode.second);

		//对交流分量的表示为：vector<std::tuple<uint8_t符号，对值分组编码后的bit数，按分组编码后的值>>
		for (unsigned int i = 0; i < ac_symbol_bitnums_groupcode.size(); i++)
		{
			uint8_t symbol = std::get<0>(ac_symbol_bitnums_groupcode[i]);
			uint8_t bitnums= std::get<1>(ac_symbol_bitnums_groupcode[i]);
			uint16_t groupcode = std::get<2>(ac_symbol_bitnums_groupcode[i]);

			std::pair<int, uint16_t > hufcode = ac_each_uint8_t_huffman_code[symbol];

			write_bits(hufcode.second,hufcode.first,ofs);  //RLE前续0数+当前值得分组编码长度构成的字节的huffman
			
			write_bits(groupcode, bitnums,ofs);   //当前值得分组编码值

			printf("交流分量的RLE前续0数+当前值得分组编码长度构成的字节的huffman编码为(长度=%d，码值=%d),  "
				"当前值分组编码为(长度=%d, 码值=%d) \n", (int)hufcode.first, (int)hufcode.second, (int)bitnums, (int)groupcode);
			
			//注意：对于游程编码后的EOB即(0,0)，第一个0为前续0数，第二个0的分组编码为(长度=0,码值0)，所以
			//组成的字节符号为0，码值为0，写入文件时，会写入字节符号0的huffman编码，码值0则不写入
		}
	};

	//开始向文件写内容
	std::ofstream ofs(dst_path, std::ios::binary);
	if (!ofs) { cout << "Couldn't open" << dst_path << endl; return; }
	write_marker(0xD8,ofs); //图片开始

	write_marker(0xE0,ofs);
	write16(16,ofs);
	write({ 0x4A, 0x46, 0x49, 0x46, 0x00, 0x01, 0x01 }, ofs);
	write({0x00,0x00, 0x01, 0x00, 0x01, 0x00, 0x00 },ofs);

	write_marker(0xDB, ofs); //Y量化表
	write16(67, ofs); 
	write({0x00},ofs);
	write(std::initializer_list<uint8_t>(&q_table_y[0], &q_table_y[q_table_y.size()-1]+1), ofs);
	write_marker(0xDB, ofs);//UV量化表
	 write16(67, ofs); 
	 write({ 0x01 }, ofs);	
	write(std::initializer_list<uint8_t>(&q_table_uv[0], &q_table_uv[q_table_uv.size() - 1]+1), ofs);
	
	write_marker(0xC0,ofs); //帧图像开始
	write({ 0x00, 0x11, 0x08 },ofs); 
	write16(height,ofs);
	write16(width, ofs);
	write({ 0x03}, ofs);
	//颜色分量ID、水平垂直采样因子、量化表ID
	write({0x00, 0x22, 0x00, 0x01, 0x11, 0x01, 0x02, 0x11, 0x01 },ofs);
	//在一个2x2像素的图像块中，Y逐点采样，U, V分别是采样一点。那么三个分量的采样因子记为：
	//	Y：SHY = 2； SVY = 2；
	//	U：SHU = 1；SVU = 1：
	//	SHV = 1； SVV = 1

	write_marker(0xC4,ofs );//哈夫曼表--直流的0号表，Y分量
	write16(2+1+16+ y_dc_symbol_num, ofs);
	write({0x00}, ofs); //直流和交流分开编ID
	write_huffman_table(y_dc_each_uint8_t_huffman_code, ofs);
	
	write_marker(0xC4, ofs);//哈夫曼表--交流的0号表，Y分量
	write16(2 + 1 + 16 + y_ac_symbol_num, ofs);
	write({ 0x10 }, ofs);
	write_huffman_table(y_ac_each_uint8_t_huffman_code, ofs);

	write_marker(0xC4, ofs);//哈夫曼表--直流的1号表，UV分量
	write16(2 + 1 + 16 + uv_dc_symbol_num, ofs);
	write({ 0x01 }, ofs);
	write_huffman_table(uv_dc_each_uint8_t_huffman_code, ofs);
	
	write_marker(0xC4, ofs);//哈夫曼表--交流的为1号表，UV分量
	write16(2 + 1 + 16 + uv_ac_symbol_num, ofs);
	write({ 0x11 }, ofs);
	write_huffman_table(uv_ac_each_uint8_t_huffman_code, ofs);


	write_marker(0xDA,ofs);  //行扫描开始
	write16(12, ofs);
	write({3},ofs);
	write({0,0x00}, ofs);//Y, 0x直流分量使用的哈夫曼树编号+交流分量使用的哈夫曼树编号
	write({1,0x11}, ofs);  //U,颜色分量ID可以是0、1、2
	write({2,0x11}, ofs); //V
	write({ 0x00, 0x3F, 0x00 }, ofs);

	//之后为真正的图像信息
	//由于Y、U、V水平采样因子分别为2、1、1，垂直采样因子分别为2、1、1，因此一个最
	//小编码单元MCU对应16*16像素块，即4个8*8像素块。一个MCU包括四个YDU、1个UDU
	//和一个VDU，数据流组织形式为YDU,YDU,YDU,YDU,UDU,VDU
	const int ydu_num = (height / 8)*(width / 8);
	const int udu_num= (height /2/ 8)*(width/2 / 8);
	const int vdu_num = (height / 2 / 8)*(width / 2 / 8);
	const int mcu_num = (height / 2 / 8)*(width / 2 / 8);
	assert(y_dc_group_code.size()== y_ac_symbol_bitnums_groupcode.size());
	assert(y_dc_group_code.size() == ydu_num);
	assert(uv_dc_group_code.size() == uv_ac_symbol_bitnums_groupcode.size());
	assert(uv_dc_group_code.size() == (udu_num+vdu_num));
	for (int i = 0; i < mcu_num; i++)
	{
		//第i号MCU包含的YDU的索引为idx[4]，包含的UDU的索引为i
		//写YDU,YDU,YDU,YDU
		const int mcu_num_per_row = width / 16;
		const int r = i/mcu_num_per_row;
		const int c = i%mcu_num_per_row;
		int idx[4];
		idx[0] = r*mcu_num_per_row * 4 + c * 2;
		idx[1] = idx[0] + 1;
		idx[2] = idx[0] + width / 8;
		idx[3] = idx[2] + 1;

		cout << endl << "********************************************" << endl;
		cout << i << "号MCU的四个Block的索引是:" ;
		for (auto i : idx)cout << i << "  ";
		cout << endl << "*********************************************" << endl << endl;

		for (int j = 0; j < 4; j++) 
		{
			write_a_block(y_dc_group_code[idx[j]], y_dc_each_uint8_t_huffman_code,
				y_ac_symbol_bitnums_groupcode[idx[j]], y_ac_each_uint8_t_huffman_code,ofs);
		}


		cout << endl << "-------写UDU-------------" << endl << endl;
		//写UDU
		write_a_block(uv_dc_group_code[i], uv_dc_each_uint8_t_huffman_code,
			uv_ac_symbol_bitnums_groupcode[i], uv_ac_each_uint8_t_huffman_code, ofs);
		cout << endl << "-------写VDU-------------" << endl << endl;
		//写VDU  
		int v_idx = i + uv_dc_group_code.size() / 2;
		write_a_block(uv_dc_group_code[v_idx], uv_dc_each_uint8_t_huffman_code,
			uv_ac_symbol_bitnums_groupcode[v_idx], uv_ac_each_uint8_t_huffman_code, ofs);
	}

	//填充1
	if (off != 0)
	{
		cout << "末尾填充" << 8 - off << "个1" << endl;
		write_bits(0xFF,8-off,ofs);
	}

	write_marker(0xD9,ofs);  //EOI(End of Image) 图像结束标记
}

int main()
{
	//test_syntex();

	//test_get_depth_in_huffman_tree();
	//test_canonical_huffman_encode();
	//test_rgb24_to_yuv411P();
	//test_padto();
	//test_dct();
	//test_zigzag_sort();
	//test_quantify();
	//test_rle_encode();
	//test_group_encode();
	//test_get_y_blocks_get_uv_blocks();
	//test_get_each_uint8_t_frequency_write_huffman_table();


	//const std::string path = "C:/Users/myth/Pictures/a.jpg";
	//Mat  src = imread(path);
	//if (src.empty()) { cout << "Couldn't read " << path << endl; return -1; }
	//cv::imshow("src", src); waitKey(0);
	//vector<uint8_t> rgbdata(src.rows*src.cols * 3);
	//memcpy_s(&rgbdata[0], src.rows*src.cols * 3,src.data, src.rows*src.cols * 3);

	constexpr int width = 16;
	constexpr  int height = 16;
	uint8_t data[width *height * 3];
#if 0
	uint8_t* p = data;
	for (int i = 0; i <height; ++i)
		for (int j = 0; j < width; ++j)
			*p++ = 255, *p++ = 0, *p++ = 255;     //G B R
#else
	for (int i = 0; i < width *height ; i++)
	{
		data[3*i] = i;
		data[3*i + 1] = i;
		data[3*i + 2] = i;
	}
#endif
	vector<uint8_t> rgbdata(height *width * 3);
	memcpy_s(&rgbdata[0], height *width * 3, data, height *width * 3);
	const char *output = "./test.jpg";
	avjpeg(rgbdata, width, height, output);
	cout << "完成编码，输出保存在" << output << endl;
}