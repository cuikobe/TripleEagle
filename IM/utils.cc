#include "utils.h"

namespace _Cide{
	std::vector< std::vector<int> > graphT;
    //std::vector<bool> is_selected;//标记该结点是否已经被选

    float lambda;
    float B_for_Twin;
    float Budget;
    int max_node;
    int node_num;
    int edge_num;

    std::vector<float> node_cost;
    std::default_random_engine e_RandomTM;
    std::default_random_engine e_BFM_RAN;
}


using namespace std;

float getCurrentMemoryUsage() { //megabyte verii
	string pid = intToStr(unsigned(getpid()));
	string outfile = "temp/tmp_" + pid + ".txt";
	string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;	
	system(command.c_str());
	
	string mem_str;
	ifstream ifs(outfile.c_str());
	std::getline(ifs, mem_str);
	ifs.close();
	
	mem_str = mem_str.substr(0, mem_str.size()-1);
	float mem = (float)strToInt(mem_str);
	
	return mem/1024; // in MB
}

float getRunningTime(time_t startTime) {	//in seconds
    time_t curTime;
    time(&curTime);
    float duration = ((float)(curTime - startTime)); // seconds icin
    return duration;
}

void stringTokenizer(string& str, float *tokens, int size, string& delimiters) {
	//Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);//跳过最初的分隔符
	//Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);//找到下一个分隔符
	
	int i = 0;
	while (string::npos != pos || string::npos != lastPos)
	{
		if(i >= size) {
			cout << "something is wrong" << endl;
			exit(1);
		}
		
		tokens[i++] = strToFloat(str.substr(lastPos, pos - lastPos));//通过分隔符间断来读取这一行的数据
		lastPos = str.find_first_not_of(delimiters, pos);		
		pos = str.find_first_of(delimiters, lastPos);
	}
}

void stringTokenizer(string& str, double *tokens, int size, string& delimiters) {
	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);
	
	int i = 0;
	while (string::npos != pos || string::npos != lastPos)
	{
		if(i >= size) {
			cout << "something is wrong" << endl;
			exit(1);
		}
			
		tokens[i++] = strToDouble(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);		
		pos = str.find_first_of(delimiters, lastPos);
	}
}


void rtrim(char *str) {
	size_t n;
	n = strlen(str);
	while (n > 0 && isspace((unsigned char)str[n - 1])) {
		n--;
	}
	str[n] = '\0';
}

void ltrim(char *str) {
	size_t n;
	n = 0;
	while (str[n] != '\0' && isspace((unsigned char)str[n])) {
		n++;
	}
	memmove(str, str + n, strlen(str) - n + 1);
}

void trim(char *str) {
	rtrim(str);
	ltrim(str);
}
