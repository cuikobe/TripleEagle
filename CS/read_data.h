//
// Created by CLAKERS on 2020/8/22.
//
#ifndef IMAGE_READ_DATA_H
#define IMAGE_READ_DATA_H
#include <iostream>
#include "fstream"
#include "time.h"
#include "vector"
#include "random"
#include "set"
#include "algorithm"
#include "list"
#include "map"
using namespace std;
const int node_num=3000;
const int m=3072;
const int category=10;
const int per_category_max=10;
const int k=1;
vector<int> lable;
vector<vector<int>> feature;
vector<double> contrast_cost(node_num,0.0);
vector<vector<double>> similarity;
int K;
const int d=1;
double min_cost=999999999;
double max_cost=-999999999;
string cost_text="./rdata/category/rd_contrast.txt";
string feature_text="./rdata/category/rd_feature.txt";
const double ave_num=10;
const int max_image=10;
void read()
{
    feature.resize(node_num, vector<int>(m));
    ifstream in1(feature_text);
    //ifstream in1("feature.txt");
    int i=0;
    int templable=-1;
    while(!in1.eof())
    {
        in1>>templable;
        if (in1.fail())
            break;
        lable.push_back(templable);
        for(int j=0;j<m;j++)
            in1 >> feature[i][j];

        i++;
    }
    in1.close();

    cout<<"lable num: "<<lable.size()<<" | "<<i<<endl;

    ifstream in2(cost_text);
    //ifstream in2("contrast.txt");
    i=0;
    double temp_cost;
    double sum_cost=0.0;
    while(!in2.eof())
    {
        in2>>temp_cost;
        if (in2.fail())
            break;

        if(temp_cost<min_cost)
            min_cost=temp_cost;
        if(temp_cost>max_cost)
            max_cost=temp_cost;

        contrast_cost[i]=temp_cost;
        sum_cost+=temp_cost;
        i++;
        //cout<<temp_cost<<endl;
    }
    in2.close();

    //cout<<"normalize: "<<endl;
    //double ave_cost=sum_cost/node_num;
    //*
    for(auto &x:contrast_cost)
    {
        x=(node_num*x)/(sum_cost*ave_num);
        //cout<<x<<endl;
    }
    //*/

}
void cal_similarity()
{
    similarity.resize(node_num, vector<double>(node_num));
    //inner product
    /*
    for(int i = 0;i < feature.size();i++)
    {
        for(int j = i;j < feature.size();j++)
        {
            double product = 0.0;
            for(int k = 0;k < feature[0].size();k++)
            {
                product+=(feature[i][k]*feature[j][k]);
            }
            similarity[i][j] = product;
            similarity[j][i] = product;
        }
    }
     */
    //distance
    /*
    for(int i = 0;i < feature.size();i++)
    {
        for(int j = i;j < feature.size();j++)
        {
            double distance = 0.0;
            for(int k = 0;k < feature[0].size();k++)
            {
                distance += pow((feature[i][k]-feature[j][k]),2);
            }
            distance = sqrt(distance);
            similarity[i][j] = distance;
            similarity[j][i] = distance;
        }
    }
     //*/
    //consine
    //*
    for(int i = 0;i < feature.size();i++)
    {
        for(int j = i;j < feature.size();j++)
        {
            double product = 0.0;
            double module_a=0.0;
            double module_b=0.0;
            for(int k = 0;k < feature[0].size();k++)
            {
                product+=(feature[i][k]*feature[j][k]);
                module_a+=pow(feature[i][k],2);
                module_b+=pow(feature[j][k],2);
            }
            module_a=sqrt(module_a);
            module_b=sqrt(module_b);
            similarity[i][j] = product/(module_a*module_b);
            similarity[j][i] = product/(module_a*module_b);
        }
    }
    // */
}

//fixed
double f_u(int node)
{
    double sum_value=0.0;
    for(int i=0;i<node_num;i++)
    {
        sum_value += similarity[i][node];
    }
    return sum_value-similarity[node][node]/(node_num);
}
bool node_budget_feasible(const int &node,double Budget)
{
    /**********check single node Budget constraint**********/
    if(contrast_cost[node]>Budget)
        return false;
    else
        return true;
}
class S_class {
public:
    S_class() {
        max_marginal=-999999999.0;
        max_element=-1;
        S_cost=0.0;
        S_revenue=0.0;
        for(int iter=0;iter<category;iter++)
            category_sum.push_back(0);

        selected.resize(node_num,0);
        node_price.resize(node_num,0.0);
    }
    double S_cost;
    double S_revenue;
    vector<int> Set;

    double max_marginal;
    int max_element;
    vector<int> selected;
    vector<int> category_sum;
    /***********only for SIM***************/
    vector<double> node_price;
    bool is_feasible(const int &e,int max_image)
    {
        if(Set.size()>=max_image||category_sum[lable[e]]>=per_category_max)
            return false;
        else
            return true;
    }
    void clear()
    {
        S_revenue = 0.0;
        S_cost=0.0;

        Set.clear();
        fill(selected.begin(), selected.end(), 0);
        fill(node_price.begin(), node_price.end(), 0.0);
        for(int iter=0;iter<category;iter++)
            category_sum[iter]=0;

        max_marginal=-999999999.0;
        max_element=-1;
    }
    void replace_with_singleton(const double &marginal,const int &node)
    {
        clear();
        add_element(marginal,node);
    }
    void add_element(const double &marginal,const int &node)
    {
        selected[node] = 1;
        Set.push_back(node);
        S_revenue += marginal;
        S_cost+=contrast_cost[node];

        category_sum[lable[node]]++;//S+e\in I
    }
    void add_element_truth(const double &marginal,const int &node,const double &price)
    {
        selected[node] = 1;
        Set.push_back(node);
        S_revenue += marginal;
        S_cost+=price;

        category_sum[lable[node]]++;//S+e\in I
        node_price[node]=price;
    }
    bool budget_feasible(const int &node,double Budget)
    {
        if(S_cost+contrast_cost[node]>Budget)
            return false;
        else
            return true;
    }
    double f_S()
    {
        if(Set.empty())
            return 0.0;

        double v1=0.0;
        double v2=0.0;
        for(int i=0;i<node_num;i++)
        {
            double max_temp=-999999999;
            for(int j=0;j<Set.size();j++)
            {
                if(similarity[i][Set[j]]>max_temp)
                    max_temp=similarity[i][Set[j]];
            }
            v1+=max_temp;
        }
        for(int i=0;i<Set.size();i++)
        {
            for(int j=0;j<Set.size();j++)
            {
                v2+=similarity[Set[i]][Set[j]];
            }
        }
        return v1-v2/(node_num);
    }
    double marginal(int e)
    {
        if(Set.empty())
        {
            return f_u(e);
        }
        //slow
        double v1=0.0;
        double v2=0.0;

        for(int i=0;i<Set.size();i++)
        {
            if(e==Set[i]) return 0.0;

            v2+=similarity[e][Set[i]];
            v2+=similarity[Set[i]][e];
        }
        v2+=similarity[e][e];

        for(int i=0;i<node_num;i++)
        {
            double max_similarity=-999999999;
            for(auto j:Set)
            {
                if(similarity[i][j]>max_similarity)
                    max_similarity=similarity[i][j];
            }

            if(similarity[i][e]>max_similarity)
            {
                v1+=(similarity[i][e]-max_similarity);
            }
        }

        return v1-v2/(node_num);
    }
    double S_sub_u(int e)
    {
        double v1=0.0;
        double v2=0.0;
        for(int i=0;i<node_num;i++)
        {
            double max_temp=-999999999;
            for(int j=0;j<Set.size();j++)
            {
                if(Set[j]==e)
                    continue;
                if(similarity[i][Set[j]]>max_temp)
                    max_temp=similarity[i][Set[j]];
            }
            v1+=max_temp;
        }
        for(int i=0;i<Set.size();i++)
        {
            if(Set[i]==e)
                continue;
            for(int j=0;j<Set.size();j++)
            {
                if(Set[j]==e)
                    continue;
                v2+=similarity[Set[i]][Set[j]];
            }
        }
        return v1-v2/(node_num);
    }
    void copy(const S_class &temp)
    {
        S_revenue=temp.S_revenue;
        S_cost=temp.S_cost;
        Set.assign(temp.Set.begin(),temp.Set.end());
        selected.assign(temp.selected.begin(),temp.selected.end());
        //category_sum.assign(temp.category_sum.begin(),temp.category_sum.end());
    }
};
class Result
{
public:
    Result(){}
    Result(double rev,double cos,int siz,long long int ora)
    {
        revenue=rev;
        cost=cos;
        size=siz;
        oracle=ora;
    }
    Result(double rev,long long int ora,double uti_der,double ora_der,int round)
    {
        revenue=rev;
        oracle=ora;
        utility_deviation=uti_der;
        oracle_deviation=ora_der;
    }
    double revenue;
    long long int oracle;
    long long int round;
    double cost;
    int size;
    long long int time=0;
    long long int max_query=0;
    double utility_deviation;
    double oracle_deviation;
};

#endif //IMAGE_READ_DATA_H
