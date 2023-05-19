#ifndef ALLOCATOR_H
#define ALLOCATOR_H


#include <ctime>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include "anyoption.h"
#include "utils.h"
#include "advertiser.h"
#include "TimGraph.h"
#include <cmath>
#include <algorithm>
#include <list>
#include <unordered_set>
#include "random"
#include "time.h"

namespace _Cide{
    void copy_tim(TimGraph *src,TimGraph *dest);
    typedef std::vector<advertiser*> advertiserList;//把advertiserList定义为vector<advertiser*>数组的别名
    typedef std::vector<TimGraph*> TimGraphList;

    class allocator {
        
    public:
        time_t startTime;

        float alpha;
        AnyOption *opt;//option对象
        int n, m, nrTopics, nrCompanies;

        string delim, probGraphFile;
        float epsilon;
        float theta_0;
        int b = -1;

        advertiserList *advList;
        //TimGraphList *timList1,*timList2,*timListTemp1,*timListTemp2;
        //TimGraph *timTemp,*timTempForCalVal;

        TimGraph *timCS;

        void TwinGreedyFastFinal();
        void ResidualRandom();
        void TwinGreedy();
        void SampleGreedyFinal();
        void greedy(list<pair<int,int>> M);
        void fantom();
        float pairgreedy(TimGraphList *timlist);

        void RandomTM_acc();
        void RandomTM();
        void IP();
        void IP_acc();
        void BFM_RAN();
        void BFM_DET();
        void BFM_NM();

        bool first_u1;
        bool first_u2;

        void copy_tim_fantom(TimGraph *src, TimGraph *dest);

        TimGraphList *timList;

        allocator(AnyOption* opt);
        ~allocator();

        //int windowSize;

        void readTICGraph();
        void readItemDistsFile();
        void readIncentiveCosts();

        float USM(vector<pair<int, int>> S);

        double cal_f_value_for_usm(int H_num, int u);
    };
    class Node
    {
        public:
        //Node();
        Node(int ind, float pri)
        {
            index=ind;
            price=pri;
        };
        Node(int ind)
        {
            index=ind;
        };
        int index;
        float price;
    };
    class S_class {
    public:
        S_class()
        {
           // cout<<"call the null construct function, be carefull !"<<endl;
        };
        S_class(advertiser *adv,TimGraph *timCS) {
            max_marginal=-999999999.0;
            max_element=-1;
            S_cost=0.0;
            S_price=0.0;
            S_revenue=0.0;
            selected.resize(node_num,0);

            graph = new _Cide::TimGraph(adv, 0.1, node_num, edge_num);
            //generate 每个广告商的RR集
            copy_graph(timCS,graph);
        }
        TimGraph *graph;
        float S_cost;
        float S_price;
        float S_revenue;
        vector<Node> Set;

        vector<int> selected;
        float max_marginal;
        int max_element;
        void copy_graph(TimGraph *src, TimGraph *dest)
        {
            dest->theta=src->theta;
            //dest->theta_old=src->theta;
            //清空原vector数组并复制
            dest->hyperGT_adv.assign(src->hyperGT_adv.begin(),src->hyperGT_adv.end());
            dest->isCovered.assign(src->isCovered.begin(),src->isCovered.end());
            dest->hyperG_adv.assign(src->hyperG_adv.begin(),src->hyperG_adv.end());
            dest->hyper_degree.assign(src->hyper_degree.begin(),src->hyper_degree.end());
        }
        /***********new method***************/
        void clear(advertiser *adv,TimGraph *timCS)
        {
            S_revenue = 0.0;
            S_cost=0.0;
            S_price=0.0;
            Set.clear();
            fill(selected.begin(), selected.end(), 0);
            copy_graph(timCS,graph);

            max_marginal=-999999999.0;
            max_element=-1;
        }
        void replace_with_singleton(const float &marginal,const Node &node,advertiser *adv,TimGraph *timCS)
        {
            clear(adv,timCS);
            add_element(marginal,node);
        }
        void add_element(const float &marginal,const Node &node)
        {
            graph->opim_assign_best_node(node.index);
            selected[node.index] = 1;
            Set.push_back(node);
            S_revenue += marginal;
            S_cost+=node_cost[node.index];
        }
        void add_element_truth(const float &marginal,const Node &node,const double &price)
        {
            graph->opim_assign_best_node(node.index);
            selected[node.index] = 1;
            Set.push_back(node);
            S_revenue += marginal;
            S_cost+=node_cost[node.index];
            S_price+=price;
        }
//        bool budget_feasible(const Node &node,double Budget)
//        {
//            if(S_cost+node_cost[node.index]>Budget)
//                return false;
//            else
//                return true;
//        }
        //need fix for our target fucntion, which is not sames as the nips2020
        float f_S(advertiser *adv,TimGraph *timCS) {
            TimGraph *graphTemp=new _Cide::TimGraph(adv, 0.1, node_num, edge_num);;
            copy_graph(timCS, graphTemp);

            int sum_node = 0;
            //float sum_cost=0.0;
            for (auto e:Set) {
                sum_node += graphTemp->hyper_degree[e.index];
                graphTemp->opim_help_cal_f(e.index);

                //sum_cost+=node_cost[e.index];//there need fix
            }
            return (float) node_num * (float)sum_node / graphTemp->theta;
        }
        //need fix for our target fucntion, which is not sames as the nips2020
        float marginal(const Node &e)
        {
            return ((float) node_num * ((float) graph->hyper_degree[e.index] / graph->theta));
        }
        void copy(const S_class &temp)
        {
            S_revenue=temp.S_revenue;
            S_cost=temp.S_cost;
            S_price=temp.S_price;
            Set.assign(temp.Set.begin(),temp.Set.end());
            selected.assign(temp.selected.begin(),temp.selected.end());
            copy_graph(temp.graph,graph);
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
        double revenue;
        long long int oracle;
        double cost;
        int size;
    };
    extern vector<Result> twin_result;
    extern vector<Result> randomTM_result;
    extern vector<Result> IP_result;
    extern vector<Result> BFM_RAN_result;
    extern vector<Result> BFM_DET_result;
    extern vector<Result> BFM_NM_result;
    extern vector<Result> randomTM_average;

    float f_u(const int &node,const TimGraph *timTemp);
    bool node_budget_feasible(const int &node,float Budget);
    pair<int,float> Boost(multimap<double,int> &A,S_class &Si,const vector<int> &available,const double &B);
    pair<int,float> Boost_density(multimap<double,int> &A,S_class &Si,const vector<int> &available,const double &B);
    //void INITIALIZE(vector<int> &C,int &v_star,float Budget);
    void PRICING(S_class &A,const vector<int>&C,float neta,vector<float>&pi,const int &except_element,long long int &oracle_times,vector<long long int> &ask_times);
}


#endif
