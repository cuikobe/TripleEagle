//
// Created by CLAKERS on 2022/10/4.
//

#ifndef IMAGE_GENSM_slow_H
#define IMAGE_GENSM_slow_H
#include "read_data.h"
#include "GENSM_slow_matroid.h"

class Alg1_queue
{
public:
    Alg1_queue(int node_temp,double t_weight)
    {
        node=node_temp;
        weight=t_weight;
        l=1;
    }
    int node;
    double weight;
    int l;
};
bool cmp_alg1(const Alg1_queue &a,const Alg1_queue &b)
{
    return a.weight>b.weight;
}
double ALG1(double eps,double B,const vector<int> &groundset,long long int &oracle_times)
{
    /*
    S_gamma t1,t2,t3,t4;
    vector<int> v1={227,174};
    vector<int> v2={};

    t1.Set=v1;
    t2.Set=v2;

    cout<<t1.f_S()<<endl;
    double c=0.0;
    for(auto e:t1.Set)
        c+=contrast_cost[e];
    cout<<c<<endl;

    cout<<t2.f_S()<<endl;

    c=0.0;
    for(auto e:t2.Set)
        c+=contrast_cost[e];
    cout<<c<<endl;
    //*/


    K=floor(B/min_cost);
    double p=0.5;

    //eps/=14.0;
    eps=0.0;

    //long long int oracle_times=0;


    S_class S;

    list<Alg1_queue> A;
    list<Alg1_queue> A_augment;
    for(int iter=0;iter<node_num;iter++)
    {
        if(contrast_cost[iter]<=B&&groundset[iter]==1)
        {
            //oracle_times++;

            double value=f_u(iter);
            Alg1_queue temp(iter,value/contrast_cost[iter]);
            A.push_back(temp);

            Alg1_queue temp2(iter,value);
            A_augment.push_back(temp2);
        }
    }
    //sort(A.begin(),A.end(),cmp_smkrandom);
    A.sort(cmp_alg1);

    bernoulli_distribution u(p);
    while(!A.empty())
    {
        //argmax i
        int node_i=-1;
        double max_marginal_cost=-999999999;

        //only used for calculate oracle times
        oracle_times+=A.size();

        for(auto &u:A)
        {
            double old_wt=u.weight;

            if(max_marginal_cost>=old_wt) break;
            //update weight
            //oracle_times++;
            u.weight=(S.marginal(u.node)/contrast_cost[u.node]);

            //cout<<"now weight: "<<u.weight<<endl;
            //cout<<"old weight: "<<old_wt<<endl;

            if(u.weight>=(old_wt/(1+eps)))
            {
                node_i=u.node;
                max_marginal_cost=u.weight;
                break;
            }
            else
            {
                if(u.weight>max_marginal_cost)
                {
                    max_marginal_cost=u.weight;
                    node_i=u.node;
                }
                u.l++;
            }
            /*
            if((*u).l>(log(K/eps)/log(1+eps)))
            {
                u=A.erase(u);
                //cout<<"happen!"<<endl;
            }
            //*/
            //sort(A.begin(),A.end(),cmp_smkrandom);
        }
        if(node_i==-1)
        {
            cout<<"la ji"<<endl;
            break;
        }
        if(max_marginal_cost<0)
            break;

        double random=u(e_GENSM);
        if(random==1)
        {
            S.Set.push_back(node_i);
            S.S_cost+=contrast_cost[node_i];
            S.S_revenue+=max_marginal_cost*contrast_cost[node_i];
        }
        for(auto e=A.begin();e!=A.end();)
        {
            if( (S.S_cost+contrast_cost[(*e).node]>B) || ((*e).node==node_i) || ((*e).weight<=0)
               )
            {
                e=A.erase(e);
            }
            else
            {
                e++;
            }
        }
        A.sort(cmp_alg1);
    }
    //cout<<"now oracle: "<<oracle_times<<endl;

    //getmax for augment
    A_augment.sort(cmp_alg1);
    S_class S_star=S;
    //j=0,ie. argmax f(u)
    int i_star=-1;
    double i_star_value=-999999999;

    /*
    for(int iter=0;iter<node_num;iter++)
    {
        if(contrast_cost[iter]<=B)
        {
            //oracle_times++;
            double value=f_u(iter);
            if (value > i_star_value)
            {
                i_star_value = value;
                i_star = iter;
            }
        }
    }*/
    auto e_temp=A_augment.begin();
    i_star=(*e_temp).node;
    i_star_value=(*e_temp).weight;
    if(i_star_value>S_star.S_revenue)
    {
        S_star.Set.clear();
        S_star.Set.push_back(i_star);
        S_star.S_revenue=i_star_value;
        S_star.S_cost=contrast_cost[i_star];
    }

    //j>=1
    vector<int> Sj_selected(node_num,0);
    S_class Sj;
    for(int j=0;j<S.Set.size();j++)
    {
        //build Sj
        int temp_e=S.Set[j];//Sj-1 and temp_e = now Sj

        Sj_selected[temp_e]=1;
        Sj.S_revenue+=Sj.marginal(temp_e);//just for running quickly, so it is not included in oracle times;
        Sj.Set.push_back(temp_e);
        Sj.S_cost+=contrast_cost[temp_e];

        int uj=-1;
        double uj_value=-999999999;

        for(auto e=A_augment.begin();e!=A_augment.end();)
        {
            if( (Sj.S_cost+contrast_cost[(*e).node]>B) || (Sj_selected[(*e).node]==1) || ((*e).weight<=0)
                 )
            {
                e=A_augment.erase(e);
            }
            else
            {
                e++;
            }
        }

        //sort(A_augment.begin(),A_augment.end(),cmp_smkrandom);
        A_augment.sort(cmp_alg1);

        //only used for calculate oracle times
        oracle_times+=A_augment.size();

        for(auto &u:A_augment)
        {
            double old_wt=u.weight;
            if(uj_value>=old_wt) break;

            //update weight
            //oracle_times++;
            u.weight=Sj.marginal(u.node);
            if(u.weight>=(old_wt/(1+eps)))
            {
                uj=u.node;
                uj_value=u.weight;
                break;
            }
            else
            {
                if(u.weight>uj_value)
                {
                    uj_value=u.weight;
                    uj=u.node;
                }
                u.l++;
            }
            /*
            if((*u).l>((log((double)node_num/eps)/log(2.0))/eps))
            {
                u=A_augment.erase(u);
            }
            //*/
            //sort(A_augment.begin(),A_augment.end(),cmp_smkrandom);
        }
        if(uj==-1) continue;

        double f_Sj_and_uj=Sj.S_revenue+uj_value;
        if(f_Sj_and_uj>S_star.S_revenue)
        {
            S_star.copy(Sj);
            S_star.S_revenue=f_Sj_and_uj;
            S_star.S_cost+=contrast_cost[uj];
            S_star.Set.push_back(uj);
        }
    }

    return S_star.S_revenue;
}

Result GENSM(double eps,double B) {
    double q1=0.201;
    double q2=1.0/2.0;
    bernoulli_distribution u1(q1);
    bernoulli_distribution u2(q2);

    cout << "GENSM & Budget: " << B << "---------start---------" << endl;
    double beta = 9.185;
    long long int oracle_times = 0;
    if (u1(e_GENSM) == 1) {
        cout << "return single element ! " << endl << endl;
        double f_u_max = -1.0;
        int u_star;

        for (int iter = 0; iter < node_num; iter++) {
//            if (!node_budget_feasible(iter, B))
//                continue;

            oracle_times++;
            double value = f_u(iter);

            if (value >= f_u_max) {
                f_u_max = value;
                u_star = iter;
            }
        }
        S_class S;
        S.add_element(f_u_max, u_star);
        return Result(S.S_revenue, 0.201, S.Set.size(), oracle_times);
    }

    vector<int> available_self(node_num,0);//mark all nodes are available ie.,N_i
    vector<int> available_alg1(node_num,0);//mark all nodes are available ie.,N_i
    for (int iter = 0; iter < node_num; iter++) {
        if (u2(e_GENSM) == 1) {
            available_self[iter] = 1;
        } else {
            available_alg1[iter] = 1;
        }
    }
    double x=ALG1(eps,B,available_alg1,oracle_times);
    cout<<"x: "<<x<<endl;


    int ell=2;
    vector<S_class> S;
    vector<multimap<double,Ai_node>> A;
    vector<double> B_array;
    for(int i=0;i<ell;i++) {
        S.emplace_back();
        A.emplace_back();
        B_array.push_back(B);
    }

    //initial Ai

    for (int iter = 0; iter < node_num; iter++) {
        if (available_self[iter] == 0)
            continue;

//        if (!node_budget_feasible(iter, B)) {
//            available_self[iter] = 0;
//            continue;
//        }

        Ai_node temp(iter);
        double value = f_u(iter);
        for (int j = 0; j < ell; j++) {
            A[j].insert(pair<double,Ai_node>(value,temp));
        }
    }

    while (true) {
        for (int j = 0; j < S.size(); j++)
        {
            S[j].max_marginal = -999999999.0;
            S[j].max_element = -1;
            //A=empty, then return element=-1 and marginal gain =0
            pair<int, double> temp = Boost(A[j], S[j], available_self, B);
            S[j].max_marginal = temp.second;
            S[j].max_element = temp.first;

            //only used for calculate oracle times
            for (int iter = 0; iter < node_num; iter++) {
                if (available_self[iter] == 0)
                    continue;
                oracle_times++;
            }
        }
        double max_marginal = 0.0;
        int max_element=-1;
        int max_solution = -1;
        for (int j = 0; j < S.size(); j++) {
            if (S[j].max_marginal > max_marginal) {
                max_solution = j;
                max_marginal = S[j].max_marginal;
                max_element = S[j].max_element;
            }
        }
        if (max_solution == -1 || max_marginal <= 0.0)//non element or marginal gain<=0
            break;
        //S<-S\cup {u}
        double price=max_marginal*beta*B/x;
        if(price>=contrast_cost[max_element]&&price<=B_array[max_solution])
        {
            S[max_solution].add_element_truth(max_marginal,max_element,price);
            B_array[max_solution]-=price;
        }
        available_self[max_element] = 0;//discarded or selected
    }

    S_class S_star;
    for(const auto& it:S)
    {
        if(it.S_revenue>=S_star.S_revenue)
            S_star=it;

        S_class temp=USM(it,oracle_times);
        if(temp.S_revenue>S_star.S_revenue)
            S_star=temp;

    }

    cout<<"S*:"<<endl;
    cout<<"  revenue: "<<S_star.S_revenue<<" size: "<<S_star.Set.size()<<" cost: "<<S_star.S_cost<<endl;
    cout<<"  all nodes: "<<endl;

    for(const auto &p:S_star.Set)
        cout<<p<<'\t';
    cout<<endl;

    cout<<"oracle times: "<<oracle_times<<endl;

    cout<<"real revenue: "<<S_star.f_S()<<endl;

    cout<<"GENSM ---------end--------- "<<endl<<endl;
    return Result(S_star.S_revenue,S_star.S_cost,S_star.Set.size(),oracle_times);
}

#endif //IMAGE_GENSM_H
