//
// Created by CLAKERS on 2022/10/8.
//

#ifndef IMAGE_GENSM_SLOW_MATROID_H
#define IMAGE_GENSM_SLOW_MATROID_H
#include "RandomClockMatroid.h"
#include "Simultaneous.h"
default_random_engine e_GENSM(12345);
double q1=1.0/3.0;
double q2=1.0/2.0;
bernoulli_distribution u1(q1);
bernoulli_distribution u2(q2);
pair<int,double> GetMax(long long int &oracle_times,double eps,multimap<double,Ai_node> &A,S_class &Si,const vector<int> &available,const double &B)
{
/********check this when pop a node*************/
    eps=0.0;
    double m=0.0;//always record weight of aij
   int ai=-1;
    oracle_times+=A.size();
    while(!A.empty())
    {
        //check the end node, i.e., the node with biggest weight
        auto it=A.end();
        it--;

        double old_value=it->first;
        //if the value of maximum element is less than 0, then all element has value less than 0 due to the diminish property of submodularity, so we break and return empty element.
        if(old_value<=0)
        {
            A.clear();
            ai=-1;
            m=0.0;
            break;
        }
        //if not satisfy the k-system constraint or has been selected by other solution, then delete it and pop the next element
        //if(!Si.all_feasible(it->second.node,B)||available[it->second.node.product][it->second.node.user]==0)
        if(!Si.is_feasible(it->second.node,max_image)||available[it->second.node]==0)
        {
            A.erase(it);
            continue;
        }
        //if the element is useful, then we compute its new weight, and let its update numbers +1
        //oracle_times++;
        double new_value=Si.marginal(it->second.node);
        //the value of element in multimap can be change directly, but the key can not be change, we need to erase it and re-insert it
        it->second.l++;
        //if the value of the element diminishes not much, we can return it
        if(new_value>=old_value/(1.0+eps))
        {
            ai=it->second.node;
            m=new_value;
            break;
        }
            //or we re-insert now element and then check the next element
        else
        {
            Ai_node temp=it->second;
            //erase the element to update its weight
            A.erase(it);
            //if its update numbers is not greater than max, then we re-insert it into the queue
            if(temp.l<=(log(node_num*2.0/eps)/log(1.0+eps))) {
                A.insert(pair<double, Ai_node>(new_value, temp));
            }
        }
    }
    return pair<int,double>(ai,m);
}
double ALG4(double eps,double B,const vector<int> &groundset, long long int &oracle_times)
{
    //cout<<"ALG4 & Budget: "<<B<<"---------start---------"<<endl;

    int n=node_num;
    //double eps_pie=eps*(2.0*d+k+sqrt(2.0*d+k))/d;
    double eps_pie=eps/(2.0*d+k+sqrt(2.0*d+k))*d;
    //cout<<eps_pie<<endl;

    //random
    double p=2.0/(1.0+sqrt((double)k+2.0*d));
    //cout<<"original p: "<<p<<endl;
    //p=1.0;
    //cout<<"now p: "<<p<<endl;

    int ell=2;
    //determined
//    double p=1.0;
//    int ell=ceil(sqrt(k))+1;

    bernoulli_distribution u(p);
    //long long int oracle_times=0;

    double f_u_max=-1.0;
    int u_star=-1;

        for (int iter = 0; iter < node_num; iter++) {

            if(!node_budget_feasible(iter,B)||groundset[iter]==0)
                continue;

            //oracle_times++;
            double value=f_u(iter);

            if(value>=f_u_max)
            {
                f_u_max=value;
                u_star=iter;
            }
        }

    S_class S_star;
    double gamma_1=0.0;
    double gamma_2=2.0*(n+eps)*f_u_max;
    double gamma=gamma_2;

    while(gamma_2-gamma_1>eps_pie*f_u_max)
    {
        /*******initial T1 T2*******/
        vector<S_class> T;
        vector<multimap<double,Ai_node>> A;

        for(int i=0;i<ell;i++) {
            T.push_back(S_class());
            A.push_back(multimap<double,Ai_node>());
        }
        vector<int> available(node_num,0);//mark all nodes are available ie.,N_i
        //initial Ai

            for (int iter = 0; iter < node_num; iter++) {
                if (!node_budget_feasible(iter, B) || groundset[iter] == 0)
                    continue;

                Ai_node temp(iter);
                double value = f_u(iter);
                for (int j = 0; j < ell; j++) {
                    A[j].insert(pair<double, Ai_node>(value, temp));
                }
                available[iter] = 1;

            }

        //cout<<"gamma: "<<gamma<<endl;
        /*******initial T1 T2*******/
        bool flag=true;//whether execute Line 20
        while (true) {
            for (int j = 0; j < T.size(); j++) {
                T[j].max_marginal = -999999999.0;
                T[j].max_element = -1;
                //A=empty, then return element=-1 and marginal gain =0
                pair<int, double> temp = GetMax(oracle_times, eps, A[j], T[j], available, B);
                T[j].max_marginal = temp.second;
                T[j].max_element = temp.first;
            }
            double max_marginal = 0.0;
            int max_element=-1;
            int max_solution = -1;
            for (int j = 0; j < T.size(); j++) {
                if (T[j].max_marginal > max_marginal) {
                    max_solution = j;
                    max_marginal = T[j].max_marginal;
                    max_element = T[j].max_element;
                }
            }
            if (max_solution == -1 || max_marginal <= 0)//non element or marginal gain<=0
                break;
            //S<-S\cup {u}
            if(max_marginal>=gamma/B*contrast_cost[max_element])
            {
                if(T[max_solution].budget_feasible(max_element,B))
                {
                    if (u(e_GENSM) == 1) {
                        T[max_solution].add_element(max_marginal, max_element);
                    }
                }
                else
                {
                    gamma_1=gamma;
                    for(auto & i : T)
                    {
                        if(i.S_revenue>S_star.S_revenue)
                        {
                            S_star=i;
                        }
                    }
                    flag=false;
                    // break;
                }
            }
            available[max_element] = 0;//discarded or selected
        }
        if(flag)
        {
            gamma_2=gamma;
        }
        for(auto & i : T)
        {
            if(i.S_revenue>S_star.S_revenue)
            {
                S_star=i;
            }
        }
        gamma=(gamma_1+gamma_2)/2.0;
    }
    return S_star.S_revenue;
}
Result GENSM_matroid(double eps,double B)
{
    cout<<"GENSM_matroid & Budget: "<<B<<"---------start---------"<<endl;
    double beta=8.5;
    long long int oracle_times=0;
    if(u1(e_GENSM)==1) {
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
        return Result(S.S_revenue, 0.333, S.Set.size(), oracle_times);
    }

    vector<int> available_self(node_num,0);//mark all nodes are available ie.,N_i
    vector<int> available_alg4(node_num,0);//mark all nodes are available ie.,N_i    for(int q=0;q<product_types;q++) {
        for (int iter = 0; iter < node_num; iter++) {
            if(u2(e_GENSM)==1)
            {
                available_self[iter]=1;
            }
            else
            {
                available_alg4[iter]=1;
            }
        }

    double x=ALG4(eps,B,available_alg4,oracle_times);
    cout<<"oracle: "<<oracle_times<<endl;
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

//            if (!node_budget_feasible(iter, B)) {
//                available_self[iter] = 0;
//                continue;
//            }

            Ai_node temp(iter);
            double value = f_u(iter);
            for (int j = 0; j < ell; j++) {
                A[j].insert(pair<double, Ai_node>(value, temp));
            }
        }


    while (true) {
        for (int j = 0; j < S.size(); j++)
        {
            S[j].max_marginal = -999999999.0;
            S[j].max_element = -1;
            //A=empty, then return element=-1 and marginal gain =0
            pair<int, double> temp = Boost_matroid(A[j], S[j], available_self, B);
            S[j].max_marginal = temp.second;
            S[j].max_element = temp.first;

            //only used for calculate oracle times

                for (int iter = 0; iter < node_num; iter++)
                {
                    if (available_self[iter] == 0)
                        continue;

                    if(S[j].is_feasible(iter,max_image))
                        oracle_times++;
                    else
                        available_self[iter] = 0;
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
        if (max_solution == -1 || max_marginal <= 0)//non element or marginal gain<=0
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

    cout<<"GENSM_matroid ---------end--------- "<<endl<<endl;
    return Result(S_star.S_revenue,S_star.S_cost,S_star.Set.size(),oracle_times);
}
#endif //IMAGE_GENSM_SLOW_MATROID_H
