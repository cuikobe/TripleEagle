//
// Created by CLAKERS on 2023/5/5.
//
#ifndef IMAGE_BFM_NM_H
#define IMAGE_BFM_NM_H
#include "GENSM_slow.h"
default_random_engine e_BFM(5);
uniform_real_distribution<double> dis(0.0,1.0);
void PRICING(S_class &A,const vector<int>&C,double neta,vector<double>&pi,const int &except_element,long long int &oracle_times,const double &Budget)
{
    for(auto &u:C)
    {
        if(u==except_element)
            continue;

        oracle_times++;

        double delta_u=A.marginal(u);
        pi[u]=min(pi[u],Budget*delta_u/(A.S_revenue+neta));
        if(pi[u]>=contrast_cost[u])
            A.add_element_truth(delta_u,u,pi[u]);
    }
}
Result BFM_NM(double eps,double B)
{
    cout<<"TENM & Budget: "<<B<<"---------start---------"<<endl;
    long long int oracle_times=0;
    float alpha=3.0;

    /***********INITIALIZE*************/
    vector<float> pi(node_num,B);
    vector<int> C;//store the node that accepts the max budget
    double f_u_max=-1.0;
    int v_star;
    for (int iter = 0; iter < node_num; iter++) {
        if (!node_budget_feasible(iter, B))
            continue;

        C.push_back(iter);
        oracle_times++;
        double value = f_u(iter);
        if (value >= f_u_max) {
            f_u_max = value;
            v_star = iter;
        }
    }
    /***********INITIALIZE*************/
    vector<S_class> A;
    A.emplace_back();
    A.emplace_back();

    /********reverse access*************/
//    for (auto u = C.rbegin(); u != C.rend(); ++u)
//    {
//        if(*u==v_star)
//            continue;
//        /************argmax f(u|Ai)***********/
//        int max_solution=-1;
//        double max_marginal=-999.0;
//        for(int iter=0;iter<A.size();iter++)
//        {
//            oracle_times++;
//            double temp_marginal=A[iter].marginal(*u);
//            if(temp_marginal>max_marginal)
//            {
//                max_marginal=temp_marginal;
//                max_solution=iter;
//            }
//        }
//        /************argmax f(u|Ai)***********/
//        double price=B*max_marginal/(A[max_solution].S_revenue+alpha*f_u_max);
//        pi[*u]=price;
//        if(price>=contrast_cost[*u])
//        {
//            A[max_solution].add_element_truth(max_marginal,*u,price);
//        }
//    }
    /********reverse access*************/

    /********order access*************/
    for(auto u:C)//order access
    {
        if(u==v_star)
            continue;
        /************argmax f(u|Ai)***********/
        int max_solution=-1;
        double max_marginal=-999.0;
        for(int iter=0;iter<A.size();iter++)
        {
            oracle_times++;
            double temp_marginal=A[iter].marginal(u);
            if(temp_marginal>max_marginal)
            {
                max_marginal=temp_marginal;
                max_solution=iter;
            }
        }
        /************argmax f(u|Ai)***********/
        double price=B*max_marginal/(A[max_solution].S_revenue+alpha*f_u_max);
        pi[u]=price;
        if(price>=contrast_cost[u])
        {
            A[max_solution].add_element_truth(max_marginal,u,price);
        }
    }
    /********order access*************/

    S_class S;
    if(max(A[0].S_revenue,A[1].S_revenue)>=f_u_max)
    {
        int u=v_star;
        /************argmax f(u|Ai)***********/
        int max_solution=-1;
        double max_marginal=-999.0;
        for(int iter=0;iter<A.size();iter++)
        {
            oracle_times++;
            double temp_marginal=A[iter].marginal(u);
            if(temp_marginal>max_marginal)
            {
                max_marginal=temp_marginal;
                max_solution=iter;
            }
        }
        /************argmax f(u|Ai)***********/
        double price=B*max_marginal/(A[max_solution].S_revenue+alpha*f_u_max);
        pi[u]=price;
        if(price>=contrast_cost[u])
        {
            A[max_solution].add_element_truth(max_marginal,u,price);
        }

        /************argmax f(X)***********/
        double max_revenue=-999.0;
        for(int iter=0;iter<A.size();iter++)
        {
            if(A[iter].S_revenue>max_revenue)
            {
                max_revenue=A[iter].S_revenue;
                max_solution=iter;
            }


//            cout<<"A"<<iter<<endl;
//            cout<<"  revenue: "<<A[iter].S_revenue<<" size: "<<A[iter].Set.size()<<" price: "<<A[iter].S_cost<<endl;
//            cout<<"  all nodes: "<<endl;
//
//            for(const auto &p:A[iter].Set)
//                cout<<p<<'\t';
//            cout<<endl;
//            cout<<"real revenue: "<<A[iter].f_S()<<endl;


        }
        /************argmax f(X)***********/
        for (auto e = A[max_solution].Set.rbegin(); e != A[max_solution].Set.rend(); ++e)
        {
            if(S.S_cost+pi[(*e)]>B)
            {
                break;
            }
            else
            {
                S.add_element_truth(0.0,(*e),pi[(*e)]);
            }
        }
        S.S_revenue=S.f_S();//calculate the revenue in the end
    }
    else
    {
//        cout<<"randomness exists !"<<endl;
        double Z=dis(e_BFM);
        if(Z<=0.5)
        {
            S.add_element_truth(f_u_max,v_star,pi[v_star]);
        }
        else
        {
            /************argmax f(X)***********/
            int max_solution=-1;
            double max_revenue=-999.0;
            for(int iter=0;iter<A.size();iter++)
            {
                if(A[iter].S_revenue>max_revenue)
                {
                    max_revenue=A[iter].S_revenue;
                    max_solution=iter;
                }
            }
            /************argmax f(X)***********/
            pi[v_star]=B-A[max_solution].S_cost;
            if(pi[v_star]>=contrast_cost[v_star]&&A[max_solution].marginal(v_star)>0)
            {
                oracle_times++;
                A[max_solution].add_element_truth(A[max_solution].marginal(v_star),v_star,pi[v_star]);
            }

            max_solution=-1;
            max_revenue=-999.0;
            for(int iter=0;iter<A.size();iter++)
            {
                if(A[iter].S_revenue>max_revenue)
                {
                    max_revenue=A[iter].S_revenue;
                    max_solution=iter;
                }
            }
            S=A[max_solution];
        }
    }
    cout<<"S*:"<<endl;
    cout<<"  revenue: "<<S.S_revenue<<" size: "<<S.Set.size()<<" price: "<<S.S_cost<<endl;
    cout<<"  all nodes: "<<endl;

    for(const auto &p:S.Set)
        cout<<p<<'\t';
    cout<<endl;
    cout<<"oracle times: "<<oracle_times<<endl;

//    cout<<"real revenue: "<<S.f_S()<<endl;
    cout<<"TENM ---------end--------- "<<endl<<endl;
    return Result(S.S_revenue,S.S_cost,S.Set.size(),oracle_times);
}


#endif //IMAGE_BFM_NM_H
