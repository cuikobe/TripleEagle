//
// Created by CLAKERS on 2022/10/3.
//

#ifndef IMAGE_SIMULTANEOUS_H
#define IMAGE_SIMULTANEOUS_H
#include "RandomClock_arbitrary.h"
S_class USM(const S_class &N,long long int &oracle_times)
{
    S_class X;
    S_class Y;
    Y.Set.assign(N.Set.begin(),N.Set.end());
    //Y.selected.assign(N.selected.begin(),N.selected.end());
    Y.S_cost=N.S_cost;
    Y.S_revenue=N.S_revenue;

    default_random_engine e(1234);
    for(auto u=N.Set.begin();u!=N.Set.end();u++)
        //for(auto u=N.Set.rbegin();u!=N.Set.rend();u++)
    {
        //calculate ai
        double f_Xi_1=X.S_revenue;
        //X.selected[*u]=1;
        X.Set.push_back(*u);

        oracle_times++;

        double f_Xi_1_and_u=X.f_S();
        X.Set.pop_back();
        //X.selected[*u]=0;
        double ai=f_Xi_1_and_u-f_Xi_1;

        //calculate bi
        double f_Yi_1=Y.S_revenue;
        //Y.selected[*u]=0;

        oracle_times++;

        double f_Yi_1_sub_u=Y.S_sub_u(*u);
        //Y.selected[*u]=1;
        double bi=f_Yi_1_sub_u-f_Yi_1;

        double ai1=max(ai,0.0);
        double bi1=max(bi,0.0);
        double probability=0.0;
        if(ai1==0&&bi1==0)
            probability=1.0;
        else
        {
            probability=ai1/(ai1+bi1);
        }
        bernoulli_distribution r(probability);

        if(bi>0) cout<<"bi>0 !!!"<<endl;
        if(r(e)==1)
        {
            //X.selected[*u]=1;
            X.Set.push_back(*u);
            X.S_revenue=f_Xi_1_and_u;

            X.S_cost+=N.node_price[*u];
            X.node_price[*u]=N.node_price[*u];
        }
        else
        {
            cout<<"ai<bi !!!"<<endl;
            //Y.selected[*u]=0;
            //delete u_i
            for(vector<int>::iterator p=Y.Set.begin();p!=Y.Set.end();)
            {
                if((*p)==(*u))
                {
                    p=Y.Set.erase(p);
                    break;
                }
                else{
                    p++;
                }
            }
            Y.S_revenue=f_Yi_1_sub_u;
            Y.S_cost-=N.node_price[*u];
        }
    }
    return X;
}
Result Simultaneous(double eps,double B)
{
    cout<<"SIP & Budget: "<<B<<"---------start---------"<<endl;

    long long int oracle_times=0;
    vector<double> pi(node_num,B);

    vector<int> A(node_num,0);//mark all nodes are available
    double f_u_max=-1.0;
    int u_star;
    for (int iter = 0; iter < node_num; iter++) {

//        if (!node_budget_feasible(iter, B))
//            continue;

        oracle_times++;
        double value = f_u(iter);

        if (value >= f_u_max) {
            f_u_max = value;
            u_star = iter;
        }

        A[iter] = 1;//single element with feasible budget is set to 1
    }

    double OPT=f_u_max;

    vector<vector<S_class>> S(2,vector<S_class>(2,S_class()));
    S[0][1].add_element(f_u_max,u_star);

    int t=0;
    while(true)
    {
        t++;
        for(auto &it:S[t%2])
            it.clear();

        OPT*=2.0;

        //initialize C and Boost_array
        vector<int> now_avaliable(node_num,0);//mark all nodes are available now
        //bool empty=true;

        vector<multimap<double,Ai_node>> Boost_array(2,multimap<double,Ai_node>());

        for (int iter = 0; iter < node_num; iter++) {
            //available in D and not in S_{t-1}
            if (A[iter] == 1 && S[(t + 1) % 2][0].selected[iter] == 0&& S[(t + 1) % 2][1].selected[iter] == 0) {
                //empty=false;//has available element
                now_avaliable[iter] = 1;

                Ai_node temp(iter);
                double value = f_u(iter);
                Boost_array[0].insert(pair<double, Ai_node>(value, temp));
                Boost_array[1].insert(pair<double, Ai_node>(value, temp));
            }
        }

//        cout<<"now OPT: "<<OPT<<endl;
        bool avaliable_empty=false;
        while (S[t%2][0].S_revenue<OPT&&S[t%2][1].S_revenue<OPT) {
            for (int j = 0; j < S[t%2].size(); j++)
            {
                S[t%2][j].max_marginal=-999999999.0;
                S[t%2][j].max_element=-1;

                //A=empty, then return element=-1 and marginal gain =0
                pair<int,double> temp=Boost(Boost_array[j], S[t%2][j], now_avaliable, B);
                S[t%2][j].max_marginal=temp.second;
                S[t%2][j].max_element=temp.first;

            }
            double max_marginal = 0.0;
            int max_element=-1;
            int max_solution = -1;
            for(int j=0;j<S.size();j++)
            {
                if (S[t%2][j].max_marginal > max_marginal) {
                    max_solution=j;
                    max_marginal=S[t%2][j].max_marginal;
                    max_element=S[t%2][j].max_element;
                }
            }
            //if (max_solution == -1||max_marginal<=0)//non element or marginal gain<=0
            if (max_solution == -1)//non element or marginal gain<=0
            {
                avaliable_empty=true;
                break;
            }

            //only used for calculate oracle times
            for (int iter = 0; iter < node_num; iter++)
            {
                if (now_avaliable[iter] == 0)
                    continue;
                oracle_times+=2;
            }

            //S<-S\cup {u}
            pi[max_element]=min(pi[max_element],max_marginal*B/OPT);

            double price=pi[max_element];
            if(price>=contrast_cost[max_element])
                S[t % 2][max_solution].add_element_truth(max_marginal, max_element, price);
            else
                A[max_element]=0;//discarded the node who doesn't accept our price
            now_avaliable[max_element] = 0;//discarded or selected
        }

//        cout<<"t: "<<t<<endl;
//        for(const auto& first:S) {
//            for (const auto &it:first) {
//                cout << "S: " << endl;
//                cout << "  revenue: " << it.S_revenue << " size: " << it.Set.size() << " cost: " << it.S_cost << endl;
//                cout << "  all nodes: " << endl;
//                for (const auto &p:it.Set)
//                    cout << p << '\t';
//                cout << endl;
//            }
//        }
        if(avaliable_empty)
            break;
    }
    //call USM
    S_class S_star;
    for(const auto& i:S) {
        for (const auto &j:i) {
            S_class temp=USM(j,oracle_times);
            if(temp.S_revenue>S_star.S_revenue)
                S_star=temp;
        }
    }
    for(const auto& it:S[(t + 1)%2]) {
        if (it.S_revenue >= S_star.S_revenue)
            S_star = it;
    }

    if(S_star.S_cost>B)
    {
//        cout<<"this happen !"<<endl;
//        cout<<"S revenue 1: "<<S_star.S_sub_u(S_star.Set.back())<<endl;
        S_star.S_cost-=S_star.node_price[S_star.Set.back()];
        S_star.Set.pop_back();
        S_star.S_revenue=S_star.f_S();
//        cout<<"S revenue 2: "<<S_star.S_revenue<<endl;
    }

    cout<<"S*:"<<endl;
    cout<<"  revenue: "<<S_star.S_revenue<<" size: "<<S_star.Set.size()<<" cost: "<<S_star.S_cost<<endl;
    cout<<"  all nodes: "<<endl;

    for(const auto &p:S_star.Set)
        cout<<p<<'\t';
    cout<<endl;

    cout<<"oracle times: "<<oracle_times<<endl;

//    cout<<"real revenue: "<<S_star.f_S()<<endl;

    cout<<"SIP ---------end--------- "<<endl<<endl;
    return Result(S_star.S_revenue,S_star.S_cost,S_star.Set.size(),oracle_times);

}
#endif //IMAGE_SIMULTANEOUS_H
