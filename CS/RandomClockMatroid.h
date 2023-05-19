//
// Created by CLAKERS on 2022/10/8.
//

#ifndef IMAGE_RANDOMCLOCKMATROID_H
#define IMAGE_RANDOMCLOCKMATROID_H
#include "read_data.h"
#include "RandomClock_arbitrary.h"
pair<int,double> Boost_matroid(multimap<double,Ai_node> &A,S_class &Si,const vector<int> &available,const double &B)
{
/********check this when pop a node*************/
    double m=0.0;//always record weight of aij
    int ai=-1;
    while(!A.empty())
    {
        //check the end node, i.e., the node with biggest weight
        auto it=A.end();
        it--;

        double old_value=it->first;

        /*******no need this?*****/
        //if the value of maximum element is less than 0, then all element has value less than 0 due to the diminish property of submodularity, so we break and return empty element.
/*        if(old_value<=0)
        {
            A.clear();
            ai.user=-1;
            ai.product=-1;
            m=0.0;
            break;
        }*/

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
        if(new_value>=old_value)
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
            //if(temp.l<=(log(node_num*product_types*2.0/eps)/log(1.0+eps))) {
            A.insert(pair<double, Ai_node>(new_value, temp));
            //}
        }
    }
    return pair<int,double>(ai,m);
}
Result RandomClockMatroid(double eps,double B)
{
    cout<<"RandomClockMatroid & Budget: "<<B<<"---------start---------"<<endl;

    //random
    double p=sqrt(5.0)/(sqrt(k+10.0)+sqrt(5.0));
    cout<<"original p: "<<p<<endl;
    //p=1.0;
    cout<<"now p: "<<p<<endl;
    bernoulli_distribution u(p);

    long long int oracle_times=0;
    vector<double> l_u(node_num,B);

    vector<int> D(node_num,0);//mark all nodes are available
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

        D[iter] = 1;//single element with feasible budget is set to 1
    }
    double gamma=f_u_max;

    vector<S_class> S;
    S.emplace_back();
    S.emplace_back();
    S[0].add_element(f_u_max,u_star);
    int t=0;

    while(true)
    {
        t++;
        S[t%2].clear();
        gamma*=2.0;

        //initialize C and A
        vector<int> C(node_num,0);//mark all nodes are available
        //bool empty=true;
        multimap<double,Ai_node> A;
        for (int iter = 0; iter < node_num; iter++) {
            //available in D and not in S_{t-1}
            if (D[iter] == 1 && S[(t + 1) % 2].selected[iter] == 0) {
                //empty=false;//has available element
                C[iter] = 1;

                Ai_node temp(iter);
                double value = f_u(iter);
                A.insert(pair<double, Ai_node>(value, temp));
            }
        }
//        if(empty)
//            break;

        cout<<"now gamma: "<<gamma<<endl;
        bool C_empty=false;
        while (true) {
            //A=empty, then return element=-1 and marginal gain =0
            pair<int, double> temp = Boost_matroid(A, S[t%2], C, B);
            double max_marginal=temp.second;
            int max_element=temp.first;

//                cout<<"node: ("<<max_element.product<<","<<max_element.user<<")"<<endl;
//                cout<<"marginal gain: "<<max_marginal<<endl;

            //only used for calculate oracle times
            for (int iter = 0; iter < node_num; iter++) {

                if (C[iter] == 0)
                    continue;
                if (S[t % 2].is_feasible(iter,max_image))
                    oracle_times++;
                else
                    C[iter] = 0;
            }
            if (max_element==-1)//non element or marginal gain<=0
            {
                C_empty=true;
                break;
            }
            //S<-S\cup {u}
            l_u[max_element]=min(l_u[max_element],max_marginal*B/gamma);

            double price=l_u[max_element];
            bool accept=false;

            //cout<<"price: "<<price<<endl;

            if(price>=contrast_cost[max_element])
            {
                if(S[t%2].S_revenue+max_marginal>gamma)
                    break;
                if(u(e2)==1) {
                    S[t % 2].add_element_truth(max_marginal, max_element, price);
                    accept=true;
                }
            }
            if(!accept)
                D[max_element]=0;
            C[max_element] = 0;//discarded or selected
        }
        if(C_empty)
            break;
/*        cout<<"t: "<<t<<endl;
        for(const auto& it:S) {
            cout<<"S: "<<endl;
            cout<<"  revenue: "<<it.S_revenue<<" size: "<<it.Set.size()<<" cost: "<<it.S_cost<<endl;
            cout<<"  all nodes: "<<endl;
        }*/

    }

    S_class S_star;
    for(const auto& it:S) {
        if (it.S_revenue >= S_star.S_revenue)
            S_star = it;
    }

    cout<<"S*:"<<endl;
    cout<<"  revenue: "<<S_star.S_revenue<<" size: "<<S_star.Set.size()<<" cost: "<<S_star.S_cost<<endl;
    cout<<"  all nodes: "<<endl;

    for(const auto &p:S_star.Set)
        cout<<p<<'\t';
    cout<<endl;

    cout<<"oracle times: "<<oracle_times<<endl;

    cout<<"real revenue: "<<S_star.f_S()<<endl;

    cout<<"RandomClockMatroid ---------end--------- "<<endl<<endl;
    return Result(S_star.S_revenue,S_star.S_cost,S_star.Set.size(),oracle_times);

}

#endif //IMAGE_RANDOMCLOCKMATROID_H
