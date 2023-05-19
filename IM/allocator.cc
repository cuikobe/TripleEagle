#include "allocator.h"
#include "anyoption.h"
#include "memoryusage.h"
namespace _Cide
{
    vector<Result> twin_result;
    vector<Result> randomTM_result;
    vector<Result> IP_result;
    vector<Result> BFM_RAN_result;
    vector<Result> BFM_DET_result;
    vector<Result> BFM_NM_result;

    vector<Result> randomTM_average;
    allocator::allocator(AnyOption* opt1)
    {
        B_for_Twin=0.0;
        opt = opt1;
        delim = " \t";//分隔符是制表符TAB

        //读取算法的参数
        n = strToInt(opt->getValue("n"));//社交网络的结点数，即用户数
        m = strToInt(opt->getValue("m"));//社交网络的边数f
        node_num=n;
        edge_num=m;

        nrTopics = strToInt(opt->getValue("nrTopics"));
        nrCompanies = strToInt(opt->getValue("nrCompanies"));
        alpha = strToFloat(opt->getValue("alpha"));

        //used for the RMA algorithm by HK
        epsilon = strToFloat(opt->getValue("epsilon"));
        theta_0 = strToFloat(opt->getValue("theta_0"));
        lambda = strToFloat(opt->getValue("lambda"));

        graphT.clear();
        for(int i = 0; i < n; i++) // keeps p^z_{uv} for each topic z
            graphT.push_back(vector<int>());//全局变量，给图中全部n个结点各一个vector数组，用于保存入边的源结点

//        max_node=100;
//        cout<<max_node<<endl;
        cout<<"h: "<<nrCompanies<<endl;
        advList = new advertiserList();
        // create advertiser objects and initialize their TIC prob-vector probT, aligned with graphT
        for(int i = 0; i < nrCompanies; i++) //遍历所有广告商
        {
            advertiser *aa = new advertiser(i,nrTopics);//初始化广告商的对象
            advList->push_back(aa);//advList存储所有广告商的对象

            aa->maxCost = 0;
            aa->minCost = (float) n ;
            aa->seedUserCosts.resize(n);//最多n个结点（用户），所以重设为n
            for (int i = 0; i < n; i++)
                aa->probT.push_back(std::vector< float>());
        }
        readItemDistsFile();
        readTICGraph();
        //is_selected=std::vector<bool>(n,true);//初始化为true,表示该结点可选
        readIncentiveCosts(); // reading incentive costs as matrix nodes X ads
        node_cost.assign(advList->at(0)->seedUserCosts.begin(),advList->at(0)->seedUserCosts.end());

        sfmt_t sfmt_seed;
        //随机种子
        for(int i = 0; i < 12; i++)//随机种子
            sfmt_init_gen_rand(&sfmt_seed, i+1234);
        //generate the original RRset,then can all copies use it
        timCS = new _Cide::TimGraph(advList->at(0), epsilon, n, m);
        timCS->theta_old=theta_0;
        timCS->theta=theta_0;
        timCS->RRset_generation();
        cout<<"RRset Size: "<<theta_0<<endl;
        //IP();

        RandomTM_acc();
        IP_acc();

        BFM_RAN();
        BFM_DET();
    }
    bool cmp(pair<pair<int, int>, float> a, pair<pair<int, int>, float> b)
    {
        return a.second > b.second;
    }
    void PRICING(S_class &A,const vector<int>&C,float neta,vector<float>&pi,const int &except_element,long long int &oracle_times,vector<long long int> &ask_times)
    {
        for(auto &u:C)
        {
            if(u==except_element)
                continue;

            oracle_times++;

            float delta_u=A.marginal(Node(u,-1.0));
            pi[u]=min(pi[u],Budget*delta_u/(A.S_revenue+neta));

            ask_times[u]++;
            if(pi[u]>=node_cost[u])
                A.add_element_truth(delta_u,Node(u,pi[u]),pi[u]);
        }
    }
    void allocator::BFM_DET()
    {
        vector<long long int> ask_times(node_num,0);

        cout<<"TED & Budget: "<<Budget<<"---------start---------"<<endl;
        long long int oracle_times=0;

        float alpha=(sqrt(13.0)+1.0)/2.0;

        /***********INITIALIZE*************/
        vector<float> pi(node_num,Budget);
        vector<int> C_available(node_num,0);//store the node that accepts the max budget
        float f_u_max=-1.0;
        int v_star;
        for (int iter = 0; iter < node_num; iter++) {
            if (!node_budget_feasible(iter, Budget))
                continue;

            C_available[iter]=1;
            oracle_times++;
            float value = f_u(iter,timCS);
            if (value >= f_u_max) {
                f_u_max = value;
                v_star = iter;
            }
        }
        S_class K(advList->at(0),timCS);
        /***********INITIALIZE*************/

        for(int u=0;u<node_num;u++)
        {
            if(C_available[u]==1&&u!=v_star)
            {
                oracle_times++;

                float temp_marginal=K.marginal(Node(u,-1.0));
                pi[u]=min(pi[u],Budget*temp_marginal/(alpha*f_u_max));

                ask_times[u]++;

                if(pi[u]>=node_cost[u])
                    K.add_element_truth(temp_marginal,Node(u,pi[u]),pi[u]);
                C_available[u]=0;
                if(K.S_revenue>=f_u_max)
                    break;
            }
        }

        S_class S(advList->at(0),timCS);
        if(K.S_revenue>=f_u_max)
        {
            vector<int> set_C_cup_v;
            for(int u=0;u<node_num;u++)
            {
                if(C_available[u]==1&&u!=v_star)
                {
                    set_C_cup_v.push_back(u);
                }
            }
            set_C_cup_v.push_back(v_star);
            PRICING(K,set_C_cup_v,alpha*f_u_max,pi,-1,oracle_times,ask_times);

            for (auto e = K.Set.rbegin(); e != K.Set.rend(); ++e)
            {
                if(S.S_price+pi[(*e).index]>Budget)
                {
                    break;
                }
                else
                {
                    S.add_element_truth(0.0,(*e),pi[(*e).index]);
                }
            }
            S.S_revenue=S.f_S(advList->at(0),timCS);//calculate the revenue in the end
        }
        else
        {
            S.add_element_truth(f_u_max,Node(v_star,pi[v_star]),pi[v_star]);
        }

        cout<<"S*:"<<endl;
        cout<<"  revenue: "<<S.S_revenue<<" size: "<<S.Set.size()<<" price: "<<S.S_price<<endl;
        cout<<"  all nodes: "<<endl;

        for(const auto &p:S.Set)
            cout<<p.index<<'\t';
        cout<<endl;
        cout<<"oracle times: "<<oracle_times<<endl;

//        cout<<"real revenue: "<<S.f_S(advList->at(0),timCS)<<endl;
        cout<<"TED ---------end--------- "<<endl<<endl;
        BFM_DET_result.emplace_back(S.S_revenue,S.S_price,S.Set.size(),oracle_times);


//        vector<int> sum_ask(4,0);
//        for(int i=0;i<node_num;i++)
//        {
//            if(ask_times[i]>=3)
//            {
//                sum_ask[3]++;
//                //cout<<"this elements is access more than 3 times: "<<i<<endl;
//            }
//            else
//            {
//                sum_ask[ask_times[i]]++;
//            }
//        }
//        //cout<<"v_star is: "<<v_star<<endl;
//        for(int i=0;i<sum_ask.size()-1;i++)
//        {
//            cout<<"the number of access "<<i<<" times: "<<'\t'<<sum_ask[i]<<endl;
//        }
//        cout<<"the number of access more than 3 times: "<<'\t'<<sum_ask[3]<<endl;
    }
    void allocator::BFM_RAN()
    {
        vector<long long int> ask_times(node_num,0);

        cout<<"TER & Budget: "<<Budget<<"---------start---------"<<endl;
        long long int oracle_times=0;

        float alpha=(sqrt(13.0)+1.0)/2.0;

        /***********INITIALIZE*************/
        vector<float> pi(node_num,Budget);
        vector<int> C;//store the node that accepts the max budget
        float f_u_max=-1.0;
        int v_star;
        for (int iter = 0; iter < node_num; iter++) {
            if (!node_budget_feasible(iter, Budget))
                continue;

            C.push_back(iter);
            oracle_times++;
            float value = f_u(iter,timCS);
            if (value >= f_u_max) {
                f_u_max = value;
                v_star = iter;
            }
        }
        /***********INITIALIZE*************/
        S_class A(advList->at(0),timCS);
        PRICING(A,C,alpha*f_u_max,pi,v_star,oracle_times,ask_times);

        S_class S(advList->at(0),timCS);
        if(A.S_revenue>=f_u_max)
        {
            vector<int> set_v;
            set_v.push_back(v_star);
            PRICING(A,set_v,alpha*f_u_max,pi,-1,oracle_times,ask_times);

            for (auto e = A.Set.rbegin(); e != A.Set.rend(); ++e)
            {
                if(S.S_price+pi[(*e).index]>Budget)
                {
                    break;
                }
                else
                {
                    S.add_element_truth(0.0,(*e),pi[(*e).index]);
                }
            }
            S.S_revenue=S.f_S(advList->at(0),timCS);//calculate the revenue in the end
        }
        else
        {
//            cout<<"randomness exists !"<<endl;

            uniform_real_distribution<float> dis(0.0,1.0);
            float Z=dis(e_BFM_RAN);
            if(Z<=  (alpha/(alpha+2.0))  )
            {
                S.add_element_truth(f_u_max,Node(v_star,pi[v_star]),pi[v_star]);
            }
            else
            {
                pi[v_star]=Budget-A.S_price;

                ask_times[v_star]++;

                if(pi[v_star]>=node_cost[v_star])
                {
                    A.add_element_truth(A.marginal(Node(v_star,pi[v_star])),Node(v_star,pi[v_star]),pi[v_star]);
                }
                S=A;
            }
        }
        cout<<"S*:"<<endl;
        cout<<"  revenue: "<<S.S_revenue<<" size: "<<S.Set.size()<<" price: "<<S.S_price<<endl;
        cout<<"  all nodes: "<<endl;

        for(const auto &p:S.Set)
            cout<<p.index<<'\t';
        cout<<endl;
        cout<<"oracle times: "<<oracle_times<<endl;

//        cout<<"real revenue: "<<S.f_S(advList->at(0),timCS)<<endl;
        cout<<"TER ---------end--------- "<<endl<<endl;
        BFM_RAN_result.emplace_back(S.S_revenue,S.S_price,S.Set.size(),oracle_times);

//        vector<int> sum_ask(4,0);
//        for(int i=0;i<node_num;i++)
//        {
//            if(ask_times[i]>=3)
//            {
//                sum_ask[3]++;
//            }
//            else
//            {
//                sum_ask[ask_times[i]]++;
//            }
//        }
//        for(int i=0;i<sum_ask.size()-1;i++)
//        {
//            cout<<"the number of access "<<i<<" times: "<<'\t'<<sum_ask[i]<<endl;
//        }
//        cout<<"the number of access more than 3 times: "<<'\t'<<sum_ask[3]<<endl;
    }
    void allocator::BFM_NM()
    {
        cout<<"BRM_NM & Budget: "<<Budget<<"---------start---------"<<endl;
        long long int oracle_times=0;

        S_class S_star;
        cout<<"S*:"<<endl;
        cout<<"  revenue: "<<S_star.S_revenue<<" size: "<<S_star.Set.size()<<" price: "<<S_star.S_price<<endl;
        cout<<"  all nodes: "<<endl;

        for(const auto &p:S_star.Set)
            cout<<p.index<<'\t';
        cout<<endl;
        cout<<"oracle times: "<<oracle_times<<endl;

        cout<<"real revenue: "<<S_star.f_S(advList->at(0),timCS)<<endl;

        cout<<"BRM_NM ---------end--------- "<<endl<<endl;
        BFM_NM_result.emplace_back(S_star.S_revenue,S_star.S_price,S_star.Set.size(),oracle_times);
    }
    void allocator::IP_acc()
    {
        vector<long long int> ask_times(node_num,0);

        cout<<"IP & Budget: "<<Budget<<"---------start---------"<<endl;

        long long int oracle_times=0;
        vector<int> A(node_num,0);//mark all nodes are available
        float f_u_max=-1.0;
        int u_star;

        for (int iter = 0; iter < node_num; iter++) {
            if (!node_budget_feasible(iter, Budget))
                continue;

            oracle_times++;
            float value = f_u(iter,timCS);
            if (value >= f_u_max) {
                f_u_max = value;
                u_star = iter;
            }
            A[iter] = 1;//single element with feasible budget is set to 1
        }

        vector<S_class> S;
        S.emplace_back(advList->at(0),timCS);
        S.emplace_back(advList->at(0),timCS);
        S[1].add_element_truth(f_u_max,Node(u_star,Budget),Budget);
//        cout<<"element star: "<<u_star<<" revenue: "<<f_u_max<<endl;
        float OPT=f_u_max;

        int t=1;
        vector<float> pi(node_num,Budget);
        while(true)
        {
            t++;
            OPT*=2.0;

            //int now_index=t%2;
            S[t%2].clear(advList->at(0),timCS);
            vector<int> now_avaliable(A);//mark all nodes are available now
            multimap<double,int> Boost_array;
            for (int iter = 0; iter < node_num; iter++) {
                //the element in S_{t-1} is not available
                if (S[(t + 1) % 2].selected[iter] == 1) {
                    now_avaliable[iter] = 0;
                }
                else
                {
                    Boost_array.insert(pair<double, int>(f_u(iter,timCS), iter));
                }
            }

//            cout<<"now OPT: "<<OPT<<endl;
            bool avaliable_empty=false;
            while (S[t%2].S_revenue<OPT)
            {
                pair<int,double> temp=Boost(Boost_array, S[t%2], now_avaliable, Budget);

//                if(Boost_array.size()%1000==0)
//                    cout<<Boost_array.size()<<endl;

                S[t%2].max_marginal=temp.second;
                S[t%2].max_element=temp.first;
                //only used for calculate oracle times
                for (int iter = 0; iter < node_num; iter++)
                {
                    if (now_avaliable[iter] == 0)
                        continue;
                    oracle_times++;
                }
/*                for(int e=0;e<node_num;e++)
                {
                    if(now_avaliable[e]==1)
                    {
                        oracle_times++;
                        float temp_marginal=S[t%2].marginal(Node(e,-1));
                        if(temp_marginal>S[t%2].max_marginal)
                        {
                            S[t%2].max_marginal=temp_marginal;
                            S[t%2].max_element=e;
                        }
                    }
                }*/

                //if (S[t%2].max_element == -1||S[t%2].max_marginal<=0)//non element or marginal gain<=0
                if (S[t%2].max_element == -1)//ONLY QUIT FOR NO ELEMENTS
                {
                    avaliable_empty=true;
                    break;
                }

                //save time
/*                cout<<S[t%2].max_marginal<<endl;
                if (S[t%2].max_marginal<=0.0)
                {
                    cout<<"the maximum marginal gain is 0.0"<<endl;
                    int remain_element=0;
                    for (int iter = 0; iter < node_num; iter++)
                    {
                        if (now_avaliable[iter] == 0)
                            continue;
                        remain_element++;
                    }
                    oracle_times+=(remain_element+1)*remain_element/2;
                    avaliable_empty=true;
                    break;
                }*/


                //S<-S\cup {u}
                pi[S[t%2].max_element]=min(pi[S[t%2].max_element],S[t%2].max_marginal*Budget/OPT);
                float price=pi[S[t%2].max_element];
                //cout<<"element: "<<S[t%2].max_element<<" price: "<<price<<" marginal: "<<S[t%2].max_marginal<<endl;

                ask_times[S[t%2].max_element]++;

                if(price>=node_cost[S[t%2].max_element]) {
                    S[t % 2].add_element_truth(S[t % 2].max_marginal, S[t % 2].max_element, price);
                }
                else
                    A[S[t%2].max_element]=0;//discarded the node who doesn't accept our price
                now_avaliable[S[t%2].max_element] = 0;//discarded or selected
            }

//            cout<<"t: "<<t<<endl;
//            for (const auto &it:S)
//            {
//                cout << "S: " << endl;
//                cout << "  revenue: " << it.S_revenue << " size: " << it.Set.size() << " price: " << it.S_price << endl;
//                cout << "  all nodes: " << endl;
//                for (const auto &p:it.Set)
//                    cout << p.index << '\t';
//                cout << endl;
//            }

            if(avaliable_empty)
                break;
        }
        if(S[(t+1)%2].S_price>Budget) {
//            cout << "The price of the W1 is more than Budget !" << endl;
            int j_index = S[(t + 1) % 2].Set.back().index;
            float j_price=pi[j_index];

            float temp_marginal = S[t % 2].marginal(Node(j_index, -1));
            pi[j_index] = min(pi[j_index], temp_marginal * Budget / OPT);
            double price = pi[j_index];

            ask_times[j_index]++;

            if (price >= node_cost[j_index]) {
                S[t % 2].add_element_truth(temp_marginal, Node(j_index, -1.0), price);
            }

            S[(t+1)%2].Set.pop_back();
            S[(t+1)%2].S_revenue=S[(t+1)%2].f_S(advList->at(0),timCS);
            S[(t+1)%2].S_cost-=j_price;
        }
        //Maximize-Value
        S_class W2(advList->at(0),timCS);
        for(auto e:S[t % 2].Set)
        {
            if(W2.S_price+pi[e.index]>Budget)
            {
                break;
            }
            else
            {
                W2.add_element_truth(0.0,e,pi[e.index]);
            }
            W2.S_revenue=W2.f_S(advList->at(0),timCS);
        }
        S_class W3=W2;
        //cout<<"W2 price: "<<W2.S_price<<endl;
        //cout<<"initialize W3 price: "<<W3.S_price<<endl;
        for(auto e:S[(t+1) % 2].Set)
        {
            if(W3.S_price+pi[e.index]>Budget)
            {
                break;
            }
            else
            {
                W3.add_element_truth(0.0,e,pi[e.index]);
            }
            W3.S_revenue=W3.f_S(advList->at(0),timCS);
        }
        S_class S_star;
        if(S[(t+1) % 2].S_revenue>W3.S_revenue)
        {
            S_star=S[(t+1) % 2];
        }
        else
        {
            S_star=W3;
        }

        cout<<"S*:"<<endl;
        cout<<"  revenue: "<<S_star.S_revenue<<" size: "<<S_star.Set.size()<<" price: "<<S_star.S_price<<endl;
        cout<<"  all nodes: "<<endl;

        for(const auto &p:S_star.Set)
            cout<<p.index<<'\t';
        cout<<endl;
        cout<<"oracle times: "<<oracle_times<<endl;

//        cout<<"real revenue: "<<S_star.f_S(advList->at(0),timCS)<<endl;
        cout<<"IP ---------end--------- "<<endl<<endl;
        IP_result.emplace_back(S_star.S_revenue,S_star.S_price,S_star.Set.size(),oracle_times);

//        time_t nowtime;
//        struct tm* p;;
//        time(&nowtime);
//        p = localtime(&nowtime);
//        string outtext="./result/ask_times_result_Budget="+to_string(Budget)+"_"+to_string(p->tm_mon+1)+"."+to_string(p->tm_mday)+"_"+to_string(p->tm_hour)+"_"+to_string(p->tm_min)+"_"+to_string(p->tm_sec)+".txt";
//        ofstream out(outtext);
//        vector<int> sum_ask(4,0);
//        for(int i=0;i<node_num;i++)
//        {
//            if(ask_times[i]>=3)
//            {
//                sum_ask[3]++;
//            }
//            else
//            {
//                sum_ask[ask_times[i]]++;
//            }
////            out<<i<<'\t'<<ask_times[i]<<endl;
//        }
////        out.close();


//        vector<int> sum_ask(4,0);
//        for(int i=0;i<node_num;i++)
//        {
//            if(ask_times[i]==2)
//                cout<<"this elements is accessed 2 times: "<<i<<endl;
//
//
//            if(ask_times[i]>=3)
//            {
//                sum_ask[3]++;
//            }
//            else
//            {
//                sum_ask[ask_times[i]]++;
//            }
//        }
//        for(int i=0;i<sum_ask.size()-1;i++)
//        {
//            cout<<"the number of access "<<i<<" times: "<<'\t'<<sum_ask[i]<<endl;
//        }
//        cout<<"the number of access more than 3 times: "<<'\t'<<sum_ask[3]<<endl;

    }

    void allocator::IP()
    {
        cout<<"IP & Budget: "<<Budget<<"---------start---------"<<endl;

        long long int oracle_times=0;
        vector<int> A(node_num,0);//mark all nodes are available
        float f_u_max=-1.0;
        int u_star;

        for (int iter = 0; iter < node_num; iter++) {
            if (!node_budget_feasible(iter, Budget))
                continue;

            oracle_times++;
            float value = f_u(iter,timCS);
            if (value >= f_u_max) {
                f_u_max = value;
                u_star = iter;
            }
            A[iter] = 1;//single element with feasible budget is set to 1
        }

        vector<S_class> S;
        S.emplace_back(advList->at(0),timCS);
        S.emplace_back(advList->at(0),timCS);
        S[1].add_element_truth(f_u_max,Node(u_star,Budget),Budget);
        cout<<"element star: "<<u_star<<" revenue: "<<f_u_max<<endl;
        float OPT=f_u_max;

        int t=1;
        vector<float> pi(node_num,Budget);
        while(true)
        {
            t++;
            OPT*=2.0;

            //int now_index=t%2;
            S[t%2].clear(advList->at(0),timCS);
            vector<int> now_avaliable(A);//mark all nodes are available now
            for (int iter = 0; iter < node_num; iter++) {
                //the element in S_{t-1} is not available
                if (S[(t + 1) % 2].selected[iter] == 1) {
                    now_avaliable[iter] = 0;
                }
            }

            cout<<"now OPT: "<<OPT<<endl;
            bool avaliable_empty=false;
            while (S[t%2].S_revenue<OPT)
            {
                S[t%2].max_marginal=-999999999.0;
                S[t%2].max_element=-1;
                for(int e=0;e<node_num;e++)
                {
                    if(now_avaliable[e]==1)
                    {
                        oracle_times++;
                        float temp_marginal=S[t%2].marginal(Node(e,-1));
                        if(temp_marginal>S[t%2].max_marginal)
                        {
                            S[t%2].max_marginal=temp_marginal;
                            S[t%2].max_element=e;
                        }
                    }
                }
                if (S[t%2].max_element == -1||S[t%2].max_marginal<=0)//non element or marginal gain<=0
                {
                    avaliable_empty=true;
                    break;
                }
                //S<-S\cup {u}
                pi[S[t%2].max_element]=min(pi[S[t%2].max_element],S[t%2].max_marginal*Budget/OPT);
                double price=pi[S[t%2].max_element];
                //cout<<"element: "<<S[t%2].max_element<<" price: "<<price<<" marginal: "<<S[t%2].max_marginal<<endl;
                if(price>=node_cost[S[t%2].max_element])
                    S[t % 2].add_element_truth(S[t%2].max_marginal, S[t%2].max_element, price);
                else
                    A[S[t%2].max_element]=0;//discarded the node who doesn't accept our price
                now_avaliable[S[t%2].max_element] = 0;//discarded or selected
            }

            cout<<"t: "<<t<<endl;
                for (const auto &it:S)
                {
                    cout << "S: " << endl;
                    cout << "  revenue: " << it.S_revenue << " size: " << it.Set.size() << " price: " << it.S_price << endl;
                    cout << "  all nodes: " << endl;
                    for (const auto &p:it.Set)
                        cout << p.index << '\t';
                    cout << endl;
                }
            if(avaliable_empty)
                break;
        }
        if(S[(t+1)%2].S_price>Budget)
        {
            cout << "The price of the W1 is more than Budget !" << endl;
            int j_index = S[(t + 1) % 2].Set.back().index;
            float j_price=pi[j_index];

            float temp_marginal=S[t%2].marginal(Node(j_index,-1));
            pi[j_index]=min(pi[j_index],temp_marginal*Budget/OPT);
            double price=pi[j_index];
            if(price>=node_cost[j_index])
                S[t % 2].add_element_truth(temp_marginal, Node(j_index,-1.0), price);

            S[(t+1)%2].Set.pop_back();
            S[(t+1)%2].S_revenue=S[(t+1)%2].f_S(advList->at(0),timCS);
            S[(t+1)%2].S_cost-=j_price;
        }
        //Maximize-Value
        S_class W2(advList->at(0),timCS);
        for(auto e:S[t % 2].Set)
        {
            if(W2.S_price+pi[e.index]>Budget)
            {
                break;
            }
            else
            {
                W2.add_element_truth(0.0,e,pi[e.index]);
            }
            W2.S_revenue=W2.f_S(advList->at(0),timCS);
        }
        S_class W3=W2;
        cout<<"W2 price: "<<W2.S_price<<endl;
        cout<<"initialize W3 price: "<<W3.S_price<<endl;
        for(auto e:S[(t+1) % 2].Set)
        {
            if(W3.S_price+pi[e.index]>Budget)
            {
                break;
            }
            else
            {
                W3.add_element_truth(0.0,e,pi[e.index]);
            }
            W3.S_revenue=W3.f_S(advList->at(0),timCS);
        }
        S_class S_star;
        if(S[(t+1) % 2].S_revenue>W3.S_revenue)
        {
            S_star=S[(t+1) % 2];
        }
        else
        {
            S_star=W3;
        }

        cout<<"S*:"<<endl;
        cout<<"  revenue: "<<S_star.S_revenue<<" size: "<<S_star.Set.size()<<" price: "<<S_star.S_price<<endl;
        cout<<"  all nodes: "<<endl;

        for(const auto &p:S_star.Set)
            cout<<p.index<<'\t';
        cout<<endl;
        cout<<"oracle times: "<<oracle_times<<endl;

        cout<<"real revenue: "<<S_star.f_S(advList->at(0),timCS)<<endl;

        cout<<"IP ---------end--------- "<<endl<<endl;
        IP_result.emplace_back(S_star.S_revenue,S_star.S_price,S_star.Set.size(),oracle_times);
    }

//    void allocator::RandomTM()
//    {
//        cout<<"RandomTM & Budget: "<<Budget<<"---------start---------"<<endl;
//        long long int oracle_times=0;
//        //vector<int> feasible_node(node_num,0);//mark all nodes are available
//        vector<Node> N;
//        float f_u_max=-1.0;
//        int i_star;
//        for (int iter = 0; iter < node_num; iter++) {
//            if (node_budget_feasible(iter, Budget))
//            {
//                N.emplace_back(iter,-1);
//
//                oracle_times++;
//                float value = f_u(iter,timCS);
//                if (value > f_u_max) {
//                    f_u_max = value;
//                    i_star = iter;
//                }
//            }
//        }
//        //random
//        float gamma=0.5;
//        float p=1.0-(gamma+1.0)/(gamma+2.0);
//
//        bernoulli_distribution u(p);
//        if(u(e_RandomTM))
//        {
//            cout << "return single element ! " << endl << endl;
//
//            randomTM_result.emplace_back(f_u_max,node_cost[i_star],1,oracle_times);
//            cout<<"S*:"<<endl;
//            cout<<"  revenue: "<<f_u_max<<" size: "<<1<<" cost: "<<node_cost[i_star]<<endl;
//            cout<<"  all nodes: "<<endl;
//            cout<<i_star<<endl;
//            cout<<"oracle times: "<<oracle_times<<endl;
//
//            S_class S(advList->at(0),timCS);
//            S.add_element(f_u_max,Node(i_star,-1));
//            cout<<"real F(S): "<<S.f_S(advList->at(0),timCS)<<endl;
//            cout<<"RandomTM ---------end--------- "<<endl<<endl;
//
//            return ;
//        }
//        S_class S(advList->at(0),timCS);
//        for(auto &e:N)
//        {
//            oracle_times++;
//            float temp_marginal=S.marginal(e);
//            if(node_cost[e.index]<=gamma*Budget*temp_marginal/S.S_revenue)
//            {
//                S.add_element(temp_marginal,e);
//            }
//            else break;
//        }
//
//        randomTM_result.emplace_back(S.S_revenue,S.S_price,S.Set.size(),oracle_times);
//        cout<<"S*:"<<endl;
//        cout<<"  revenue: "<<S.S_revenue<<" size: "<<S.Set.size()<<" price: "<<S.S_price<<endl;
//        cout<<"  all nodes: "<<endl;
//
//        for(const auto &e:S.Set)
//            cout<<e.index<<'\t';
//        cout<<endl;
//
//        cout<<"oracle times: "<<oracle_times<<endl;
//
//        cout<<"real F(S): "<<S.f_S(advList->at(0),timCS)<<endl;
//        cout<<"RandomTM ---------end--------- "<<endl;
//    }

    void allocator::RandomTM_acc()
    {
        cout<<"RTM & Budget: "<<Budget<<"---------start---------"<<endl;
        long long int oracle_times=0;
        //vector<int> feasible_node(node_num,0);//mark all nodes are available
        vector<Node> N;
        float f_u_max=-1.0;
        int i_star;
        for (int iter = 0; iter < node_num; iter++) {
            if (node_budget_feasible(iter, Budget))
            {
                N.emplace_back(iter,-1);

                oracle_times++;
                float value = f_u(iter,timCS);
                if (value > f_u_max) {
                    f_u_max = value;
                    i_star = iter;
                }
            }
        }
        //random
        float gamma=0.5;
        float p=1.0-(gamma+1.0)/(gamma+2.0);

        bernoulli_distribution u(p);
        if(u(e_RandomTM))
        {
//            cout << "return single element ! " << endl << endl;

            randomTM_result.emplace_back(f_u_max,node_cost[i_star],1,oracle_times);
            cout<<"S*:"<<endl;
            cout<<"  revenue: "<<f_u_max<<" size: "<<1<<" price: "<<node_cost[i_star]<<endl;
            cout<<"  all nodes: "<<endl;
            cout<<i_star<<endl;
            cout<<"oracle times: "<<oracle_times<<endl;

            S_class S(advList->at(0),timCS);
            S.add_element(f_u_max,Node(i_star,-1));
//            cout<<"real F(S): "<<S.f_S(advList->at(0),timCS)<<endl;
            cout<<"RTM ---------end--------- "<<endl<<endl;

            return ;
        }
        S_class S(advList->at(0),timCS);
        vector<int> now_avaliable(node_num,1);//mark all nodes are available now
        multimap<double,int> Boost_array;
        for (int iter = 0; iter < node_num; iter++) {
                Boost_array.insert(pair<double, int>(f_u(iter,timCS)/node_cost[iter], iter));
        }
        bool avaliable_empty=false;
        while(!Boost_array.empty()) {
            pair<int, double> temp = Boost_density(Boost_array, S, now_avaliable, Budget);
            S.max_marginal = temp.second*node_cost[temp.first];//return density, so we need to multiply cost to get marginal gain
            S.max_element = temp.first;
            //only used for calculate oracle times
            for (int iter = 0; iter < node_num; iter++) {
                if (now_avaliable[iter] == 0)
                    continue;
                oracle_times++;
            }
            if (S.max_element == -1)//ONLY QUIT FOR NO ELEMENTS
            {
                avaliable_empty = true;
                break;
            }
            if (node_cost[S.max_element] <= gamma * Budget * S.max_marginal / S.S_revenue) {
                S.add_element(S.max_marginal, S.max_element);
            }
            else break;
        }
        randomTM_result.emplace_back(S.S_revenue,S.S_price,S.Set.size(),oracle_times);
        cout<<"S*:"<<endl;
        cout<<"  revenue: "<<S.S_revenue<<" size: "<<S.Set.size();
        //<<" price: "<<S.S_price<<endl;
        cout<<"  all nodes: "<<endl;

        for(const auto &e:S.Set)
            cout<<e.index<<'\t';
        cout<<endl;

        cout<<"oracle times: "<<oracle_times<<endl;

//        cout<<"real F(S): "<<S.f_S(advList->at(0),timCS)<<endl;
        cout<<"RTM ---------end--------- "<<endl;
    }
//    void allocator::TwinGreedy()
//    {
//        S_class S1(advList->at(0),timCS);
//        S_class S2(advList->at(0),timCS);
//
//        B_for_Twin= B_for_Twin * nrCompanies;
//
//        first_u1=true;
//        first_u2=true;
//
////        //S_1
////            //generate 每个广告商的RR集
////            copy_tim(timCS,S1.graph);
////        //S_2
////            //generate 每个广告商的RR集
////            copy_tim(timCS,S2.graph);
//
//        //TimGraph *tim;
//        S1.graph->node_aval.resize(n,true);
//        S2.graph->node_aval.resize(n,true);
//
//        int aval_num1=nrCompanies*n;
//        int aval_num2=nrCompanies*n;
//        //int random_num;
//
//        long long int mariginal_count=0;
//        long long int martoid_count=0;
//        int sum_node1=0;
//        int sum_node2=0;
//
//        float max_revenue1=-999999999;
//        pair<int,int> argmax_node1;
//        bool need_argmax1=true;
//        float max_revenue2=-999999999;
//        pair<int,int> argmax_node2;
//        bool need_argmax2=true;
//
//        bool S1full=false;
//        bool S2full=false;
//        time(&startTime);
//        while((aval_num1>0)||(aval_num2>0))
//        {
//            if(sum_node1<max_node)
//            {
//                if (need_argmax1)
//                {
//                    max_revenue1 = -999999999;
//                    if(aval_num1>0)
//                    {
//                            for (int j = 0; j < n; j++) {
//                                martoid_count++;
//                                if (S1.graph->node_aval[j]) {
//                                    mariginal_count++;
//
//                                    float temp = S1.marginal(Node(j,-1.0));
//                                    if (first_u1) {
//                                        temp += lambda * B_for_Twin;
//                                    }
//                                    if (temp > max_revenue1) {
//                                        max_revenue1 = temp;
//                                        argmax_node1.first = j;
//                                        argmax_node1.second = 0;
//                                    }
//                                }
//                            }
//
//                    }
//                }
//            }
//            else
//            {
//                S1full= true;
//                max_revenue1 = -999999999;
//            }
//            if(sum_node2<max_node)
//            {
//                if (need_argmax2)
//                {
//                    max_revenue2 = -999999999;
//                    if(aval_num2>0)
//                    {
//                            for (int j = 0; j < n; j++) {
//                                martoid_count++;
//                                if (S2.graph->node_aval[j]) {
//                                    mariginal_count++;
//
//                                    float temp = S2.marginal(Node(j,-1.0));
//                                    if (first_u2) {
//                                        temp += lambda * B_for_Twin;
//                                    }
//                                    if (temp > max_revenue2) {
//                                        max_revenue2 = temp;
//                                        argmax_node2.first = j;
//                                        argmax_node2.second = 0;
//                                    }
//                                }
//                            }
//                    }
//                }
//            }
//            else
//            {
//                S2full= true;
//                max_revenue2 = -999999999;
//            }
//            if(S1full&&S2full) break;
//            if((max_revenue1<=0)&&(max_revenue2<=0))
//            {
//                cout<<"quit for argmax<=0!"<<endl;
//                break;
//            }
//
//            if(max_revenue1>max_revenue2)
//            {
//                //S1.graph->opim_assign_best_node(argmax_node1.first);
//                S1.add_element(max_revenue1,Node(argmax_node1.first,-1.0));
//
//                sum_node1++;
//                first_u1=false;
//                //S1.graph->currentRevenue+=max_revenue1;
//                //cout<<"S1 choose ("<<argmax_node1.first<<" , "<<argmax_node1.second<<")"<<" currentRev"<<tim->currentRevenue<<endl;
//
//                S1.graph->node_aval[argmax_node1.first] = false;
//                    aval_num1--;
//
//
//                S2.graph->node_aval[argmax_node1.first] = false;
//                    aval_num2--;
//
//
//                need_argmax1=true;
//                need_argmax2=false;
//            }
//            else{
//                //S2.graph->opim_assign_best_node(argmax_node2.first);
//                S2.add_element(max_revenue2,Node(argmax_node2.first,-1.0));
//
//                sum_node2++;
//                first_u2=false;
//                //S2.graph->currentRevenue+=max_revenue2;
//                //cout<<"S2 choose ("<<argmax_node2.first<<" , "<<argmax_node2.second<<")"<<" currentRev"<<tim->currentRevenue<<endl;
//
//
//                S2.graph->node_aval[argmax_node2.first] = false;
//                    aval_num2--;
//
//                S1.graph->node_aval[argmax_node2.first] = false;
//                    aval_num1--;
//
//                need_argmax2=true;
//                need_argmax1=false;
//            }
//
//        }
//
//        float totalDuration = getRunningTime(startTime); // in seconds
//        float totalMemory = disp_mem_usage(); // in MB
//        cout << "总时间开销: " << totalDuration << " s 总内存开销: " << totalMemory <<" MB "<< endl;
//
//        // write results to master and adv specific files
//        float total_revenue1 = 0.0;
//        float total_revenue2 = 0.0;
//        int total_seedsize1=0;
//        int total_seedsize2=0;
//
//        cout << "writing the results to output files.." << endl;
//
//        for(int i = 0; i < nrCompanies; i++) //遍历每个广告商，输出到对应的txt文件中
//        {
////            total_revenue1 += S1.graph->currentRevenue;//累计整个分配的收益等
////            total_revenue2 += S2.graph->currentRevenue;//累计整个分配的收益等
//
//            total_seedsize1+=S1.graph->seedSet.size();
//            total_seedsize2+=S2.graph->seedSet.size();
//        }
//        float total1 = 0.0;
//        float total2 = 0.0;
//        total1=S1.S_revenue;
//        total2=S2.S_revenue;
//        cout<<"S1:"<<endl;
//        cout<<"f(S1): "<<total1<<endl;
//        cout<<"Seed Size: "<<total_seedsize1<<endl;
//        cout<<"S2:"<<endl;
//        cout<<"f(S2): "<<total2<<endl;
//        cout<<"Seed Size: "<<total_seedsize2<<endl;
//
//        if(total1>total2)
//        {
//            cout<<"S*:"<<endl;
//            cout<<"f(S*): "<<total1<<endl;
//            cout<<"Seed Size: "<<total_seedsize1<<endl;
//        }
//        else{
//            cout<<"S*:"<<endl;
//            cout<<"f(S*): "<<total2<<endl;
//            cout<<"Seed Size: "<<total_seedsize2<<endl;
//        }
//        cout<<"mariginal times: "<<mariginal_count<<endl;
//        cout<<"martoid times: "<<martoid_count<<endl;
//        cout<<"max node: "<<max_node<<endl;
//
//        cout<<"real f(S1):"<<S1.f_S(advList->at(0),timCS)<<endl;
//        cout<<"real f(S2):"<<S2.f_S(advList->at(0),timCS)<<endl;
//
//        twin_result.emplace_back(total1,-1.0,total_seedsize1,mariginal_count);
//    }

    void copy_tim(TimGraph *src, TimGraph *dest)
    {
        dest->theta=src->theta;
        //dest->theta_old=src->theta;
        //清空原vector数组并复制
        dest->hyperGT_adv.assign(src->hyperGT_adv.begin(),src->hyperGT_adv.end());
        dest->isCovered.assign(src->isCovered.begin(),src->isCovered.end());
        dest->hyperG_adv.assign(src->hyperG_adv.begin(),src->hyperG_adv.end());
        dest->hyper_degree.assign(src->hyper_degree.begin(),src->hyper_degree.end());
    }

    void allocator::readTICGraph() {
        
        string probGraphFile = opt->getValue("probGraphFile");//社交网络图，每行是一条边，和在对应主题上生效的概率
        cout << "Reading file " << probGraphFile << endl;
        ifstream myfile (probGraphFile.c_str(), ios::in);
        
        float *dists;
        float p;
        advertiser *advTemp;
        
        int nrEdges = 0;
        set<int> nodes; // kontrol amacli simdi
        
        if (myfile.is_open()) {
            while (! myfile.eof() )	{
                std::string line;
                getline (myfile,line);
                if (line.empty()) continue;
                nrEdges++;
                
                std::string::size_type pos = line.find_first_of(delim);
                int prevpos = 0;
                
                //first user
                string str = line.substr(prevpos, pos-prevpos);
                int u1 = strToInt(str);
                
                //second user
                prevpos = line.find_first_not_of(delim, pos);
                pos = line.find_first_of(delim, prevpos);
                int u2 = strToInt(line.substr(prevpos, pos-prevpos));
                
                if (u1 == u2)//如果边是结点自身到自身，说明有误，直接跳过
                    continue;
                
                graphT[u2].push_back(u1); //insert to the transposed graph，为每个目标结点保存入边的源结点
                
                // kontrol amacli
                nodes.insert(u1);//nodes是一个set集合，因此用来统计所有结点
                nodes.insert(u2);
                
                prevpos = line.find_first_not_of(delim, pos);
                
                str = line.substr(prevpos);
                dists = new float[nrTopics];
                stringTokenizer(str, dists, nrTopics, delim);//把在所有主题上的概率保存下来
                
                for(int i = 0; i < nrCompanies; i++)
                {
                    advTemp = advList->at(i);//临时广告商对象，遍历所有广告商
                    p = 0.0;
                    for(int j = 0; j < nrTopics; j++)
                        p += (dists[j] * advTemp->gamma[j]);//遍历所有主题，边的概率乘以广告商在主题上的分布概率才是真正的影响概率p^i_{uv}
                    //累加起来得到的p就是这只广告在这个边上可能会产生影响力的概率（没错，一只广告所有主题上概率的累加才是合理的）。
                    advTemp->probT[u2].push_back(p);//保存一个广告商，在选择一个用户激活宣传后，该用户通过不同的边产生影响力的概率
                }
            }
            
            myfile.close();
        }
        
        else
            cout << "Can't open friendship graph file " << probGraphFile << endl;
        
        cout << "Built transposed adj matrix from file " << endl;
        cout << "number of nodes " << nodes.size() << endl;
        cout << "number of edges " << nrEdges << endl;
        
    }
    
    void allocator::readItemDistsFile() {
        
        cout << "reading item distributions file " << endl;
        string itemDistsFile = opt->getValue("itemDistsFile");
        ifstream myfile(itemDistsFile.c_str(), ios::in);
        float *tokens;
        
        int advIndex = 0;
        
        if(myfile.is_open()) {
            while(!myfile.eof()) {
                std::string line;
                getline(myfile, line);//逐行读取
                if(line.empty())
                    continue;
                tokens = new  float[nrTopics];//主题数大小的float数组
                stringTokenizer(line, tokens, nrTopics, delim);//读取每行的广告商在不同主题上的概率分布
                advList->at(advIndex++)->setItemDist(tokens, nrTopics);
                if(advIndex >= nrCompanies)//超过了公司数，文件中剩余的概率分布就不需要了
                    break;
            }
            
            myfile.close();
        }
        
        else {
            cout << "problem opening the item distributions file, exiting... " << itemDistsFile <<  endl;
            exit(1);
        }
        
    }

    // this reads a cost file in the form of a cost matrix in the form of nodes X ads -- might be different for the scalability version
    void allocator::readIncentiveCosts() {
        //每行表示一个用户被不同广告商收买所需的费用
        //该费用通过不同的模型由影响力估得
        string readIncentiveCostsFile = opt->getValue("incentiveCostsFile");
        ifstream myfile(readIncentiveCostsFile.c_str(), ios::in);
        
        int lineIndex = 0;
        float *tokens;
        float tempCostToken = 0.0;


        if(myfile.is_open()) {
            while(!myfile.eof())
            {
                std::string line;
                getline(myfile, line);
                if(line.empty())
                    continue;
                tokens = new float[nrCompanies];
                stringTokenizer(line, tokens, nrCompanies, delim);//按广告商数目切割每行
                /*for (int i = 0; i < nrCompanies; i++)
                {
                    if(tokens[i] < 1.0)
                        tokens[i]= 1.1;
                }*/
                for(int i = 0; i < nrCompanies; i++)
                {
                    //采用不同的激励模型，alpha是控制大小的系数
                    //线性
                    if(string(opt->getValue("costFunctionType")).compare("l") == 0) { // linear
                        tempCostToken = tokens[i] * alpha;
                    }
                    //所有用户相同的，constan
                    else if(string(opt->getValue("costFunctionType")).compare("u") == 0) { // uniform -- reads uniform input
                        tempCostToken = tokens[i] * alpha;
                    }
                    //二次函数
                    else if(string(opt->getValue("costFunctionType")).compare("q") == 0) { // quadratic
                        tempCostToken = tokens[i] * tokens[i] * alpha;
                    }
                    //次线性，即Log
                    else if(string(opt->getValue("costFunctionType")).compare("s") == 0) { // sublinear
                        tempCostToken = log(tokens[i]) * alpha;
                        if(tempCostToken==0) tempCostToken=log(tokens[i]+0.1) * alpha;
                    }
                    //随机的
                    else if(string(opt->getValue("costFunctionType")).compare("r") == 0) { // random
                        
                        tempCostToken = tokens[i] * alpha; // daha sonra shuffle edilecek bu
                    }

                    //依次存储到list数据结构中，每次循环对应存储到当前广告商所对应于每个用户标号的cost
                    advList->at(i)->seedUserCosts[lineIndex] = tempCostToken;
                    //保存对于一个广告商来说，最贵的代言人和最便宜的代言人
                    if(advList->at(i)->maxCost < tempCostToken)
                        advList->at(i)->maxCost = tempCostToken;
                    
                    if(advList->at(i)->minCost > tempCostToken)
                        advList->at(i)->minCost = tempCostToken;

                }
                B_for_Twin+= tempCostToken;
                lineIndex++;
            }
            
            myfile.close();
            
            //            cout << " coost - kontrol : " << advList->at(0)->seedUserCosts[0] << endl;
            // for random cost function shuffle the vectors here
            //如果选择激励费用是随机的，则之前正常读入cost，然后在这里打乱
            if(string(opt->getValue("costFunctionType")).compare("r") == 0)
            { // random
                for (int i =0 ; i < nrCompanies; i++)
                {
                    //cout << "adv " << i << " cost kontrol " << advList->at(i)->maxCost << endl;
                    //                    cout << "randomness kontrol b4 : " << advList->at(i)->seedUserCosts[0] << endl;
			        std::srand(1);
                    std::random_shuffle(advList->at(i)->seedUserCosts.begin(), advList->at(i)->seedUserCosts.end());
                    //                    cout << "randomness kontrol after : " << advList->at(i)->seedUserCosts[0] << endl;
                }
            }
            //输出对于每个广告商来说，最贵的代言人和最便宜的代言人
            //*
            for (int i = 0; i < nrCompanies; i++)
            {
                cout << "max cost for adv " << i << " " << advList->at(i)->maxCost << endl;
                cout << "min cost for adv " << i << " " << advList->at(i)->minCost << endl;
            }//*/

        }
        
        else {
            cout << "problem opening the incentive costs file, exiting... " << readIncentiveCostsFile <<  endl;
            exit(1);
        }
        
    }
    
    allocator::~allocator()
    {
        cout << "destructor called " << endl;

//        delete advList;
//        for(int i = 0; i < nrCompanies; i++)
//        {
//            delete advList->at(i);
//        }
//        delete timCS;
    }
    pair<int,float> Boost(multimap<double,int> &A,S_class &Si,const vector<int> &available,const double &B)
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
            if(available[it->second]==0)
            {
                A.erase(it);
                continue;
            }
            //if the element is useful, then we compute its new weight, and let its update numbers +1
            //oracle_times++;
            double new_value=Si.marginal(it->second);
            //the value of element in multimap can be change directly, but the key can not be change, we need to erase it and re-insert it
            //it->second.l++;
            //if the value of the element diminishes not much, we can return it
            if(new_value>=old_value)
            {
                ai=it->second;
                m=new_value;
                break;
            }
                //or we re-insert now element and then check the next element
            else
            {
                int temp=it->second;
                //erase the element to update its weight
                A.erase(it);
                //if its update numbers is not greater than max, then we re-insert it into the queue
                //if(temp.l<=(log(node_num*product_types*2.0/eps)/log(1.0+eps))) {
                A.insert(pair<double, int>(new_value, temp));
                //}
            }
        }
        return pair<int,float>(ai,m);
    }
    pair<int,float> Boost_density(multimap<double,int> &A,S_class &Si,const vector<int> &available,const double &B)
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
            if(available[it->second]==0)
            {
                A.erase(it);
                continue;
            }
            //if the element is useful, then we compute its new weight, and let its update numbers +1
            //oracle_times++;
            double new_value=Si.marginal(it->second)/node_cost[it->second];
            //the value of element in multimap can be change directly, but the key can not be change, we need to erase it and re-insert it
            //it->second.l++;
            //if the value of the element diminishes not much, we can return it
            if(new_value>=old_value)
            {
                ai=it->second;
                m=new_value;
                break;
            }
                //or we re-insert now element and then check the next element
            else
            {
                int temp=it->second;
                //erase the element to update its weight
                A.erase(it);
                //if its update numbers is not greater than max, then we re-insert it into the queue
                //if(temp.l<=(log(node_num*product_types*2.0/eps)/log(1.0+eps))) {
                A.insert(pair<double, int>(new_value, temp));
                //}
            }
        }
        return pair<int,float>(ai,m);
    }
    float f_u(const int &node,const TimGraph *timTemp)
    {
        return (( float) node_num * (( float) timTemp->hyper_degree[node] / timTemp->theta));
    }
    bool node_budget_feasible(const int &node,float Budget)
    {
        if(node_cost[node]>Budget)
            return false;
        else
            return true;
    }
}
