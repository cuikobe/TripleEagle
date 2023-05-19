#include "generate_data.h"
#include "read_data.h"
#include "time.h"
#include "BFM_NM.h"
int main(int argc,char *argv[]) {

    //generate_data_with_lable();
    //return 0;
    //generate_data();
    string::size_type pos1, pos2,posend;
    pos1=cost_text.find_last_of("/");
    pos2=cost_text.rfind("/",pos1-1);

    //posend=cost_text.find_last_not_of("/");
    string name1=cost_text.substr(pos2+1,pos1-pos2-1);
    //string name2=cost_text.substr(pos1+1,posend);
    //string result_name=name1+"_"+name2;
    string result_name="engine_"+name1;

    read();
    cal_similarity();

    //test
    /*
    cout<<"f(u)"<<endl;
    for(int iter=0;iter<node_num;iter++)
    {
        cout<<f_u(iter)<<endl;
    }
    cout<<"f(u)/cost"<<endl;
    double sum_density=0.0;
    for(int iter=0;iter<node_num;iter++)
    {
        cout<<f_u(iter)/contrast_cost[iter]<<endl;
        sum_density+=f_u(iter)/contrast_cost[iter];
    }
    sum_density/=node_num;
    cout<<"ave density:"<<sum_density<<endl;
    //*/

    //double B=atof(argv[1]);
    double eps=0.1;
//    cout<<"eps: "<<eps<<endl;
    double default_p=sqrt(2)-1;

    time_t nowtime;
    struct tm* p;;
    time(&nowtime);
    p = localtime(&nowtime);
    string outtext="./result/image_result_normalize"+to_string((int)ave_num)+"_"+result_name+"_"+to_string(p->tm_mon+1)+"."+to_string(p->tm_mday)+"_"+to_string(p->tm_hour)+"_"+to_string(p->tm_min)+"_"+to_string(p->tm_sec)+".txt";

    vector<Result> GENSM_ksystem_result;
    vector<Result> random_clock_ksystem_result;

    vector<Result> simultaneous_result;
    vector<Result> GENSM_knapsack_result;
    vector<Result> random_clock_knapsack_result;
    vector<Result> BFM_NM_result;

    double B_start=1.0;
    double B_end=3.01;
    double B_step=0.2;


//    ofstream rcout("RC.txt");
//    ofstream tenmout("TENM.txt");
    for(double B=B_start;B<=B_end;B+=B_step)
    {
        //GENSM_ksystem_result.push_back(GENSM_matroid(eps,B));
        //random_clock_ksystem_result.push_back(RandomClockMatroid(eps,B));
        //random_clock_knapsack_result.push_back(RandomClockArbitrary(eps,B));
        // GENSM_result.push_back(GENSM_matroid(eps,B));

        /*********single run***************/

       // GENSM_knapsack_result.push_back(GENSM(eps,B));

        simultaneous_result.push_back(Simultaneous(eps,B));
      //  random_clock_knapsack_result.push_back(RandomClockOrder(eps,B));
        //BFM_NM_result.push_back(BFM_NM(eps,B));
        /*********multiple run***************/

        /**********calculate average for random algorithm************/
//        double revenue_ave=0.0;
//        long long int oracle_ave=0;
//        int repeat_time=50;
//        for(int i=0;i<repeat_time;i++) {
//            Result temp_result = GENSM(eps,B);
//            revenue_ave+=temp_result.revenue;
//            oracle_ave+=temp_result.oracle;
//        }
//        Result pricing_result(revenue_ave/repeat_time,-1,-1,oracle_ave/repeat_time);
//        GENSM_knapsack_result.emplace_back(pricing_result);
        /**********calculate average for random algorithm************/

        /**********calculate average for random algorithm************/
//        double revenue_ave2=0.0;
//        long long int oracle_ave2=0;
//        int repeat_time2=50;
//        for(int i=0;i<repeat_time2;i++) {
//            Result temp_result2 = RandomClockArbitrary(eps,B);
//            revenue_ave2+=temp_result2.revenue;
//            oracle_ave2+=temp_result2.oracle;
//        }
//        Result pricing_result2(revenue_ave2/repeat_time2,-1,-1,oracle_ave2/repeat_time2);
//        random_clock_knapsack_result.emplace_back(pricing_result2);
        /**********calculate average for random algorithm************/

        /**********calculate average for random algorithm************/
//        double revenue_ave3=0.0;
//        long long int oracle_ave3=0;
//        int repeat_time3=50;
//        for(int i=0;i<repeat_time3;i++) {
//            Result temp_result3 = BFM_NM(eps,B);
//            revenue_ave3+=temp_result3.revenue;
//            oracle_ave3+=temp_result3.oracle;
//        }
//        Result pricing_result3(revenue_ave3/repeat_time3,-1,-1,oracle_ave3/repeat_time3);
//        BFM_NM_result.emplace_back(pricing_result3);
        /**********calculate average for random algorithm************/

        /**********calculate average for random algorithm************/
        int repeat_time4=50;

        vector<double> utility4;
        double sum_utility4=0.0;
        vector<long long int> oracle4;
        long long int sum_oracle4=0;

        for(int i=0;i<repeat_time4;i++) {
            Result temp_result4 = BFM_NM(eps,B);

            sum_utility4+=temp_result4.revenue;
            utility4.push_back(temp_result4.revenue);

            sum_oracle4+=temp_result4.oracle;
            oracle4.push_back(temp_result4.oracle);
        }
        double average_utility4=sum_utility4/repeat_time4;
//        cout<<"average utility: "<<average_utility4<<endl;
        long long int average_oracle4=sum_oracle4/repeat_time4;
//        cout<<"average oracle: "<<average_oracle4<<endl;

//        tenmout<<"Budget="<<B<<endl;
        double utility_deviation4=0.0;
        for(auto x:utility4)
        {
            double temp=x-average_utility4;
            utility_deviation4+=pow(temp,2);
//            cout<<x<<'\t';
//
//            tenmout<<x<<'\t';
        }
//        tenmout<<endl;
//
//        cout<<endl;
        utility_deviation4=sqrt(utility_deviation4/(repeat_time4-1));
//        cout<<"utility standard_deviation: "<<utility_deviation4<<endl;

        double oracle_deviation4=0.0;
        for(auto y:oracle4)
        {
            long long int temp=y-average_oracle4;
            oracle_deviation4+=pow(temp,2);
        }
        oracle_deviation4=sqrt(oracle_deviation4/(repeat_time4-1));
//        cout<<"oracle standard_deviation: "<<oracle_deviation4<<endl;

        Result pricing_result4(average_utility4,average_oracle4,utility_deviation4,oracle_deviation4,-1);
        BFM_NM_result.emplace_back(pricing_result4);
        /**********calculate average for random algorithm************/

        /**********calculate average for random algorithm************/
        int repeat_time5=50;

        vector<double> utility5;
        double sum_utility5=0.0;
        vector<long long int> oracle5;
        long long int sum_oracle5=0;

        for(int i=0;i<repeat_time5;i++) {
            Result temp_result5 = RandomClockArbitrary(eps,B);

            sum_utility5+=temp_result5.revenue;
            utility5.push_back(temp_result5.revenue);

            sum_oracle5+=temp_result5.oracle;
            oracle5.push_back(temp_result5.oracle);
        }
        double average_utility5=sum_utility5/repeat_time5;
//        cout<<"average utility: "<<average_utility5<<endl;
        long long int average_oracle5=sum_oracle5/repeat_time5;
//        cout<<"average oracle: "<<average_oracle5<<endl;

//        rcout<<"Budget="<<B<<endl;
        double utility_deviation5=0.0;
        for(auto x:utility5)
        {
            double temp=x-average_utility5;
            utility_deviation5+=pow(temp,2);
//            cout<<x<<'\t';
//
//            rcout<<x<<'\t';
        }
//        cout<<endl;
//        rcout<<endl;

        utility_deviation5=sqrt(utility_deviation5/(repeat_time5-1));
//        cout<<"utility standard_deviation: "<<utility_deviation5<<endl;

        double oracle_deviation5=0.0;
        for(auto y:oracle5)
        {
            long long int temp=y-average_oracle5;
            oracle_deviation5+=pow(temp,2);
        }
        oracle_deviation5=sqrt(oracle_deviation5/(repeat_time5-1));
//        cout<<"oracle standard_deviation: "<<oracle_deviation5<<endl;

        Result pricing_result5(average_utility5,average_oracle5,utility_deviation5,oracle_deviation5,-1);
        random_clock_knapsack_result.emplace_back(pricing_result5);
        /**********calculate average for random algorithm************/

//        int equal_or_better=0;
//        for(int i=0;i<utility4.size();i++)
//        {
//            if(utility4[i]>=utility5[i])
//                equal_or_better++;
//        }
//        cout<<"the number of times our algorithm achieves equal or better results: "<<equal_or_better<<endl;
//
//        int better=0;
//        for(int i=0;i<utility4.size();i++)
//        {
//            if(utility4[i]>utility5[i])
//                better++;
//        }
//        cout<<"the number of times our algorithm achieves better results: "<<better<<endl;
    }

    ofstream out(outtext);
    out<<"eps: "<<eps<<endl;
    out<<"Budget: "<<endl;
    for(double B=B_start;B<=B_end;B+=B_step)
    {
        out<<B<<"\t";
    }
    out<<endl;

    out<<"RCA"<<endl;
    out<<"revenue: "<<endl;
    for(auto &p:random_clock_knapsack_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto &p:random_clock_knapsack_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
//    out<<"deviation revenue: "<<endl;
//    for(auto &p:random_clock_knapsack_result)
//    {
//        out<<p.utility_deviation<<"\t";
//    }
//    out<<endl;
//    out<<"deviation oracle times: "<<endl;
//    for(auto &p:random_clock_knapsack_result)
//    {
//        out<<p.oracle_deviation<<"\t";
//    }
//    out<<endl;


/*    out<<"RandomClock with ksystem"<<endl;
    out<<"revenue: "<<endl;
    for(auto &p:random_clock_ksystem_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto &p:random_clock_ksystem_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;*/

    out<<"SIP "<<endl;
    out<<"revenue: "<<endl;
    for(auto &p:simultaneous_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto &p:simultaneous_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;

//    out<<"GENSM without knapsack"<<endl;
//    out<<"revenue: "<<endl;
//    for(auto &p:GENSM_knapsack_result)
//    {
//        out<<p.revenue<<"\t";
//    }
//    out<<endl;
//    out<<"oracle times: "<<endl;
//    for(auto &p:GENSM_knapsack_result)
//    {
//        out<<p.oracle<<"\t";
//    }
//    out<<endl;

    out<<"TENM"<<endl;
    out<<"revenue: "<<endl;
    for(auto &p:BFM_NM_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto &p:BFM_NM_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
//    out<<"deviation revenue: "<<endl;
//    for(auto &p:BFM_NM_result)
//    {
//        out<<p.utility_deviation<<"\t";
//    }
//    out<<endl;
//    out<<"deviation oracle times: "<<endl;
//    for(auto &p:BFM_NM_result)
//    {
//        out<<p.oracle_deviation<<"\t";
//    }
//    out<<endl;

 /*   out<<"GENSM ksystem"<<endl;
    out<<"revenue: "<<endl;
    for(auto &p:GENSM_ksystem_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto &p:GENSM_ksystem_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    */


    return 0;
}
