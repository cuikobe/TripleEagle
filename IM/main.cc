#include <cstdlib>
#include <iostream>
#include "anyoption.h"
#include "allocator.h"

//test

AnyOption* readOptions(int argc, char* argv[]);

using namespace _Cide;

int main(int argc, char* argv[]) {

    AnyOption *opt = readOptions(argc,argv);//原作者自己写的方便读取option的函数

    string configFilename = opt->getValue("c");
    time_t nowtime;
    struct tm* p;;
    time(&nowtime);
    p = localtime(&nowtime);
    string outtext="./result/result_"+configFilename+"_"+to_string(p->tm_mon+1)+"."+to_string(p->tm_mday)+"_"+to_string(p->tm_hour)+"_"+to_string(p->tm_min)+"_"+to_string(p->tm_sec)+".txt";

    float B_start=5.0;
    float B_end=10.01;
    float B_step=0.5;

    e_RandomTM.seed(time(NULL));
    e_BFM_RAN.seed(time(NULL));

//    e_RandomTM.seed(11);
//    e_BFM_RAN.seed(2);

    for(Budget=B_start;Budget<=B_end;Budget+=B_step)
    {
        _Cide::allocator *alloc = new _Cide::allocator(opt);//论文核心算法的类
        delete alloc;//因为是指针型，所以要delete
    }


    ofstream out(outtext);
    out<<"eps: "<<strToFloat(opt->getValue("epsilon"))<<endl;
    out<<"Budget: "<<endl;
    for(double B=B_start;B<=B_end;B+=B_step)
    {
        out<<B<<"\t";
    }
    out<<endl;

//    out<<"RandomTM_average"<<endl;
//    out<<"revenue: "<<endl;
//    for(auto p:randomTM_average)
//    {
//        out<<p.revenue<<"\t";
//    }
//    out<<endl;
//    out<<"oracle times: "<<endl;
//    for(auto p:randomTM_average)
//    {
//        out<<p.oracle<<"\t";
//    }
//    out<<endl;

/*    out<<"BFM_RAN_average"<<endl;
    out<<"revenue: "<<endl;
    for(auto p:BFM_RAN_average)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:BFM_RAN_average)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;*/

    out<<"RTM"<<endl;
    out<<"revenue: "<<endl;
    for(auto p:randomTM_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:randomTM_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;

    out<<"IP"<<endl;
    out<<"revenue: "<<endl;
    for(auto p:IP_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:IP_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;

    out<<"TER"<<endl;
    out<<"revenue: "<<endl;
    for(auto p:BFM_RAN_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:BFM_RAN_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;

    out<<"TED"<<endl;
    out<<"revenue: "<<endl;
    for(auto p:BFM_DET_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"oracle times: "<<endl;
    for(auto p:BFM_DET_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
}

AnyOption* readOptions(int argc, char* argv[]) {//内荣全部是直接引用一个option库

    AnyOption *opt = new AnyOption();//首先创建一个option对象

    // ignore POSIX style options
    opt->noPOSIX();
    //设置报错信息
    opt->setVerbose(); /* print warnings about unknown options */
    opt->autoUsagePrint(true); /* print usage for bad options */

    opt->addUsage("");
    opt->addUsage("Usage: ");
    opt->addUsage("");
    opt->addUsage("-help Prints this help ");
    opt->addUsage(" -c <config_file> Specify config file ");
    opt->addUsage("");

    //设置option的名字，这样在optionfile或者命令行中就可以输入这些option的实际内容
    opt->setOption("probGraphFile");
    opt->setOption("n");
    opt->setOption("m");
    opt->setOption("itemDistsFile");
    opt->setOption("nrTopics");
    opt->setOption("nrCompanies");
    opt->setOption("costFunctionType");
    opt->setOption("alpha");
    opt->setOption("epsilon");
    opt->setOption("theta_0");
    opt->setOption("lambda");
    opt->setOption("max_node");
    opt->setOption("incentiveCostsFile"); // in the form of spread

    //设置option的名字，只能从命令行读取
    opt->setCommandFlag("help");//设置flag，与option不一样，这只是个标志，不会去读取option的内容
    opt->setCommandOption("c");//设置命令行中的option

    opt->processCommandArgs(argc, argv);//读取执行程序时的命令行

    if(opt->getFlag( "help" )) //如果是输入的flag：help，那就调出使用方法介绍
    {
        opt->printUsage();
        delete opt;
        exit(0);
    }

    const char* configFile = opt->getValue("c");//通过获取option：c的内容来获取参数文件的名称
    if (configFile == NULL) {//没输入则报错
        cout << "Config file not mentioned" << endl;
        opt->printUsage();
        delete opt;
        exit(0);
    }

    cout << "Config file : " << configFile << endl;//输出参数文件名称
    opt->processFile(configFile);//让option对象读取参数文件中的option
    cout << endl << endl;
    return opt;

}
