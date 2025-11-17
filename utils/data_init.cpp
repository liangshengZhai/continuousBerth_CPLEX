#include "data_init.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <cstring>

using namespace std;
// 简单的 mkdir -p 实现：逐级创建目录
static bool mkdir_p(const std::string& dirPath) {
    if (dirPath.empty()) return true;
    std::string path;
    for (size_t i = 0; i < dirPath.size(); ++i) {
        char c = dirPath[i];
        path.push_back(c);
        if (c == '/' || i == dirPath.size() - 1) {
            if (!path.empty() && path != "/" && path != "./") {
                struct stat st;
                if (stat(path.c_str(), &st) != 0) {
                    if (mkdir(path.c_str(), 0755) != 0 && errno != EEXIST) {
                        std::cerr << "创建目录失败: " << path << ", 错误: " << std::strerror(errno) << std::endl;
                        return false;
                    }
                } else if (!S_ISDIR(st.st_mode)) {
                    std::cerr << "路径存在但不是目录: " << path << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

// 获取 baseName 的目录部分（去掉最后一个'/'之后的内容）
static std::string dirname_of(const std::string& path) {
    auto pos = path.find_last_of('/') ;
    if (pos == std::string::npos) return std::string();
    return path.substr(0, pos);
}


// 读取模型参数（从 verify.cpp 提取）
ModelParams setParams(int numBerths,int numShips) {
    ModelParams params;

    params.numBerths = numBerths;
    params.numRows = 20;
    params.numSlotsPerRow = 24;
    params.numShips = numShips;
    params.planningHorizon = 168; // 一周（小时）
    params.numShipK=3;
    
    params.width = 50.0;
    params.relativeHeight = 0.5;
    params.alpha = 1.0;
    params.beta = 20000.0;

    std::random_device rd;  // 用于获取种子
    std::mt19937 gen(rd()); // 随机数生成器

    std::uniform_int_distribution<> weight(60*1000, 300*1000);
    std::uniform_int_distribution<> saveCos(900,1200);
    std::uniform_int_distribution<> R(15,20);
    std::uniform_int_distribution<> unloadSpeed(5*1000,11*1000);
    std::uniform_real_distribution<> angel(35 * M_PI / 180,50* M_PI / 180);
    std::uniform_int_distribution<> density(2,5);
    std::uniform_int_distribution<> Tau_s(0,params.planningHorizon);


    //随机生成到达时间
    params.arrivalTime.resize(params.numShips);
    for(int s = 0 ; s< params.numShips;s++){
        params.arrivalTime[s] = Tau_s(gen);
        cout<<"船"<<s<<"到达时间："<<params.arrivalTime[s]<<endl;
    }

    // 生成随机卸载速度
    params.unloadingSpeed.resize(params.numShips);
    for (int s = 0; s < params.numShips; ++s) {
        params.unloadingSpeed[s].resize(params.numBerths);
        for (int b = 0; b < params.numBerths; ++b) {
        params.unloadingSpeed[s][b].resize(params.numShipK);
            for(int k =0;k < params.numShipK;k++){
                params.unloadingSpeed[s][b][k] = unloadSpeed(gen);  // 为每个船舶-泊位组合生成随机速度
            }
        }
    }
    //随机生成 转运成本
    params.transshipmentCost.resize(params.numBerths);
    for(int b = 0; b < params.numBerths ;b++){
        params.transshipmentCost[b].resize(params.numRows);
        for(int r = 0; r < params.numRows;r++){
            params.transshipmentCost[b][r].resize(params.numSlotsPerRow);
            for(int v = 0 ; v < params.numSlotsPerRow;v++){
                //由于程序中 v 是从0开始的 不需要-1
                params.transshipmentCost[b][r][v] = R(gen) + 0.5*v;
            }
        }
    }
    //同一行对同一艘船存储成本一致
    params.storageCost.resize(params.numShips);
    for(int s = 0; s < params.numShips ; s++){
        params.storageCost[s].resize(params.numShipK);
        for(int k =0; k<params.numShipK;k++){
            params.storageCost[s][k].resize(params.numRows);
            for(int r = 0; r < params.numRows;r++){
                params.storageCost[s][k][r] = saveCos(gen);
            }   
        }  
    }

    // 初始化货物重量
    params.cargoWeight.resize(params.numShips);
    for(int s = 0 ; s < params.numShips;s++){
        params.cargoWeight[s] = weight(gen);
    }
   
    //初始化货物密度
    params.cargoDensity.resize(params.numShips);
    for(int s =0; s < params.numShips ; s++){
        params.cargoDensity[s].resize(params.numShipK);
        for(int k =0;k<params.numShipK;k++){
            params.cargoDensity[s][k] = density(gen);
        }
    }
    
    //初始化安息角
    params.maxResponseAngle.resize(params.numShips);
    for(int s = 0 ;s < params.numShips; s++){
        params.maxResponseAngle[s].resize(params.numShipK);
        for(int k =0; k< params.numShipK;k++){
            params.maxResponseAngle[s][k]=angel(gen);
        }
        
    }
    
    
    // 计算每艘船需要的槽数 n_s（根据文档公式）
    params.requiredSlots.resize(params.numShips);
    for (int s = 0; s < params.numShips; s++) {
        params.requiredSlots[s].resize(params.numShipK);
        for(int k = 0 ;k <params.numShipK;k++){
            double volume = params.cargoWeight[s] / (params.cargoDensity[s][k]*params.numShipK);
            double term = 4.0 / (pow(params.width, 3) * params.relativeHeight * 
                            (2 - params.relativeHeight) * tan(params.maxResponseAngle[s][k]));
            double inside = volume + (1.0/12.0) * pow(params.width, 3) * 
                        pow(params.relativeHeight, 2) * (3 - 2 * params.relativeHeight) * 
                        tan(params.maxResponseAngle[s][k]);

            params.requiredSlots[s][k] = ceil(term * inside);
        }
    }
    
    //槽数检验
    int totalRequiredSlots = 0;
    for (int s = 0; s < params.numShips; s++) {
        for(int k =0 ; k< params.numShipK;k++){
            cout<<params.requiredSlots[s][k]<<endl;
            totalRequiredSlots += params.requiredSlots[s][k];
        }
    }   
    cout<<"一共需要槽数：" <<totalRequiredSlots<<endl;
    int totalAvailableSlots = params.numRows * params.numSlotsPerRow;
    if (totalRequiredSlots > totalAvailableSlots) {
            cerr << "错误: 总需求槽数(" << totalRequiredSlots 
                << ")超过可用槽数(" << totalAvailableSlots << ")" << endl;
        exit(1);
    }

    return params;
}

// 将不同参数导出为多个 CSV 文件，Excel 可直接打开这些 CSV
void writeParamsToCSV(const ModelParams& params, const std::string& baseName) {
    // 确保输出目录存在
    {
        std::string dir = dirname_of(baseName);
        if (!dir.empty()) {
            mkdir_p(dir + "/");
        }
    }
    // 1) general
    {
        std::ofstream ofs(baseName + "_general.csv");
        ofs << "key,value\n";
        ofs << "numBerths," << params.numBerths << "\n";
        ofs << "numRows," << params.numRows << "\n";
        ofs << "numSlotsPerRow," << params.numSlotsPerRow << "\n";
        ofs << "numShips," << params.numShips << "\n";
        ofs << "planningHorizon," << params.planningHorizon << "\n";
        ofs << "numShipK," << params.numShipK << "\n";
        ofs << "width," << params.width << "\n";
        ofs << "relativeHeight," << params.relativeHeight << "\n";
        ofs << "alpha," << params.alpha << "\n";
        ofs << "beta," << params.beta << "\n";
        ofs.close();
    }

    // 2) arrivalTime
    {
        std::ofstream ofs(baseName + "_arrival.csv");
        ofs << "ship,arrivalTime\n";
        for (int s = 0; s < (int)params.arrivalTime.size(); ++s) {
            ofs << s << "," << params.arrivalTime[s] << "\n";
        }
        ofs.close();
    }

    // 3) cargoWeight
    {
        std::ofstream ofs(baseName + "_cargoWeight.csv");
        ofs << "ship,weight\n";
        for (int s = 0; s < (int)params.cargoWeight.size(); ++s) {
            ofs << s << "," << params.cargoWeight[s] << "\n";
        }
        ofs.close();
    }

    // 4) cargoDensity (s,k,value)
    {
        std::ofstream ofs(baseName + "_cargoDensity.csv");
        ofs << "ship,k,value\n";
        for (int s = 0; s < (int)params.cargoDensity.size(); ++s) {
            for (int k = 0; k < (int)params.cargoDensity[s].size(); ++k) {
                ofs << s << "," << k << "," << params.cargoDensity[s][k] << "\n";
            }
        }
        ofs.close();
    }

    // 5) maxResponseAngle (s,k,value)
    {
        std::ofstream ofs(baseName + "_maxResponseAngle.csv");
        ofs << "ship,k,value\n";
        for (int s = 0; s < (int)params.maxResponseAngle.size(); ++s) {
            for (int k = 0; k < (int)params.maxResponseAngle[s].size(); ++k) {
                ofs << s << "," << k << "," << params.maxResponseAngle[s][k] << "\n";
            }
        }
        ofs.close();
    }

    // 6) requiredSlots (s,k,value)
    {
        std::ofstream ofs(baseName + "_requiredSlots.csv");
        ofs << "ship,k,value\n";
        for (int s = 0; s < (int)params.requiredSlots.size(); ++s) {
            for (int k = 0; k < (int)params.requiredSlots[s].size(); ++k) {
                ofs << s << "," << k << "," << params.requiredSlots[s][k] << "\n";
            }
        }
        ofs.close();
    }

    // 7) unloadingSpeed (s,b,k,value)
    {
        std::ofstream ofs(baseName + "_unloadingSpeed.csv");
        ofs << "ship,berth,k,value\n";
        for (int s = 0; s < (int)params.unloadingSpeed.size(); ++s) {
            for (int b = 0; b < (int)params.unloadingSpeed[s].size(); ++b) {
                for (int k = 0; k < (int)params.unloadingSpeed[s][b].size(); ++k) {
                    ofs << s << "," << b << "," << k << "," << params.unloadingSpeed[s][b][k] << "\n";
                }
            }
        }
        ofs.close();
    }

    // 8) transshipmentCost (b,r,v,value)
    {
        std::ofstream ofs(baseName + "_transshipmentCost.csv");
        ofs << "berth,row,slot,value\n";
        for (int b = 0; b < (int)params.transshipmentCost.size(); ++b) {
            for (int r = 0; r < (int)params.transshipmentCost[b].size(); ++r) {
                for (int v = 0; v < (int)params.transshipmentCost[b][r].size(); ++v) {
                    ofs << b << "," << r << "," << v << "," << params.transshipmentCost[b][r][v] << "\n";
                }
            }
        }
        ofs.close();
    }

    // 9) storageCost (s,k,r,value)
    {
        std::ofstream ofs(baseName + "_storageCost.csv");
        ofs << "ship,k,row,value\n";
        for (int s = 0; s < (int)params.storageCost.size(); ++s) {
            for (int k = 0; k < (int)params.storageCost[s].size(); ++k) {
                for (int r = 0; r < (int)params.storageCost[s][k].size(); ++r) {
                    ofs << s << "," << k << "," << r << "," << params.storageCost[s][k][r] << "\n";
                }
            }
        }
        ofs.close();
    }
}

// 合并导出：每行一个 ship，列展开为各个可对应的参数，方便单表分析
void writeParamsCombinedCSV(const ModelParams& params, const std::string& baseName) {
    // 确保输出目录存在
    {
        std::string dir = dirname_of(baseName);
        if (!dir.empty()) {
            mkdir_p(dir + "/");
        }
    }
    std::string file = baseName + "_combined.csv";
    std::ofstream ofs(file);
    if (!ofs.is_open()) {
        std::cerr << "无法打开文件写入: " << file << std::endl;
        return;
    }

    // 构建表头
    std::vector<std::string> headers;
    headers.push_back("ship");
    headers.push_back("arrivalTime");
    headers.push_back("cargoWeight");

    // requiredSlots per k
    for (int k = 0; k < params.numShipK; ++k) {
        headers.push_back("requiredSlots_k" + std::to_string(k));
    }

    // cargoDensity per k
    for (int k = 0; k < params.numShipK; ++k) {
        headers.push_back("cargoDensity_k" + std::to_string(k));
    }

    // maxResponseAngle per k
    for (int k = 0; k < params.numShipK; ++k) {
        headers.push_back("maxResponseAngle_k" + std::to_string(k));
    }

    // unloadingSpeed: for each berth and k
    for (int b = 0; b < params.numBerths; ++b) {
        for (int k = 0; k < params.numShipK; ++k) {
            headers.push_back("unloadingSpeed_b" + std::to_string(b) + "_k" + std::to_string(k));
        }
    }

    // storageCost: for each row and k
    for (int r = 0; r < params.numRows; ++r) {
        for (int k = 0; k < params.numShipK; ++k) {
            headers.push_back("storageCost_r" + std::to_string(r) + "_k" + std::to_string(k));
        }
    }

    // 写表头
    for (size_t i = 0; i < headers.size(); ++i) {
        if (i) ofs << ',';
        ofs << headers[i];
    }
    ofs << '\n';

    // 写每一行（每艘船）
    for (int s = 0; s < params.numShips; ++s) {
        bool first = true;
        auto put = [&](const std::string &v){ if (!first) ofs<<','; ofs<<v; first=false; };

        put(std::to_string(s));
        put(std::to_string(params.arrivalTime[s]));
        put(std::to_string(params.cargoWeight[s]));

        // requiredSlots
        for (int k = 0; k < params.numShipK; ++k) {
            int val = 0;
            if (s < (int)params.requiredSlots.size() && k < (int)params.requiredSlots[s].size()) val = params.requiredSlots[s][k];
            put(std::to_string(val));
        }

        // cargoDensity
        for (int k = 0; k < params.numShipK; ++k) {
            double v = 0.0;
            if (s < (int)params.cargoDensity.size() && k < (int)params.cargoDensity[s].size()) v = params.cargoDensity[s][k];
            put(std::to_string(v));
        }

        // maxResponseAngle
        for (int k = 0; k < params.numShipK; ++k) {
            double v = 0.0;
            if (s < (int)params.maxResponseAngle.size() && k < (int)params.maxResponseAngle[s].size()) v = params.maxResponseAngle[s][k];
            put(std::to_string(v));
        }

        // unloadingSpeed per berth
        for (int b = 0; b < params.numBerths; ++b) {
            for (int k = 0; k < params.numShipK; ++k) {
                double v = 0.0;
                if (s < (int)params.unloadingSpeed.size() && b < (int)params.unloadingSpeed[s].size() && k < (int)params.unloadingSpeed[s][b].size()) v = params.unloadingSpeed[s][b][k];
                put(std::to_string(v));
            }
        }

        // storageCost per row
        for (int r = 0; r < params.numRows; ++r) {
            for (int k = 0; k < params.numShipK; ++k) {
                double v = 0.0;
                if (s < (int)params.storageCost.size() && k < (int)params.storageCost[s].size() && r < (int)params.storageCost[s][k].size()) v = params.storageCost[s][k][r];
                put(std::to_string(v));
            }
        }

        ofs << '\n';
    }

    ofs.close();
}

// 当作独立可执行使用的入口（合并原 data_init_runner 功能）
int main(int argc, char** argv) {
    int numBerths = 9;
    int numShips = 50;
    if (argc >= 3) {
        try {
            numBerths = std::stoi(argv[1]);
            numShips = std::stoi(argv[2]);
        } catch (...) {
            std::cerr << "参数解析失败，使用默认值 numBerths=2 numShips=3" << std::endl;
        }
    }

    std::cout << "Running data_init (merged) with numBerths=" << numBerths << " numShips=" << numShips << std::endl;
    ModelParams params = setParams(numBerths, numShips);

    // Ensure output directory exists: caller created cpp/data earlier; just write into it.
    writeParamsToCSV(params, "data/example_L12/params_output");
    // writeParamsCombinedCSV(params, "data/example_2/params_output_combined");

    std::cout << "Data export completed to cpp/data/*.csv" << std::endl;
    return 0;
}
