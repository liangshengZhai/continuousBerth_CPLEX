#include <ilcplex/ilocplex.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <cmath>
#include <random>
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <cstring>
// #include "utils/modelParam.h"

#include "utlis/modelParam.h"
// 集中管理数据输入和输出路径，只需改这里即可
static const std::string INPUT_BASE = "data/example_1/params_output"; // 不带扩展名的前缀
static const std::string OUTPUT_DIR = "output/output_1";               // 输出目录

// 递归创建目录（等价于 mkdir -p）
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


//读取数据


int main() {
	IloEnv env;

		// 1. 读取模型参数（从 data/ 下的 CSV 文件）
    auto loadParamsFromCSV = [](const std::string &baseName) -> ModelParams {
        ModelParams params;
        // try to read general
        std::ifstream ifs_gen(baseName + "_general.csv");
        if (!ifs_gen.is_open()) {
            std::cerr << "无法打开 general 参数文件: " << baseName + "_general.csv" << std::endl;
            // Return empty/default params (caller should handle)
            return params;
        }
        std::string line;
        // skip header
        std::getline(ifs_gen, line);
        while (std::getline(ifs_gen, line)) {
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string key, val;
            if (!std::getline(ss, key, ',')) continue;
            if (!std::getline(ss, val)) continue;
            try {
                if (key == "numCrane") params.numCrane = std::stoi(val);
                // if (key == "crane") params.numcCrane = std::stoi(val);
				else if (key == "numRows") params.numRows = std::stoi(val);
                else if (key == "numSlotsPerRow") params.numSlotsPerRow = std::stoi(val);
                else if (key == "numShips") params.numShips = std::stoi(val);
                else if (key == "planningHorizon") params.planningHorizon = std::stod(val);
                else if (key == "numShipK") params.numShipK = std::stoi(val);
                else if (key == "width") params.width = std::stod(val);
                else if (key == "relativeHeight") params.relativeHeight = std::stod(val);
                else if (key == "alpha") params.alpha = std::stod(val);
                else if (key == "beta") params.beta = std::stod(val);
            } catch (...) {
                // ignore parse errors per-value
            }
        }
        ifs_gen.close();

        // allocate containers based on general
        if (params.numShips <= 0) params.numShips = 0;
        if (params.numRows <= 0) params.numRows = 0;
        if (params.numSlotsPerRow <= 0) params.numSlotsPerRow = 0;
        if (params.numShipK <= 0) params.numShipK = 0;

        params.arrivalTime.assign(params.numShips, 0.0);
        params.cargoWeight.assign(params.numShips, 0.0);
        params.cargoDensity.assign(params.numShips, std::vector<double>(params.numShipK, 0.0));
        params.maxResponseAngle.assign(params.numShips, std::vector<double>(params.numShipK, 0.0));
        params.requiredSlots.assign(params.numShips, std::vector<int>(params.numShipK, 0));
        // params.unloadingSpeed.assign(params.numShips, std::vector<std::vector<double>>(params.numBerths, std::vector<double>(params.numShipK, 0.0)));
        params.storageCost.assign(params.numShips, std::vector<std::vector<double>>(params.numShipK, std::vector<double>(params.numRows, 0.0)));
        params.unloadingSpeed.assign(params.numShips, std::vector<double>(params.numShipK, 0.0));
        // params.transshipmentCost.assign(params.numBerths, std::vector<std::vector<double>>(params.numRows, std::vector<double>(params.numSlotsPerRow, 0.0)));

        // helper to parse CSV lines
        auto parse_two = [&](const std::string &filename, std::function<void(int,double)> fn){
            std::ifstream ifs(filename);
            if(!ifs.is_open()) return;
            std::string h; std::getline(ifs,h);
            std::string l;
            while(std::getline(ifs,l)){
                if(l.empty()) continue;
                std::istringstream ss(l);
                std::string a,b;
                if(!std::getline(ss,a,',')) continue;
                if(!std::getline(ss,b)) continue;
                try{ fn(std::stoi(a), std::stod(b)); }catch(...){}
            }
        };

        // arrival
        parse_two(baseName + "_arrival.csv", [&](int s, double v){ if(s>=0 && s < params.numShips) params.arrivalTime[s]=v; });
        // cargoWeight
        parse_two(baseName + "_cargoWeight.csv", [&](int s, double v){ if(s>=0 && s < params.numShips) params.cargoWeight[s]=v; });

        // cargoDensity (ship,k,value)
        {
            std::ifstream ifs(baseName + "_cargoDensity.csv");
            if(ifs.is_open()){
                std::string h; std::getline(ifs,h);
                std::string l;
                while(std::getline(ifs,l)){
                    if(l.empty()) continue;
                    std::istringstream ss(l);
                    std::string s,k,v;
                    if(!std::getline(ss,s,',')) continue;
                    if(!std::getline(ss,k,',')) continue;
                    if(!std::getline(ss,v)) continue;
                    try{ int si=std::stoi(s), ki=std::stoi(k); double dv=std::stod(v); if(si>=0 && si<params.numShips && ki>=0 && ki<params.numShipK) params.cargoDensity[si][ki]=dv; }catch(...){}
                }
            }
        }

        // maxResponseAngle
        {
            std::ifstream ifs(baseName + "_maxResponseAngle.csv");
            if(ifs.is_open()){
                std::string h; std::getline(ifs,h);
                std::string l;
                while(std::getline(ifs,l)){
                    if(l.empty()) continue;
                    std::istringstream ss(l);
                    std::string s,k,v;
                    if(!std::getline(ss,s,',')) continue;
                    if(!std::getline(ss,k,',')) continue;
                    if(!std::getline(ss,v)) continue;
                    try{ int si=std::stoi(s), ki=std::stoi(k); double dv=std::stod(v); if(si>=0 && si<params.numShips && ki>=0 && ki<params.numShipK) params.maxResponseAngle[si][ki]=dv; }catch(...){}
                }
            }
        }

        // requiredSlots
        {
            std::ifstream ifs(baseName + "_requiredSlots.csv");
            if(ifs.is_open()){
                std::string h; std::getline(ifs,h);
                std::string l;
                while(std::getline(ifs,l)){
                    if(l.empty()) continue;
                    std::istringstream ss(l);
                    std::string s,k,v;
                    if(!std::getline(ss,s,',')) continue;
                    if(!std::getline(ss,k,',')) continue;
                    if(!std::getline(ss,v)) continue;
                    try{ int si=std::stoi(s), ki=std::stoi(k); int iv=std::stoi(v); if(si>=0 && si<params.numShips && ki>=0 && ki<params.numShipK) params.requiredSlots[si][ki]=iv; }catch(...){}
                }
            }
        }

        // unloadingSpeed (ship,k,value)
        {
            std::ifstream ifs(baseName + "_unloadingSpeed.csv");
            if(ifs.is_open()){
                std::string h; std::getline(ifs,h);
                std::string l;
                while(std::getline(ifs,l)){
                    if(l.empty()) continue;
                    std::istringstream ss(l);
                    std::string s,k,v;
                    if(!std::getline(ss,s,',')) continue;
                    if(!std::getline(ss,k,',')) continue;
                    if(!std::getline(ss,v)) continue;
                    try{ int si=std::stoi(s), ki=std::stoi(k); double dv=std::stod(v); if(si>=0 && si<params.numShips && ki>=0 && ki<params.numShipK) params.unloadingSpeed[si][ki]=dv; }catch(...){}
                }
            }
        }

        // storageCost (ship,k,row,value)
        {
            std::ifstream ifs(baseName + "_storageCost.csv");
            if(ifs.is_open()){
                std::string h; std::getline(ifs,h);
                std::string l;
                while(std::getline(ifs,l)){
                    if(l.empty()) continue;
                    std::istringstream ss(l);
                    std::string s,k,r,v;
                    if(!std::getline(ss,s,',')) continue;
                    if(!std::getline(ss,k,',')) continue;
                    if(!std::getline(ss,r,',')) continue;
                    if(!std::getline(ss,v)) continue;
                    try{ int si=std::stoi(s), ki=std::stoi(k), ri=std::stoi(r); double dv=std::stod(v); if(si>=0 && si<params.numShips && ki>=0 && ki<params.numShipK && ri>=0 && ri<params.numRows) params.storageCost[si][ki][ri]=dv; }catch(...){}
                }
            }
        }

        return params;
    };
	ModelParams params = loadParamsFromCSV(INPUT_BASE);
    // Diagnostic print: verify that params were loaded correctly
    std::cout << "[DEBUG] Loaded params: numCrane=" << params.numCrane
              << " numRows=" << params.numRows
              << " numSlotsPerRow=" << params.numSlotsPerRow
              << " numShips=" << params.numShips
              << " numShipK=" << params.numShipK
              << " planningHorizon=" << params.planningHorizon << std::endl;
    std::cout << "[DEBUG] arrivalTime.size=" << params.arrivalTime.size()
              << " cargoWeight.size=" << params.cargoWeight.size()
              << " unloadingSpeed.size=" << params.unloadingSpeed.size()
              << " storageCost.size=" << params.storageCost.size() << std::endl;
	

}

     // Closing the lambda function properly

