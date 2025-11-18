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


int main() {
	//初始化CPLEX环境和模型
    IloEnv env;
    IloModel model(env);
	try{
		//读取模型参数
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
		
		// 3. 定义决策变量
		// x_skrv: 船舶s货舱k的货物是否分配到行r的槽v
		IloArray<IloArray<IloArray<IloArray<IloBoolVar>>>> x(env, params.numShips);
		// h_skrv: 船舶s的货物是否结束于行r的槽v
		IloArray<IloArray<IloArray<IloArray<IloBoolVar>>>> h(env, params.numShips);
		// f_skr: 船舶s货舱的货物是否分配到行r
		IloArray<IloArray<IloArray<IloBoolVar>>> f(env, params.numShips);
		// y_st: 船舶s的停靠位置是否在船舶t之前
		IloArray<IloArray<IloBoolVar>> y(env, params.numShips);
		// z_skq: 船舶s货舱k由起重机q负责
		IloArray<IloArray<IloArray<IloBoolVar>>> z(env, params.numShips);
		// q_skt 船舱卸货的顺序变量
		IloArray<IloArray<IloArray<IloBoolVar>>> q(env, params.numShips);
		// e_s: 船舶s的卸载开始时间
		IloArray<IloNumVar> e(env, params.numShips);
		//e_sk:船舶s 货舱k的卸货时间
		IloArray<IloArray<IloNumVar>> e_sk(env, params.numShips);
		// a_s: 船舶s停靠的起始坐标
		IloArray<IloNumVar> a_s(env, params.numShips);

		// 初始化变量
        for (int s = 0; s < params.numShips; s++) {
			// z_skq: 船舶s货舱k由起重机q负责
			IloArray<IloArray<IloBoolVar>> z_s(env, params.numShipK);
			for (int k = 0; k < params.numShipK; k++) {
				z_s[k] = IloArray<IloBoolVar>(env, params.numCrane);
				for (int q = 0; q < params.numCrane; q++) {
					std::string z_name = "z_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(q);
					z_s[k][q] = IloBoolVar(env, z_name.c_str());
				}
			}
			z[s] = z_s;

			// q_skt: 船舱卸货顺序变量
			IloArray<IloArray<IloBoolVar>> q_s(env, params.numShipK);
			for (int k = 0; k < params.numShipK; k++) {
				q_s[k] = IloArray<IloBoolVar>(env, params.numShipK);
				for (int t = 0; t < params.numShipK; t++) {
					std::string q_name = "q_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(t);
					q_s[k][t] = IloBoolVar(env, q_name.c_str());
				}
			}
			q[s] = q_s;

			// y_st: 船舶s的停靠位置是否在船舶t之前
			IloArray<IloBoolVar> y_s(env, params.numShips);
			for (int t = 0; t < params.numShips; t++) {
				std::string y_name = "y_" + std::to_string(s) + "_" + std::to_string(t);
				y_s[t] = IloBoolVar(env, y_name.c_str());
			}
			y[s] = y_s;

			// f_skr: 船舶s货舱的货物是否分配到行r
			IloArray<IloArray<IloBoolVar>> f_s(env, params.numShipK);
			for (int k = 0; k < params.numShipK; k++) {
				f_s[k] = IloArray<IloBoolVar>(env, params.numRows);
				for (int r = 0; r < params.numRows; r++) {
					std::string f_name = "f_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(r);
					f_s[k][r] = IloBoolVar(env, f_name.c_str());
				}
			}
			f[s] = f_s;

			// e_s: 船舶s的卸载开始时间
			e[s] = IloNumVar(env, params.arrivalTime[s], params.planningHorizon);

			// e_sk: 船舶s 货舱k的卸货时间
			IloArray<IloNumVar> e_s(env, params.numShipK);
			for (int k = 0; k < params.numShipK; k++) {
				e_s[k] = IloNumVar(env, params.arrivalTime[s], IloInfinity);
			}
			e_sk[s] = e_s;

			// x_skrv: 船舶s货舱k的货物是否分配到行r的槽v
			IloArray<IloArray<IloArray<IloBoolVar>>> x_s(env, params.numShipK);
			for (int k = 0; k < params.numShipK; k++) {
				x_s[k] = IloArray<IloArray<IloBoolVar>>(env, params.numRows);
				for (int r = 0; r < params.numRows; r++) {
					x_s[k][r] = IloArray<IloBoolVar>(env, params.numSlotsPerRow);
					for (int v = 0; v < params.numSlotsPerRow; v++) {
						std::string x_name = "x_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(r) + "_" + std::to_string(v);
						x_s[k][r][v] = IloBoolVar(env, x_name.c_str());
					}
				}
			}
			x[s] = x_s;

			// h_skrv: 船舶s的货物是否结束于行r的槽v
			IloArray<IloArray<IloArray<IloBoolVar>>> h_s(env, params.numShipK);
			for (int k = 0; k < params.numShipK; k++) {
				h_s[k] = IloArray<IloArray<IloBoolVar>>(env, params.numRows);
				for (int r = 0; r < params.numRows; r++) {
					h_s[k][r] = IloArray<IloBoolVar>(env, params.numSlotsPerRow);
					for (int v = 0; v < params.numSlotsPerRow; v++) {
						std::string h_name = "h_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(r) + "_" + std::to_string(v);
						h_s[k][r][v] = IloBoolVar(env, h_name.c_str());
					}
				}
			}
			h[s] = h_s;
			// a_s: 船舶s停靠的起始坐标
			a_s[s] = IloNumVar(env, 0, IloInfinity);
		}

		// 4. 构建目标函数：最小化总转运成本、存储成本和靠泊时间
        IloExpr objExpr(env);
        // 总转运成本
        for (int s = 0; s < params.numShips; s++) {
			for(int k =0 ;k < params.numShipK;k++){
				for (int r = 0; r < params.numRows; r++) {
					for (int v = 0; v < params.numSlotsPerRow; v++) {
						// 槽位中心坐标 (x_rv, y_rv)
						double x_rv = v * params.width + params.width / 2.0;
						double y_rv = r * params.width + params.width / 2.0;
						// 船中心坐标 (a_s + l_s/2, y0)
						IloExpr ship_center = (a_s[s] + params.shipLength[s]) / 2.0;
						double y0 = 0.0; // 你可以根据实际情况设定
						IloExpr dist2 = (ship_center - x_rv) * (ship_center - x_rv) + (y0 - y_rv) * (y0 - y_rv);
						IloExpr dist = IloPower(dist2, 0.5);
						objExpr += dist * params.cargoWeight[s] / (params.requiredSlots[s][k]*params.numShipK) * x[s][k][r][v];
					}
				}
			}
        }
        
        // 总存储成本
        for (int s = 0; s < params.numShips; s++) {
            for(int k =0 ; k <params.numShipK;k++){
                for (int r = 0; r < params.numRows; r++) {
                    for (int v = 0; v < params.numSlotsPerRow; v++) {
                        objExpr += params.storageCost[s][k][r] * x[s][k][r][v];
                    }
                }
            }
        }


       IloExpr berthTime(env);
        // 正确计算每艘船的靠泊时间（考虑泊位分配）
        for (int s = 0; s < params.numShips; s++) {
            IloExpr singleBerth(env);
            singleBerth += e[s] - params.arrivalTime[s];
    
            for(int k =0 ; k< params.numShipK;k++){
                // 仅累加分配泊位的卸载时间
				double speed = params.unloadingSpeed[s][k];
				if (speed <= 0) speed = 1.0; // 防除零
				singleBerth += (params.cargoWeight[s] / (speed * params.numShipK));
                
            }
            berthTime += singleBerth;
            singleBerth.end();
        }
        // 应用权重
        objExpr = params.alpha * objExpr + params.beta * berthTime; // 注意：目标函数公式需根据文档调整权重应用方式
        
        model.add(IloMinimize(env, objExpr));
        berthTime.end();
        objExpr.end();




	}catch (IloException& e) {
		std::cerr << "CPLEX异常: " << e.getMessage() << std::endl;
	} catch (std::exception& e) {
		std::cerr << "标准异常: " << e.what() << std::endl;
	} catch (...) {
		std::cerr << "未知异常发生" << std::endl;
	}
	env.end();
	return 0;
        

}

     // Closing the lambda function properly

