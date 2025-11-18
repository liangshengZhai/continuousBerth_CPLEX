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
			else if (key == "safe_distance") params.safe_distance = std::stoi(val);
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
	params.shipLength.assign(params.numShips, 0); // 初始化船长
	params.storageCost.assign(params.numShips, std::vector<std::vector<double>>(params.numShipK, std::vector<double>(params.numRows, 0.0)));
	params.unloadingSpeed.assign(params.numShips, std::vector<double>(params.numShipK, 0.0));

	// shipLength.csv (ship,value)
	{
		std::ifstream ifs(baseName + "_shipLength.csv");
		if(ifs.is_open()){
			std::string h; std::getline(ifs,h);
			std::string l;
			while(std::getline(ifs,l)){
				if(l.empty()) continue;
				std::istringstream ss(l);
				std::string s,v;
				if(!std::getline(ss,s,',')) continue;
				if(!std::getline(ss,v)) continue;
				try{ int si=std::stoi(s); int iv=std::stoi(v); if(si>=0 && si<params.numShips) params.shipLength[si]=iv; }catch(...){}}
		}
	}

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
		// 健壮性检查：如有任何关键参数为0，直接报错退出，防止空变量被加入模型
		if (params.numShips <= 0 || params.numShipK <= 0 || params.numRows <= 0 || params.numSlotsPerRow <= 0 || params.numCrane <= 0) {
			std::cerr << "[ERROR] 参数维度为0，无法建立模型。" << std::endl;
			return 1;
		}
		// x_skrv: 船舶s货舱k的货物是否分配到行r的槽v
		IloArray<IloArray<IloArray<IloArray<IloBoolVar>>>> x(env);
		// h_skrv: 船舶s的货物是否结束于行r的槽v
		IloArray<IloArray<IloArray<IloArray<IloBoolVar>>>> h(env);
		// f_skr: 船舶s货舱的货物是否分配到行r
		IloArray<IloArray<IloArray<IloBoolVar>>> f(env);
		// y_st: 船舶s的停靠位置是否在船舶t之前
		IloArray<IloArray<IloBoolVar>> y(env);
		// z_skq: 船舶s货舱k由起重机q负责
		IloArray<IloArray<IloArray<IloBoolVar>>> z(env);
		// q_skt 船舱卸货的顺序变量
		IloArray<IloArray<IloArray<IloBoolVar>>> q(env);
		// e_s: 船舶s的卸载开始时间
		IloArray<IloNumVar> e(env);
		//e_sk:船舶s 货舱k的卸货时间
		IloArray<IloArray<IloNumVar>> e_sk(env);
		// a_s: 船舶s停靠的起始坐标
		IloArray<IloNumVar> a_s(env);

		// 初始化变量
		for (int s = 0; s < params.numShips; s++) {
			// z_skq: 船舶s货舱k由起重机q负责
			IloArray<IloArray<IloBoolVar>> z_s(env);
			for (int k = 0; k < params.numShipK; k++) {
				IloArray<IloBoolVar> z_sk(env);
				for (int q_ = 0; q_ < params.numCrane; q_++) {
					std::string z_name = "z_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(q_);
					z_sk.add(IloBoolVar(env, z_name.c_str()));
				}
				z_s.add(z_sk);
			}
			z.add(z_s);

			// q_skt: 船舱卸货顺序变量
			IloArray<IloArray<IloBoolVar>> q_s(env);
			for (int k = 0; k < params.numShipK; k++) {
				IloArray<IloBoolVar> q_sk(env);
				for (int t = 0; t < params.numShipK; t++) {
					std::string q_name = "q_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(t);
					q_sk.add(IloBoolVar(env, q_name.c_str()));
				}
				q_s.add(q_sk);
			}
			q.add(q_s);

			// y_st: 船舶s的停靠位置是否在船舶t之前
			IloArray<IloBoolVar> y_s(env);
			for (int t = 0; t < params.numShips; t++) {
				std::string y_name = "y_" + std::to_string(s) + "_" + std::to_string(t);
				y_s.add(IloBoolVar(env, y_name.c_str()));
			}
			y.add(y_s);

			// f_skr: 船舶s货舱的货物是否分配到行r
			IloArray<IloArray<IloBoolVar>> f_s(env);
			for (int k = 0; k < params.numShipK; k++) {
				IloArray<IloBoolVar> f_sk(env);
				for (int r = 0; r < params.numRows; r++) {
					std::string f_name = "f_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(r);
					f_sk.add(IloBoolVar(env, f_name.c_str()));
				}
				f_s.add(f_sk);
			}
			f.add(f_s);

			// e_s: 船舶s的卸载开始时间
			e.add(IloNumVar(env, params.arrivalTime[s], params.planningHorizon));

			// e_sk: 船舶s 货舱k的卸货时间
			IloArray<IloNumVar> e_s(env);
			for (int k = 0; k < params.numShipK; k++) {
				e_s.add(IloNumVar(env, params.arrivalTime[s], IloInfinity));
			}
			e_sk.add(e_s);

			// x_skrv: 船舶s货舱k的货物是否分配到行r的槽v
			IloArray<IloArray<IloArray<IloBoolVar>>> x_s(env);
			for (int k = 0; k < params.numShipK; k++) {
				IloArray<IloArray<IloBoolVar>> x_sk(env);
				for (int r = 0; r < params.numRows; r++) {
					IloArray<IloBoolVar> x_skr(env);
					for (int v = 0; v < params.numSlotsPerRow; v++) {
						std::string x_name = "x_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(r) + "_" + std::to_string(v);
						x_skr.add(IloBoolVar(env, x_name.c_str()));
					}
					x_sk.add(x_skr);
				}
				x_s.add(x_sk);
			}
			x.add(x_s);

			// h_skrv: 船舶s的货物是否结束于行r的槽v
			IloArray<IloArray<IloArray<IloBoolVar>>> h_s(env);
			for (int k = 0; k < params.numShipK; k++) {
				IloArray<IloArray<IloBoolVar>> h_sk(env);
				for (int r = 0; r < params.numRows; r++) {
					IloArray<IloBoolVar> h_skr(env);
					for (int v = 0; v < params.numSlotsPerRow; v++) {
						std::string h_name = "h_" + std::to_string(s) + "_" + std::to_string(k) + "_" + std::to_string(r) + "_" + std::to_string(v);
						h_skr.add(IloBoolVar(env, h_name.c_str()));
					}
					h_sk.add(h_skr);
				}
				h_s.add(h_sk);
			}
			h.add(h_s);

			// a_s: 船舶s停靠的起始坐标
			a_s.add(IloNumVar(env, 0, IloInfinity));
		}

		// 4. 构建目标函数：最小化总转运成本、存储成本和靠泊时间
		// 若任何参数为0，直接跳过模型构建
		if (params.numShips == 0 || params.numShipK == 0 || params.numRows == 0 || params.numSlotsPerRow == 0 || params.numCrane == 0) {
			std::cerr << "[ERROR] 参数维度为0，跳过模型构建。" << std::endl;
			return 1;
		}
		IloExpr objExpr(env);
		int objTermCount = 0;
		// 总转运成本（去除IloPower的0.5次幂，直接用距离平方dist2，且跳过requiredSlots为0的项）
		for (int s = 0; s < params.numShips; s++) {
			for(int k =0 ;k < params.numShipK;k++){
				if(params.requiredSlots[s][k] == 0 || params.unloadingSpeed[s][k] == 0) continue;
				for (int r = 0; r < params.numRows; r++) {
					for (int v = 0; v < params.numSlotsPerRow; v++) {
						// 槽位中心坐标 (x_rv, y_rv)
						double x_rv = v * params.width + params.width / 2.0;
						double y_rv = r * params.width + params.width / 2.0;
						// 船中心坐标 (a_s + l_s/2, y0)
						IloExpr ship_center = (a_s[s] + params.shipLength[s]) / 2.0;
						double y0 = 0.0; // 你可以根据实际情况设定
						// 曼哈顿距离：横向+纵向绝对值
						objExpr += (IloAbs(ship_center - x_rv) + IloAbs(y0 - y_rv)) * params.cargoWeight[s] / (params.requiredSlots[s][k]*params.numShipK) * x[s][k][r][v];
						objTermCount++;
					}
				}
			}
		}
		std::cout << "[DEBUG] 目标函数有效项数: " << objTermCount << std::endl;
		if(objTermCount == 0) {
			std::cerr << "[ERROR] 目标函数无有效项，模型无法建立。请检查参数数据。" << std::endl;
			return 1;
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
				if (speed <= 0) continue; // 跳过非法速度
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
		// 5. 添加约束条件
        
		// 约束: 船泊停靠的起始位置大于y0，且终止位置小于码头总长度
		double y0 = 0.0; // 泊位起点
		double berthEnd = params.Long; // 码头总长度
		for (int s = 0; s < params.numShips; s++) {
			// 起始位置大于y0
			model.add(a_s[s] >= y0);
			// 终止位置小于码头总长度
			model.add(a_s[s] + params.shipLength[s] <= berthEnd);
		}

        
        // 约束: 船舶间停靠需要保持安全距离，防止重叠
		double M = 1e6; // 足够大的常数
		for (int s = 0; s < params.numShips; s++) {
			for (int t = 0; t < params.numShips; t++) {
				if (s == t) continue;
				// y[s][t] == 1 表示s在t左侧
				model.add(a_s[s] + params.shipLength[s] + params.safe_distance <= a_s[t] + M * (1 - y[s][t]));
				model.add(a_s[t] + params.shipLength[t] + params.safe_distance <= a_s[s] + M * y[s][t]);
			}
		}

        //约束：每艘船的每个货舱只能使用一个起重机
		for (int s = 0; s < params.numShips; s++) {
			for (int k = 0; k < params.numShipK; k++) {
				IloExpr con(env);
				for (int q = 0; q < params.numCrane; q++) {
					con += z[s][k][q];
				}
				model.add(con == 1);
				con.end();
			}
		}

		//约束：一个起重机在同一时刻，只能服务于一艘船的一个船舱
		// 时间冲突Big-M约束：同一时刻同一台起重机不能服务于不同船舱
		double Mbig = 1e6;
		for (int q = 0; q < params.numCrane; q++) {
			for (int s = 0; s < params.numShips; s++) {
				for (int k = 0; k < params.numShipK; k++) {
					for (int t = 0; t < params.numShips; t++) {
						for (int m = 0; m < params.numShipK; m++) {
							if (s == t && k == m) continue;
							// 作业时长 = params.cargoWeight[s] / (params.unloadingSpeed[s][k] * params.numShipK)
							double duration_sk = params.cargoWeight[s] / (params.unloadingSpeed[s][k] * params.numShipK);
							// 非线性乘积z[s][k][q]*z[t][m][q]，Cplex不支持，需线性化。先注释掉防止模型报错。
							// model.add(e_sk[s][k] + duration_sk <= e_sk[t][m] + Mbig * (1 - z[s][k][q] * z[t][m][q]));
						}
					}
				}
			}
		}

		//约束同一艘船内任意两个不同货舱的卸货顺序互斥：q[s][k][t] + q[s][t][k] == 1
		for (int s = 0; s < params.numShips; s++) {
			for (int k = 0; k < params.numShipK; k++) {
				for (int t = 0; t < params.numShipK; t++) {
					if (k == t) continue;
					model.add(q[s][k][t] + q[s][t][k] == 1);
				}
			}
		}

		//约束：前一个货舱作业完成之后，后一个货舱才能开始作业
		for (int s = 0; s < params.numShips; s++) {
			for (int k = 0; k < params.numShipK; k++) {
				for (int t = 0; t < params.numShipK; t++) {
					if (k == t) continue;
					if(params.unloadingSpeed[s][t] > 0 && params.unloadingSpeed[s][k] > 0) {
						double duration_st = params.cargoWeight[s] / (params.unloadingSpeed[s][t] * params.numShipK);
						double duration_sk = params.cargoWeight[s] / (params.unloadingSpeed[s][k] * params.numShipK);
						// e_st + duration_st <= e_sk * q[s][k][t]
						model.add(e_sk[s][t] + duration_st <= e_sk[s][k] + (1 - q[s][k][t]) * 1e6);
						// e_sk + duration_sk <= e_st * q[s][t][k]
						model.add(e_sk[s][k] + duration_sk <= e_sk[s][t] + (1 - q[s][t][k]) * 1e6);
					}
				}
			}
		}
		
		// 约束(3.11): 每艘船的每个舱占用足够的槽数
		for (int s = 0; s < params.numShips; s++) {
			for(int k = 0 ; k< params.numShipK;k++){
				if(params.requiredSlots[s][k] == 0) continue;
				IloExpr con(env);
				for (int r = 0; r < params.numRows; r++) {
					for (int v = 0; v < params.numSlotsPerRow; v++) {
						con += x[s][k][r][v];
					}
				}
				model.add(con == params.requiredSlots[s][k]);
				con.end();
			}
		}

        // 约束(3.12): 每个槽最多放一种货物（跨船舶 s 和货舱 k，总和 <= 1）
        // 之前的实现仅对每个 k 单独约束，允许相同槽被不同 k 的货物占用，造成重复占用的问题。
        for (int r = 0; r < params.numRows; r++) {
            for (int v = 0; v < params.numSlotsPerRow; v++) {
                IloExpr con(env);
                for (int s = 0; s < params.numShips; s++) {
                    for (int k = 0; k < params.numShipK; k++) {
                        con += x[s][k][r][v];
                    }
                }
                model.add(con <= 1);
                con.end();
            }
        }
        
        // 约束(3.13): 每艘船的货物存储在同一行
        for (int s = 0; s < params.numShips; s++) {
            for(int k = 0 ; k < params.numShipK;k++){
                IloExpr con(env);
                for (int r = 0; r < params.numRows; r++) {
                    con += f[s][k][r];
                }
                model.add(con == 1);
                con.end();
            }
        }
        
        // 约束(3.14): x_srv与f_sr的关联
		for (int s = 0; s < params.numShips; s++) {
			for(int k =0 ; k <params.numShipK;k++){
				for (int r = 0; r < params.numRows; r++) {
					if(params.numSlotsPerRow == 0) continue;
					IloExpr con(env);
					for (int v = 0; v < params.numSlotsPerRow; v++) {
						con += x[s][k][r][v];
					}
					model.add(con <= params.numSlotsPerRow * f[s][k][r]);
					con.end();
				}
			}
		}
        
        // 约束(12)-(14): 存储槽的连续性（简化实现，完整逻辑需按文档详细处理）
        for (int s = 0; s < params.numShips; s++) {
            for(int k =0; k< params.numShipK;k++){
                for (int r = 0; r < params.numRows; r++) {
                    IloExpr con(env);
                    for (int v = 0; v < params.numSlotsPerRow; v++) {
                        con += h[s][k][r][v];
                    }
                    model.add(con == f[s][k][r]);
                    con.end();
                    
                    // 约束(13): 最后一个槽的h_srv约束
                    model.add(x[s][k][r][params.numSlotsPerRow-1] <= h[s][k][r][params.numSlotsPerRow-1]);
                    
                    // 约束(14): 中间槽的连续性约束
                    for (int v = 0; v < params.numSlotsPerRow-1; v++) {
                        model.add(x[s][k][r][v] - x[s][k][r][v+1] <= h[s][k][r][v]);
                    }
                }
            }
        }
        
        //约束船舱卸货顺序
        for(int s = 0; s< params.numShips;s++){
            for(int k = 0 ; k< params.numShipK;k++){
                // IloExpr con(env);
                model.add(e[s] <= e_sk[s][k]);
            }
        }
		// 7. 求解模型
        IloCplex cplex(model);
        cout <<"导出模型"<<endl;
        // cout <<"导出模型"<<endl;
        // cplex.setOut(env.getNullStream()); // 关闭输出
        cplex.setParam(IloCplex::TiLim, 3600); // 设置时间限制为1小时
		 // 计时：使用 CPLEX 的计时（与当前 ClockType 一致：CPU/WallClock/Deterministic）
        double t0 = cplex.getCplexTime();
        bool solved = cplex.solve();
        double solveSeconds = cplex.getCplexTime() - t0;
        if(solved){
			std::cout << "求解状态: " << cplex.getStatus() << std::endl;
			std::cout << "目标值: " << cplex.getObjValue() << std::endl;
			std::cout << "求解时间 (秒): " << solveSeconds << std::endl;
			// 输出部分决策变量结果作为示例
			for (int s = 0; s < params.numShips; s++) { // 仅输出前两艘船
				std::cout << "船舶 " << s << " 的卸载开始时间: " << cplex.getValue(e[s]) << std::endl;
				std::cout << "船舶 " << s << " 的停靠起始位置: " << cplex.getValue(a_s[s]) << std::endl;
			}
		} else {
			std::cout << "未找到可行解，状态: " << cplex.getStatus() << std::endl;
		}


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

