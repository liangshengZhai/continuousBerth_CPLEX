#ifndef MODELPARAM_H
#define MODELPARAM_H


#include <vector>
#include <map>
#include <string>
using namespace std;
// 模型参数结构体
struct ModelParams {
    int numRows;           // 存储行数量
    int numSlotsPerRow;    // 每行槽数量
    int numShips;          // 船舶数量
    double planningHorizon; // 计划 horizon
    double Length;          //岸线长度
    int y0;              //岸线起点
    //船的船舱个数
    int numShipK;
    int safe_distance=10; //船舶间最小安全距离
    
    vector<int> shipLength; //船舶长度
    
    // 模型参数
    vector<vector<double>> unloadingSpeed;  // 卸载速度 p_sk
    //具体转运成本由分配了泊位后，通过距离得到
    //vector<vector<vector<double>>> transshipmentCost; // 转运成本 d_brv
    vector<vector<vector<double>>> storageCost;       // 存储成本 phi_skrv. //
    vector<double> arrivalTime;        // 到达时间 tau_s
    vector<double> cargoWeight;        // 货物重量 gamma_s
    vector<vector<double>> cargoDensity;       // 货物密度 rho_s_k
    vector<vector<double>> maxResponseAngle;   // 最大响应角 delta_s_k
    
    // 计算得到的参数
    vector<vector<int>> requiredSlots;         // 每艘船需要的槽数 n_sk
    double width;                      // 存储行宽度 w
    double relativeHeight;             // 相对高度 c
    
    // 权重系数
    double alpha;                      // 转运成本权重
    double beta;                       // 靠泊时间权重
    double gamma;              // 存储成本权重

};

#endif