#ifndef DATA_INIT_H
#define DATA_INIT_H

#include "modelParam.h"

//test
// 将原 verify.cpp 中的数据初始化提取到单独文件中
ModelParams setParams(int numBerths, int numShips);

// 将 ModelParams 写出为 Excel 可读的 CSV 文件。
// baseName 会作为文件名前缀生成若干 CSV（如 baseName_arrival.csv）
void writeParamsToCSV(const ModelParams& params, const std::string& baseName = "params_output");

// 生成一个合并的 CSV（每行对应一艘船，列为各参数），便于在单一表格中查看对应关系
// 输出文件: <baseName>_combined.csv
void writeParamsCombinedCSV(const ModelParams& params, const std::string& baseName = "params_output");

// (已移除：XLSX 输出策略) 若需再次启用，请集成 libxlsxwriter 并在此处添加声明。

#endif // DATA_INIT_H
