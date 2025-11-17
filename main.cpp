#include <ilcplex/ilocplex.h>
#include <iostream>

int main() {
	IloEnv env;
	try {
		IloModel model(env);
		IloNumVar x(env, 0.0, IloInfinity, ILOFLOAT, "x");
		IloNumVar y(env, 0.0, IloInfinity, ILOFLOAT, "y");
		// 约束: x + y <= 1
		model.add(x + y <= 1);
		// 目标: maximize x + 2y
		model.add(IloMaximize(env, x + 2 * y));

		IloCplex cplex(model);
		if (!cplex.solve()) {
			std::cout << "CPLEX 求解失败!" << std::endl;
			env.end();
			return 1;
		}
		std::cout << "最优目标值: " << cplex.getObjValue() << std::endl;
		std::cout << "x = " << cplex.getValue(x) << std::endl;
		std::cout << "y = " << cplex.getValue(y) << std::endl;
	} catch (IloException& e) {
		std::cerr << "Concert 异常: " << e << std::endl;
	} catch (...) {
		std::cerr << "未知异常!" << std::endl;
	}
	env.end();
	return 0;
}
