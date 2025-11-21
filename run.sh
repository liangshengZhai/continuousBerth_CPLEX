g++ -std=c++11 \
-I/Applications/CPLEX_Studio2211/cplex/include \
-I/Applications/CPLEX_Studio2211/concert/include \
-L/Applications/CPLEX_Studio2211/cplex/lib/arm64_osx/static_pic \
-L/Applications/CPLEX_Studio2211/concert/lib/arm64_osx/static_pic \
main.cpp -o main \
-lilocplex -lcplex -lconcert -lm -lpthread -ldl