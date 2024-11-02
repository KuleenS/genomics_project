#include <vector>
#include <map>

long smith_waterman(std::vector<std::vector<long>>& dp, char x, char y, int i, int j, std::map<std::tuple<char,char>, long>& penalty);
long needleman_wunsch(std::vector<std::vector<long>>& dp, char x, char y, int i,int j, std::map<std::tuple<char,char>, long>& penalty);