#include <vector>
#include <tuple>
#include <map>

#include <bits/stdc++.h>

long smith_waterman(std::vector<std::vector<long>>& dp, char x, char y, int i, int j, std::map<std::tuple<char,char>, long>& penalty){
    return std::max(
        {
            dp[i-1][j-1] + penalty[std::make_tuple(x, y)],
            dp[i-1][j] + penalty[std::make_tuple(x, '-')],
            dp[i][j-1] + penalty[std::make_tuple('-', y)],
            0L
        }
    );
}

long needleman_wunsch(std::vector<std::vector<long>>& dp, char x, char y, int i,int j, std::map<std::tuple<char,char>, long>& penalty){
    return std::min(
        {
            dp[i-1][j-1] + penalty[std::make_tuple(x, y)],
            dp[i-1][j] + penalty[std::make_tuple(x, '-')],
            dp[i][j-1] + penalty[std::make_tuple('-', y)],
        }
    );
}