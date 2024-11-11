#include <vector>
#include <tuple>
#include <map>
#include <Kokkos_Core.hpp>

#include <bits/stdc++.h>

long smith_waterman(Kokkos::View<long**> dp, char x, char y, int i, int j, std::map<std::tuple<char,char>, long> penalty){
    return std::max(
        {
            dp(i-1,j-1) + penalty.at(std::make_tuple(x, y)),
            dp(i-1,j) + penalty.at(std::make_tuple(x, '-')),
            dp(i,j-1) + penalty.at(std::make_tuple('-', y)),
            0L
        }
    );
}

long needleman_wunsch(Kokkos::View<long**> dp, char x, char y, int i,int j, std::map<std::tuple<char,char>, long> penalty){
    return std::min(
        {
            dp(i-1,j-1) + penalty.at(std::make_tuple(x, y)),
            dp(i-1,j) + penalty.at(std::make_tuple(x, '-')),
            dp(i,j-1) + penalty.at(std::make_tuple('-', y)),
        }
    );
}