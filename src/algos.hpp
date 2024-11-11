#include <vector>
#include <map>
#include <Kokkos_Core.hpp>

long smith_waterman(Kokkos::View<long**> dp, char x, char y, int i, int j, std::map<std::tuple<char,char>, long> penalty);
long needleman_wunsch(Kokkos::View<long**> dp, char x, char y, int i,int j, std::map<std::tuple<char,char>, long> penalty);