#include <iostream>
#include <fstream> 
#include <string> 
#include <vector>
#include <sstream>
#include <map>
#include <tuple>

#include "algos.hpp"

#include "rapidcsv.h"

enum Traversal {
  LEFT_RIGHT = 0,
  UP_DOWN = 1,
  DIAGONAL = 2,
};

using DP_TABLE = std::vector<std::vector<long>>;
using PENALTY_MAP = std::map<std::tuple<char,char>, long>;

void affine_dp_step(DP_TABLE& V, DP_TABLE& G, DP_TABLE& E, DP_TABLE& F, PENALTY_MAP& penalty, char x, char y, int i, int j, long gap_open, long gap_extension){

    if (x == y){
        G[i][j] = V[i-1][j-1] + penalty[std::make_tuple(x,y)];
    } else {
        G[i][j] = V[i-1][j-1] - penalty[std::make_tuple(x,y)];
    }
    
    E[i][j] = std::max({E[i][j-1], V[i][j-1] - gap_open}) - gap_extension;
    F[i][j] = std::max({F[i-1][j], V[i-1][j] - gap_open}) - gap_extension;

    V[i][j] = std::max({E[i][j], F[i][j],G[i][j]});

}


void affine(std::string x, std::string y, int iterating_mode, PENALTY_MAP& penalty, long gap_open, long gap_extension){
    long length_x = x.length();

    long length_y = y.length();

    DP_TABLE V;
    V.resize(length_x+1, std::vector<long>(length_y+1, 0));

    DP_TABLE G;
    G.resize(length_x+1, std::vector<long>(length_y+1, 0));

    DP_TABLE E;
    E.resize(length_x+1, std::vector<long>(length_y+1, 0));

    DP_TABLE F;
    F.resize(length_x+1, std::vector<long>(length_y+1, 0));

    for (int i = 0; i < x.length() + 1; ++i){
        V[i][0] = -gap_open - i * gap_extension;
        E[i][0] = -gap_open - i * gap_extension;
    }

    for (int i = 0; i < y.length() + 1; ++i){
        V[0][i] = -gap_open - i * gap_extension;
        F[0][i] = -gap_open - i * gap_extension;
    }

    switch (iterating_mode){
        case LEFT_RIGHT:

            for (int i = 1; i < x.length()+1; i++){
                for (int j = 1; j < y.length()+1; j++){
                    affine_dp_step(V,G,E,F,penalty, x[i],y[j], i,j, gap_open, gap_extension);
                }
            }

            break;
        case UP_DOWN:
            for (int j = 1; j < y.length()+1; j++){
                for (int i = 1; i < x.length()+1; i++){
                    affine_dp_step(V,G,E,F,penalty, x[i],y[j], i,j, gap_open, gap_extension);
                }
            }

            break;
        case DIAGONAL:
            for( int k = 0 ; k <= y.length()+1 + x.length()+1 - 2; k++ ) {
                for( int j = 0 ; j <= k ; j++ ) {
                    int i = k - j;
                    if( i < x.length()+1 && j < y.length()+1  && i != 0 && j != 0) {
                        affine_dp_step(V,G,E,F,penalty, x[i],y[j], i,j, gap_open, gap_extension);
                    }
                }
            }
            break;
    }

    std::cout << V[x.length()][y.length()] << std::endl;
}


void global_local(std::string x, std::string y, int iterating_mode, bool global, PENALTY_MAP& penalty){
    long (*dp_step)(DP_TABLE&, char, char, int, int, PENALTY_MAP&);

    if (!global){ 
        dp_step = smith_waterman;
    } else {
        dp_step = needleman_wunsch;
    }

    long length_a = x.length();

    long length_b = y.length();

    DP_TABLE dp;
    dp.resize(length_a+1, std::vector<long>(length_b+1, 0));

    if (global){
        for (int i = 0; i < y.length() + 1; ++i){
            dp[0][i] = i*penalty[std::make_tuple('-', y[i])];
        }
        for (int i = 0; i < x.length() + 1; ++i){
            dp[i][0] = i*penalty[std::make_tuple('-', x[i])];
        }
    }
    else{
        for (int i = 0; i < y.length() + 1; ++i){
            dp[0][i] = 0;
        }
        for (int i = 0; i < x.length() + 1; ++i){
            dp[i][0] = 0;
        }
    }

    long best_item = 0;

    switch (iterating_mode){
        case LEFT_RIGHT:

            for (int i = 1; i < x.length()+1; i++){
                for (int j = 1; j < y.length()+1; j++){
                    dp[i][j] = dp_step(dp, x[i-1],y[j-1], i,j, penalty);

                    if (!global){
                        best_item = std::max(best_item, dp[i][j]);
                    }
                }
            }

            break;
        case UP_DOWN:
            for (int j = 1; j < y.length()+1; j++){
                for (int i = 1; i < x.length()+1; i++){
                    dp[i][j] = dp_step(dp, x[i-1],y[j-1], i,j, penalty);

                    if (!global){
                        best_item = std::max(best_item, dp[i][j]);
                    }
                }
            }

            break;
        case DIAGONAL:
            for( int k = 0 ; k <= y.length()+1 + x.length()+1 - 2; k++ ) {
                for( int j = 0 ; j <= k ; j++ ) {
                    int i = k - j;
                    if( i < x.length()+1 && j < y.length()+1  && i != 0 && j != 0) {

                        dp[i][j] = dp_step(dp, x[i-1],y[j-1], i,j, penalty);

                        if (!global){
                            best_item = std::max(best_item, dp[i][j]);
                        }
                    }
                }
            }
            break;
    }

    if (!global){
        std::cout << best_item << std::endl;
    } else{
        std::cout << dp[x.length()][y.length()] << std::endl;
    }

}


int main(int argc, char* argv[]){
    if (argc != 5 && argc != 7) {
        std::cerr << "Usage: " << argv[0] << "<input_file> affine <iterating_mode> <penalty_function_csv> <gap_open> <gap_extension>" << std::endl;
        std::cerr << "Usage: " << argv[0] << "<input_file> <global|local> <iterating_mode> <penalty_function_csv>" << std::endl;

        return 1;
    }


    std::ifstream input(argv[1]);

    if (!input.is_open()){
        std::cerr << "File " << argv[1] <<  " does not exist" << std::endl;
        return 1;
    }

    std::string x,y; 

    std::getline(input, x);
    std::getline(input, y);

    input.close();

    std::istringstream global_local_affine_ss(argv[2]);

    std::string global_local_affine;

    if (!(global_local_affine_ss >> global_local_affine)) {
        std::cerr << "Invalid string: " << argv[2] << '\n';
        return 1;
    } else if (!global_local_affine_ss.eof()) {
        std::cerr << "Trailing characters: " << argv[2] << '\n';
        return 1;
    }

    if (global_local_affine != "global" && global_local_affine != "local" && global_local_affine != "affine"){
        std::cerr << "Impossible mode: " << argv[2] << '\n';
        return 1;
    }

    std::istringstream iterating_mode_ss(argv[3]);

    int iterating_mode;

    if (!(iterating_mode_ss >> iterating_mode)) {
        std::cerr << "Invalid number: " << argv[3] << '\n';
        return 1;
    } else if (!iterating_mode_ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv[3] << '\n';
        return 1;
    }

    rapidcsv::Document doc(argv[4]);

    PENALTY_MAP penalty;

    std::vector<char> firstchar = doc.GetColumn<char>(0);
    std::vector<char> secondchar = doc.GetColumn<char>(1);
    std::vector<long> penalty_number = doc.GetColumn<long>(2);

    for (int i = 0; i < penalty_number.size(); ++i){
        penalty[std::make_tuple(firstchar[i], secondchar[i])] = penalty_number[i];
    }

    if (global_local_affine == "affine"){
        std::istringstream gap_open_ss(argv[5]);

        long gap_open;

        if (!(gap_open_ss >> gap_open)) {
            std::cerr << "Invalid number: " << argv[5] << '\n';
            return 1;
        } else if (!gap_open_ss.eof()) {
            std::cerr << "Trailing characters after number: " << argv[5] << '\n';
            return 1;
        }

        std::istringstream gap_extension_ss(argv[6]);

        long gap_extension;

        if (!(gap_extension_ss >> gap_extension)) {
            std::cerr << "Invalid number: " << argv[6] << '\n';
            return 1;
        } else if (!gap_extension_ss.eof()) {
            std::cerr << "Trailing characters after number: " << argv[6] << '\n';
            return 1;
        }

        affine(x,y,iterating_mode, penalty,gap_open, gap_extension);

    } else{
        global_local(x,y,iterating_mode,global_local_affine == "global", penalty);
    }

    return 0;
}
    
