#include <vector>
#include <tuple>
#include <iostream>
#include <fstream> 
#include <sstream>

#include "rapidcsv.h"

using DP_TABLE = std::vector<std::vector<long>>;
using PENALTY_MAP = std::map<std::tuple<char,char>, long>;

//Evaluate if diagonal outcome of Needleman-Wunsch
long match_or_mismatch(char c1, char c2)
{
    return (c1 == c2) ? 1L : -1L;
}

std::vector<long> NWScore(const std::string& X, const std::string& Y, PENALTY_MAP& penalty)
{
    const int n = X.length();
    const int m = Y.length();

    DP_TABLE Score;
    Score.resize(2, std::vector<long>(m+1, 0));

    Score[0][0]=0;
    
    for (int j=1; j<=m; j++) {
        Score[0][j] = Score[0][j-1] - penalty[std::make_tuple('-', Y[j-1])];
    }
   
    for (int i=1; i<=n; i++) {
        Score[1][0] = Score[0][0] - penalty[std::make_tuple(X[i-1], '-')];
        for (int j=1; j<=m; j++) {
            Score[1][j] = std::max({
                               Score[1][j-1] - penalty[std::make_tuple('-', Y[j-1])],
                               Score[0][j] - penalty[std::make_tuple(X[i-1], '-')],
                               Score[0][j-1] + match_or_mismatch(X[i-1],Y[j-1])*penalty[std::make_tuple(X[i-1],Y[j-1])]
                            });
        }
        Score[0] = Score[1];
    }
    
    return Score[1];
    
}



std::pair < std::string, std::string > NeedlemanWunsch (const std::string& X, const std::string& Y, PENALTY_MAP& penalty)
{
    int n = X.length(), m = Y.length();

    DP_TABLE Score;

    Score.resize(n+1, std::vector<long>(m+1, 0));
    //STEP 1: assign first row and column
    Score[0][0] = 0;
    for (int i=1;i<n+1;i++) {
        Score[i][0] = Score[i-1][0] - penalty[std::make_tuple(X[i-1], '-')];
    }

    for (int j=1;j<m+1;j++) {
        Score[0][j] = Score[0][j-1] - penalty[std::make_tuple('-', Y[j-1])];
    }
    
    //STEP 2: Needelman-Wunsch
    for (int i=1;i<n+1;i++) {
        for (int j=1;j<m+1;j++) {
            Score[i][j] = std::max({Score[i-1][j-1] + match_or_mismatch(X[i-1],Y[j-1])*penalty[std::make_tuple(X[i-1],Y[j-1])],
                          Score[i][j-1] - penalty[std::make_tuple('-', Y[j-1])],
                          Score[i-1][j] - penalty[std::make_tuple(X[i-1], '-')]});
        }
    }

    
    std::string A_1 = "";
    std::string A_2 = "";
    int i = n, j = m;
    while (i>0 || j>0){
        if (i>0 && j>0 && (Score[i][j] == Score[i-1][j-1] + match_or_mismatch(X[i-1],Y[j-1])*penalty[std::make_tuple(X[i-1],Y[j-1])])){
            A_1 = X[i-1] + A_1;
            A_2 = Y[j-1] + A_2;
            i--;
            j--;
        } else if (i>0 && (Score[i][j] == Score[i-1][j] - penalty[std::make_tuple(X[i-1], '-')])){
            A_1 = X[i-1] + A_1;
            A_2 = '-' + A_2;
            i--;
        } else {
            A_1 = '-' + A_1;
            A_2 = Y[j-1] + A_2;
            j--;
        }
    }
    
    std::pair < std::string, std::string > alignment_pair;
    alignment_pair.first = A_1;
    alignment_pair.second = A_2;
    return alignment_pair;
}


std::tuple< std::string, std::string, PENALTY_MAP> Hirschberg(const std::string& X, const std::string& Y, PENALTY_MAP& penalty){
    const int n = X.length();
    const int m = Y.length();

    std::string Z,W;
    
    if (n==0){
        for (int i=1; i<=m; i++){
            Z += '-';
            W += Y[i-1];
        }
    }
    
    else if (m==0){
        for (int i=1; i<=n; i++){
            Z += X[i-1];
            W += '-';
        }
    }
    
    else if (n==1 || m ==1){
        std::pair< std::string, std::string> result = NeedlemanWunsch(X,Y, penalty);

        Z = result.first;
        W = result.second;
    }
    
    else {
        int xmid = n/2;
        std::string X_to_xmid=X.substr(0,xmid),
        
        X_from_xmid = X.substr(xmid);

        std::string X_from_xmid_rev(X_from_xmid.rbegin(), X_from_xmid.rend());

        std::string Y_rev(Y.rbegin(), Y.rend());
        
        std::vector<long> scoreL = NWScore(X_to_xmid,Y, penalty);
        std::vector<long> scoreR = NWScore(X_from_xmid_rev,Y_rev, penalty);
        
        std::reverse(scoreR.begin(), scoreR.end());

        std::transform(scoreL.begin(), scoreL.end(), scoreR.begin(), scoreL.begin(),std::plus<long>( ));
        
        int ymid = static_cast<int>(std::distance(scoreL.begin(), max_element(scoreL.begin(), scoreL.end())));
        
        std::string Y_to_ymid=Y.substr(0,ymid);

        std::string Y_from_ymid=Y.substr(ymid);
        
        std::tuple<std::string, std::string, PENALTY_MAP> left = Hirschberg(X_to_xmid, Y_to_ymid, penalty);
        
        std::tuple<std::string, std::string, PENALTY_MAP> right = Hirschberg(X_from_xmid, Y_from_ymid, penalty);

        Z = std::get<0>(left) + std::get<0>(right);
        W = std::get<1>(left) + std::get<1>(right);   
    }
    
    return std::make_tuple(Z,W,penalty);
}

int main(int argc, char* argv[]){
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "<input_file> <penalty_function_csv>" << std::endl;

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

    rapidcsv::Document doc(argv[2]);

    PENALTY_MAP penalty;

    std::vector<char> firstchar = doc.GetColumn<char>(0);
    std::vector<char> secondchar = doc.GetColumn<char>(1);
    std::vector<long> penalty_number = doc.GetColumn<long>(2);

    for (int i = 0; i < penalty_number.size(); ++i) {
        penalty[std::make_tuple(firstchar[i], secondchar[i])] = penalty_number[i];
    }
    
    std::tuple<std::string, std::string, PENALTY_MAP> ZWpair = Hirschberg(x,y, penalty);
    std::string Z = std::get<0>(ZWpair);
    std::string W = std::get<1>(ZWpair);

    long global_alignment_value = 0;

    for (int i = 0; i < Z.length(); ++i){
        global_alignment_value += penalty[std::make_tuple(Z[i], W[i])];
    }

    std::cout << global_alignment_value << std::endl;

    return 0;
}


