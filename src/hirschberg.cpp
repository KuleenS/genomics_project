#include <vector>
#include <tuple>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <Kokkos_Core.hpp>

#include "rapidcsv.h"

using DP_TABLE = Kokkos::View<long**>;
using PENALTY_MAP = std::map<std::tuple<char,char>, long>;

//Evaluate if diagonal outcome of Needleman-Wunsch
long match_or_mismatch(char c1, char c2)
{
    return (c1 == c2) ? 1L : -1L;
}

int NWScore(const std::string& X, const std::string& Y, PENALTY_MAP& penalty) {
    const int n = X.length();
    const int m = Y.length();

    // Initialize Kokkos Views for the DP table
    DP_TABLE Score("Score", 2, m + 1);

    // Set the initial values for the first row
    Kokkos::parallel_for("InitializeFirstRow", m + 1, KOKKOS_LAMBDA(int j) {
        if (j == 0) {
            Score(0, j) = 0;
        } else {
            // std::cout << "(" << '-' << ", " << Y[j - 1] << ")" << std::endl;
            Score(0, j) = Score(0, j - 1) - penalty.at(std::make_tuple('-', Y[j - 1]));
        }
    });
    Kokkos::fence();

    // Fill in the DP table
    for (int i = 1; i <= n; ++i) {
        Score(1, 0) = Score(0, 0) - penalty.at(std::make_tuple(X[i - 1], '-'));

        Kokkos::parallel_for("FillRow", m, KOKKOS_LAMBDA(int j) {
            int col = j + 1;
            Score(1, col) = std::max({
                Score(1, col - 1) - penalty.at(std::make_tuple('-', Y[col - 1])),
                Score(0, col) - penalty.at(std::make_tuple(X[i - 1], '-')),
                Score(0, col - 1) + match_or_mismatch(X[i - 1], Y[col - 1]) * penalty.at(std::make_tuple(X[i - 1], Y[col - 1]))
            });
        });
        Kokkos::fence();

        // Swap rows for the next iteration
        Kokkos::parallel_for("SwapRows", m + 1, KOKKOS_LAMBDA(int j) {
            Score(0, j) = Score(1, j);
        });
        Kokkos::fence();
    }

    // Copy the final row into a host vector
    Kokkos::View<long*> result("Result", m + 1);

    // Copy the desired subview into the result view
    Kokkos::deep_copy(result, Kokkos::subview(Score, 1, Kokkos::ALL));

    // // Print out the result view
    // Kokkos::parallel_for("PrintResult", m + 1, KOKKOS_LAMBDA(int i) {
    //     printf("result[%d] = %ld\n", i, result(i));
    // });

    // Return the bottom-right value of the 2D view (assuming the dimensions of Score are (n, m))
    int x = Score(Score.extent(0) - 1, Score.extent(1) - 1);
    return x;
}



std::pair<std::string, std::string> NeedlemanWunsch(const std::string& X, const std::string& Y, PENALTY_MAP& penalty) {
    const int n = X.length();
    const int m = Y.length();

    // Initialize the Kokkos View for the DP table with dimensions (n+1) x (m+1)
    DP_TABLE Score("Score", n + 1, m + 1);

    // Initialize first row and first column in parallel
    Kokkos::parallel_for("InitFirstColumn", n + 1, KOKKOS_LAMBDA(int i) {
        Score(i, 0) = (i == 0) ? 0 : Score(i - 1, 0) - penalty.at(std::make_tuple(X[i - 1], '-'));
    });

    Kokkos::parallel_for("InitFirstRow", m + 1, KOKKOS_LAMBDA(int j) {
        Score(0, j) = (j == 0) ? 0 : Score(0, j - 1) - penalty.at(std::make_tuple('-', Y[j - 1]));
    });
    Kokkos::fence();

    // Fill in the rest of the DP table in parallel row by row
    for (int i = 1; i <= n; ++i) {
        Kokkos::parallel_for("FillRow", m, KOKKOS_LAMBDA(int j) {
            int col = j + 1;
            Score(i, col) = std::max({
                Score(i - 1, col - 1) + match_or_mismatch(X[i - 1], Y[col - 1]) * penalty.at(std::make_tuple(X[i - 1], Y[col - 1])),
                Score(i, col - 1) - penalty.at(std::make_tuple('-', Y[col - 1])),
                Score(i - 1, col) - penalty.at(std::make_tuple(X[i - 1], '-'))
            });
        });
        Kokkos::fence();
    }

    // Backtracking to determine the optimal alignment
    std::string A_1, A_2;
    int i = n, j = m;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 &&
            (Score(i, j) == Score(i - 1, j - 1) + match_or_mismatch(X[i - 1], Y[j - 1]) * penalty.at(std::make_tuple(X[i - 1], Y[j - 1])))) {
            A_1 = X[i - 1] + A_1;
            A_2 = Y[j - 1] + A_2;
            --i;
            --j;
        } else if (i > 0 && (Score(i, j) == Score(i - 1, j) - penalty.at(std::make_tuple(X[i - 1], '-')))) {
            A_1 = X[i - 1] + A_1;
            A_2 = '-' + A_2;
            --i;
        } else {
            A_1 = '-' + A_1;
            A_2 = Y[j - 1] + A_2;
            --j;
        }
    }

    return {A_1, A_2};
}


// Hirschberg Function using Kokkos
std::pair<std::string, std::string> Hirschberg(const std::string& X, const std::string& Y, PENALTY_MAP& penalty) {
    const int n = X.length();
    const int m = Y.length();

    std::string Z, W;

    if (n == 0) {
        for (int i = 0; i < m; ++i) {
            Z += '-';
            W += Y[i];
        }
    } else if (m == 0) {
        for (int i = 0; i < n; ++i) {
            Z += X[i];
            W += '-';
        }
    } else if (n == 1 || m == 1) {
        auto result = NeedlemanWunsch(X, Y, penalty);
        Z = result.first;
        W = result.second;
    } else {
        int xmid = n / 2;
        std::string X_to_xmid = X.substr(0, xmid);
        std::string X_from_xmid = X.substr(xmid);
        std::string X_from_xmid_rev(X_from_xmid.rbegin(), X_from_xmid.rend());
        std::string Y_rev(Y.rbegin(), Y.rend());

        // Left and right scores
        auto scoreL = NWScore(X_to_xmid, Y, penalty);
        auto scoreR = NWScore(X_from_xmid_rev, Y_rev, penalty);

        Kokkos::View<long*> scoreR_host("scoreR_host", m + 1);
        Kokkos::deep_copy(scoreR_host, scoreR);
        std::reverse(scoreR_host.data(), scoreR_host.data() + m + 1);

        Kokkos::View<long*> scoreL_host("scoreL_host", m + 1);
        Kokkos::deep_copy(scoreL_host, scoreL);

        std::vector<long> totalScores(m + 1);
        for (int j = 0; j <= m; ++j) {
            totalScores[j] = scoreL_host(j) + scoreR_host(j);
        }

        int ymid = std::distance(totalScores.begin(), std::max_element(totalScores.begin(), totalScores.end()));

        std::string Y_to_ymid = Y.substr(0, ymid);
        std::string Y_from_ymid = Y.substr(ymid);

        auto left = Hirschberg(X_to_xmid, Y_to_ymid, penalty);
        auto right = Hirschberg(X_from_xmid, Y_from_ymid, penalty);

        Z = left.first + right.first;
        W = left.second + right.second;
    }

    return {Z, W};
}

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <penalty_function_csv>" << std::endl;
        Kokkos::finalize();
        return 1;
    }

    std::ifstream input(argv[1]);
    if (!input.is_open()) {
        std::cerr << "File " << argv[1] << " does not exist" << std::endl;
        Kokkos::finalize();
        return 1;
    }

    std::string x, y;
    std::getline(input, x);
    std::getline(input, y);
    input.close();

    rapidcsv::Document doc(argv[2]);
    PENALTY_MAP penalty;
    std::vector<char> firstchar = doc.GetColumn<char>(0);
    std::vector<char> secondchar = doc.GetColumn<char>(1);
    std::vector<long> penalty_number = doc.GetColumn<long>(2);

    for (size_t i = 0; i < penalty_number.size(); ++i) {
        penalty[std::make_tuple(firstchar[i], secondchar[i])] = penalty_number[i];
    }

    // for (auto it = penalty.begin(); it != penalty.end(); ++it) {
    //     std::cout << "(" << std::get<0>(it->first) << ", " << std::get<1>(it->first) << ") => " << it->second << std::endl;
    // }


    auto result = NWScore(x, y, penalty);

    // Sum penalties based on the result
    long global_alignment_value = 0;
    
    // for (int i = 0; i < x.length(); ++i) {
    //     std::cout << "(" << x[i] << ", " << y[i] << ")" << std::endl;
    //     global_alignment_value += penalty.at(std::make_tuple(x[i], y[i]));
    // }

    global_alignment_value = result;

    std::cout << global_alignment_value << std::endl;

    Kokkos::finalize();
    return 0;
}