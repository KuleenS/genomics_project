#include <iostream>
#include <fstream> 
#include <string> 
#include <vector>
#include <sstream>
#include <map>
#include <tuple>
#include <Kokkos_Core.hpp>

#include "algos.hpp"
#include "block_based_alignment.hpp"

#include "rapidcsv.h"


using DP_TABLE = Kokkos::View<long**>;
using PENALTY_MAP = std::map<std::tuple<char,char>, long>;

void affine_dp_step(DP_TABLE V, DP_TABLE G, DP_TABLE E, DP_TABLE F, PENALTY_MAP penalty, char x, char y, int i, int j, long gap_open, long gap_extension) {
    // std::cout << "(" << x << ", " << y << ")" << '\n';
    if (x == y) {
        G(i, j) = V(i-1, j) + penalty.at(std::make_tuple(x, y));
    } else {
        G(i, j) = V(i-1, j) - penalty.at(std::make_tuple(x, y));
    }
    
    E(i, j) = std::max(E(i, j-1), V(i, j-1) - gap_open) - gap_extension;
    F(i, j) = std::max(F(i-1, j), V(i-1, j) - gap_open) - gap_extension;

    V(i, j) = std::max({E(i, j), F(i, j), G(i, j)});
}


void affine(std::string x, std::string y, int iterating_mode, PENALTY_MAP& penalty, long gap_open, long gap_extension, int block_size){
    long length_x = x.length();
    long length_y = y.length();

    DP_TABLE V("V", length_x + 1, length_y + 1);
    DP_TABLE G("G", length_x + 1, length_y + 1);
    DP_TABLE E("E", length_x + 1, length_y + 1);
    DP_TABLE F("F", length_x + 1, length_y + 1);

    Kokkos::parallel_for("Initialize_VE", Kokkos::RangePolicy<>(1, length_x + 1), KOKKOS_LAMBDA(int i) {
        V(i, 0) = -gap_open - i * gap_extension;
        E(i, 0) = -gap_open - i * gap_extension;
    });

    Kokkos::parallel_for("Initialize_VF", Kokkos::RangePolicy<>(1, length_y + 1), KOKKOS_LAMBDA(int j) {
        V(0, j) = -gap_open - j * gap_extension;
        F(0, j) = -gap_open - j * gap_extension;
    });

    if (block_size > 0) {
        int internal_iterating_mode = iterating_mode % 10;
        int block_iterating_mode = iterating_mode / 10;

        BlockBasedIterator iter(x.length(), y.length(), block_size, internal_iterating_mode, block_iterating_mode);

        for (; iter.has_more(); ++iter) {
            int i = std::get<0>(iter.get_location());
            int j = std::get<1>(iter.get_location());

            affine_dp_step(V,G,E,F,penalty, x[i+1],y[j+1], i+1,j+1, gap_open, gap_extension);
        }
    } else {
        switch (iterating_mode){
          case LEFT_RIGHT:
              Kokkos::parallel_for("Affine_Left_Right", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {x.length()+1, y.length()+1}), KOKKOS_LAMBDA(int i, int j) {
                  // std::cout << i << ',' << j << '\n';
                  affine_dp_step(V, G, E, F, penalty, x[i-1], y[j-1], i, j, gap_open, gap_extension);
              });
              break;
          case UP_DOWN:
              Kokkos::parallel_for("Affine_Up_Down", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {y.length()+1, x.length()+1}), KOKKOS_LAMBDA(int j, int i) {
                  affine_dp_step(V, G, E, F, penalty, x[i-1], y[j-1], i, j, gap_open, gap_extension);
              });
              break;
          case DIAGONAL:
             for (int k = 0; k < x.length() + y.length() + 1; k++) {
                int start_i = std::max(0, k - (int) y.length() - 1);
                int end_i = std::min((int) x.length() + 1, k);

                Kokkos::parallel_for("diagonal_loop", Kokkos::RangePolicy<>(start_i, end_i + 1), KOKKOS_LAMBDA(int i) {
                    int j = k - i;
                    if (j >= 1 && j <= y.length() && i >=1 && i <= x.length()) {  // Ensure i and j are within bounds for length_x and length_y   
                        affine_dp_step(V, G, E, F, penalty, x[i - 1], y[j - 1], i, j, gap_open, gap_extension);
                    }
                });
                Kokkos::fence();  // Ensure all work in this diagonal is completed before moving to the next
            }
            break;
        }
    }
    // for (int i = 0; i < length_x+1; ++i) {
    //     for (int j = 0; j < length_y+1; ++j) {
    //         std::cout << std::setw(10) << V(i, j) << " ";  // Format with 10 spaces
    //     }
    //     std::cout << "\n";  // Newline after each row
    // }

    std::cout << V(length_x, length_y) << std::endl;
}



void global_local(std::string x, std::string y, int iterating_mode, bool global, PENALTY_MAP penalty, int block_size){
    long int (*dp_step)(DP_TABLE, char, char, int, int, PENALTY_MAP);

    if (!global) {
        dp_step = smith_waterman;
    } else {
        dp_step = needleman_wunsch;
    }

    long length_x = x.length();
    long length_y = y.length();

    DP_TABLE dp("dp", length_x + 1, length_y + 1);

    if (global) {
        dp(0,0) = 0;
        Kokkos::parallel_for("Initialize_Global_2", Kokkos::RangePolicy<>(1, length_x + 1), KOKKOS_LAMBDA(int i) {
            // std::cout << "(" << x[i] << ", " << "-" << ")" << std::endl;
            dp(i, 0) = i * penalty.at(std::make_tuple(x[i-1], '-'));
        });
        Kokkos::parallel_for("Initialize_Global",  Kokkos::RangePolicy<>(1, length_y + 1), KOKKOS_LAMBDA(int i) {
            // std::cout << "(" << '-' << ", " << y[i] << ")" << std::endl;
            dp(0, i) = i * penalty.at(std::make_tuple('-', y[i-1]));
        });
    } 
    else{
        //DO NOTHING (Kokkos views initialize to 0)
    }

    Kokkos::View<long int> best_item("best_item");
    Kokkos::deep_copy(best_item, -1.0);

    if (block_size > 0) {
        int internal_iterating_mode = iterating_mode % 10;
        int block_iterating_mode = iterating_mode / 10;

        BlockBasedIterator iter(x.length(), y.length(), block_size, internal_iterating_mode, block_iterating_mode);

        for (; iter.has_more(); ++iter) {
            
            int i = std::get<0>(iter.get_location());
            int j = std::get<1>(iter.get_location());

            dp(i+1,j+1) = dp_step(dp, x[i], y[j], i+1, j+1, penalty);

            if (!global) {
                Kokkos::atomic_max(&best_item(), dp(i+1, j+1));
            }
        }
    } else {

        switch (iterating_mode){
            case LEFT_RIGHT:
            Kokkos::parallel_for("Global_Local_Left_Right", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {x.length()+1, y.length()+1}), KOKKOS_LAMBDA(int i, int j) {
                dp(i, j) = dp_step(dp, x[i-1], y[j-1], i, j, penalty);
                if (!global) {
                    Kokkos::atomic_max(&best_item(), dp(i, j));
                }
            });
            break;
        case UP_DOWN:
            Kokkos::parallel_for("Global_Local_Up_Down", Kokkos::MDRangePolicy<Kokkos::Rank<2>>({1, 1}, {y.length()+1, x.length()+1}), KOKKOS_LAMBDA(int j, int i) {
                dp(i, j) = dp_step(dp, x[i-1], y[j-1], i, j, penalty);
                if (!global) {
                    Kokkos::atomic_max(&best_item(), dp(i, j));
                }
            });
            break;
        case DIAGONAL:
            for (int k = 0; k < x.length() + y.length() + 1; k++) {
                int start_i = std::max(0, k - (int) y.length() - 1);
                int end_i = std::min((int) x.length() + 1, k);

                Kokkos::parallel_for("diagonal_loop", Kokkos::RangePolicy<>(start_i, end_i + 1), KOKKOS_LAMBDA(int i) {
                    int j = k - i;
                    if (j >= 1 && j <= length_y && i >=1 && i <= length_x) {  // Ensure i and j are within bounds for length_x and length_y
                        dp(i, j) = dp_step(dp, x[i-1], y[j-1], i, j, penalty);  // Replace with your function
                        // Update best_item atomically if not in global mode
                        if (!global) {
                            Kokkos::atomic_max(&best_item(), dp(i, j));
                        }
                    }
                });
                Kokkos::fence();  // Ensure all work in this diagonal is completed before moving to the next
            }
            break;
        }

    }
    // for (int i = 0; i < length_x+1; ++i) {
    //     for (int j = 0; j < length_y+1; ++j) {
    //         std::cout << std::setw(10) << dp(i, j) << " ";  // Format with 10 spaces
    //     }
    //     std::cout << "\n";  // Newline after each row
    // }
    if (!global) {
        long int host_max_val;
        Kokkos::deep_copy(host_max_val, best_item);
        std::cout << host_max_val << std::endl;
    } else {
        std::cout << dp(length_x, length_y) << std::endl;
    }
}


int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv); // Initialize Kokkos

    if (argc != 5 && argc != 7) {
        std::cerr << "Usage: " << argv[0] << "<input_file> affine <iterating_mode> <penalty_function_csv> <gap_open> <gap_extension>" << std::endl;
        std::cerr << "Usage: " << argv[0] << "<input_file> <global|local> <iterating_mode> <penalty_function_csv>" << std::endl;
        Kokkos::finalize(); // Finalize Kokkos
        return 1;
    }

    std::ifstream input(argv[1]);
    if (!input.is_open()) {
        std::cerr << "File " << argv[1] << " does not exist" << std::endl;
        Kokkos::finalize(); // Finalize Kokkos
        return 1;
    }

    std::string x, y;
    std::getline(input, x);
    std::getline(input, y);
    input.close();

    std::istringstream global_local_affine_ss(argv[2]);
    std::string global_local_affine;
    if (!(global_local_affine_ss >> global_local_affine)) {
        std::cerr << "Invalid string: " << argv[2] << '\n';
        Kokkos::finalize(); // Finalize Kokkos
        return 1;
    } else if (!global_local_affine_ss.eof()) {
        std::cerr << "Trailing characters: " << argv[2] << '\n';
        Kokkos::finalize(); // Finalize Kokkos
        return 1;
    }

    if (global_local_affine != "global" && global_local_affine != "local" && global_local_affine != "affine") {
        std::cerr << "Impossible mode: " << argv[2] << '\n';
        Kokkos::finalize(); // Finalize Kokkos
        return 1;
    }

    std::istringstream iterating_mode_ss(argv[3]);
    int iterating_mode;
    if (!(iterating_mode_ss >> iterating_mode)) {
        std::cerr << "Invalid number: " << argv[3] << '\n';
        Kokkos::finalize(); // Finalize Kokkos
        return 1;
    } else if (!iterating_mode_ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv[3] << '\n';
        Kokkos::finalize(); // Finalize Kokkos
        return 1;
    }



    int offset = 0;

    int internal_iterating_mode = 0;
    int block_iterating_mode = 0;
    int block_size = 0;

    if (iterating_mode == BLOCK_BASED) {

        std::istringstream internal_iterating_mode_ss(argv[4]);

        if (!(internal_iterating_mode_ss >> internal_iterating_mode)) {
            std::cerr << "Invalid number: " << argv[4] << '\n';
            return 1;
        } else if (!internal_iterating_mode_ss.eof()) {
            std::cerr << "Trailing characters after number: " << argv[4] << '\n';
            return 1;
        } else if (internal_iterating_mode == BLOCK_BASED) {
            std::cerr << "Nested block-based iteration is not allowed! \n";
            return 1;
        }  

        std::istringstream block_iterating_mode_ss(argv[5]);

        if (!(block_iterating_mode_ss >> block_iterating_mode)) {
            std::cerr << "Invalid number: " << argv[5] << '\n';
            return 1;
        } else if (!block_iterating_mode_ss.eof()) {
            std::cerr << "Trailing characters after number: " << argv[5] << '\n';
            return 1;
        } else if (block_iterating_mode == BLOCK_BASED) {
            std::cerr << "Nested block-based iteration is not allowed! \n";
            return 1;
        }

        std::istringstream  block_size_ss(argv[6]);

        if (!(block_size_ss >> block_size)) {
            std::cerr << "Invalid number: " << argv[6] << '\n';
            return 1;
        } else if (!block_size_ss.eof()) {
            std::cerr << "Trailing characters after number: " << argv[6] << '\n';
            return 1;
        } else if (block_size <= 0) {
            std::cerr << "Block size must be positive, got: " << argv[6] << '\n';
            return 1;
        } 

        iterating_mode = block_iterating_mode * 10 + internal_iterating_mode;
        offset = 3;

    }

    rapidcsv::Document doc(argv[4+offset]);


    PENALTY_MAP penalty;
    std::vector<char> firstchar = doc.GetColumn<char>(0);
    std::vector<char> secondchar = doc.GetColumn<char>(1);
    std::vector<long> penalty_number = doc.GetColumn<long>(2);

    for (int i = 0; i < penalty_number.size(); ++i) {
        // std::cout << firstchar[i] << ", " << secondchar[i] << '\n';
        penalty[std::make_tuple(firstchar[i], secondchar[i])] = penalty_number[i];
    } 



    if (global_local_affine == "affine"){
        std::istringstream gap_open_ss(argv[5+offset]);


        long gap_open;
        if (!(gap_open_ss >> gap_open)) {

            std::cerr << "Invalid number: " << argv[5+offset] << '\n';
            return 1;
        } else if (!gap_open_ss.eof()) {
            std::cerr << "Trailing characters after number: " << argv[5+offset] << '\n';
            return 1;
        }

        std::istringstream gap_extension_ss(argv[6+offset]);

        long gap_extension;
        if (!(gap_extension_ss >> gap_extension)) {
            std::cerr << "Invalid number: " << argv[6+offset] << '\n';
            return 1;
        } else if (!gap_extension_ss.eof()) {
            std::cerr << "Trailing characters after number: " << argv[6+offset] << '\n';
            return 1;
        }

        affine(x,y,iterating_mode, penalty,gap_open, gap_extension, block_size);

    } else{
        global_local(x,y,iterating_mode,global_local_affine == "global", penalty, block_size);
    }

    Kokkos::finalize(); // Finalize Kokkos
    return 0;
}