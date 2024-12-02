#include "block_based_alignment.hpp"

#include "algos.hpp"
#include <cassert>




BlockBasedIterator::BlockBasedIterator(int length_x, int length_y, int box_size, int micro_strategy, int macro_strategy)
    : m_length_x(length_x)
    , m_length_y(length_y)
    , m_box_size(box_size)
    , m_macro_strategy(macro_strategy)
    , m_micro_strategy(micro_strategy)
    , micro_x_pos(0)
    , micro_y_pos(0)
    , micro_k(0)
    , box_x_pos(0)
    , box_y_pos(0)
    , box_k(0)
    , can_iterate(true)
{
    assert(m_length_x % m_box_size == 0);
    assert(m_length_y % m_box_size == 0);

    num_box_x = m_length_x / m_box_size;
    num_box_y = m_length_y / m_box_size;
}

void BlockBasedIterator::macro_step() {
    
    micro_x_pos = 0;
    micro_y_pos = 0;
    micro_k = 0;

    switch (m_macro_strategy){
        case LEFT_RIGHT:
            box_x_pos ++;

            if (box_y_pos >= num_box_y) {
                box_y_pos = 0;
                box_x_pos ++;
            }

            if (box_x_pos >= num_box_x) {
                // we have reached the end of our iteration
                can_iterate = false;
            };
            return;
        
        case UP_DOWN:

            micro_x_pos ++;

            if (box_x_pos >= num_box_x) {
                box_x_pos = 0;
                box_y_pos ++;
            }

            if (box_y_pos >= num_box_y) {
                can_iterate = false;
            }

            return;

        case DIAGONAL:

            do {

                box_y_pos ++;

                if (box_y_pos > box_k) {
                    box_y_pos = 0;
                    box_k ++;
                }

                if (box_k > num_box_x + num_box_y - 2) {
                    can_iterate = false;
                    return;
                }

            box_x_pos = box_k - box_y_pos;

            } while (box_x_pos < num_box_x && box_y_pos < num_box_y);

            return;
    }
}

BlockBasedIterator &BlockBasedIterator::operator++() {

    assert(can_iterate);

    switch (m_micro_strategy){
        case LEFT_RIGHT:
            micro_y_pos ++;

            if (micro_y_pos >= m_box_size) {
                micro_y_pos = 0;
                micro_x_pos ++;
            }

            if (micro_x_pos >= m_box_size) {
                macro_step();
            };
            return *this;
        
        case UP_DOWN:

            micro_x_pos ++;

            if (micro_x_pos >= m_box_size) {
                micro_x_pos = 0;
                micro_y_pos ++;
            }

            if (micro_y_pos >= m_box_size) {
                macro_step();
            }

            return *this;

        case DIAGONAL:

            do {

                micro_y_pos ++;

                if (micro_y_pos > micro_k) {
                    micro_y_pos = 0;
                    micro_k ++;
                }

                if (micro_k > 2*m_box_size - 2) {
                    macro_step();
                    return *this;
                }

            micro_x_pos = micro_k - micro_y_pos;

            } while (micro_x_pos < m_box_size && micro_y_pos < m_box_size);


            return *this;

    }
    return *this;
}

bool BlockBasedIterator::has_more() {
    return can_iterate;
}

std::tuple<int, int> BlockBasedIterator::get_location() {
    int true_x = box_x_pos * m_box_size + micro_x_pos;
    int true_y = box_y_pos * m_box_size + micro_y_pos;

    return std::tuple<int, int>(true_x, true_y);
}