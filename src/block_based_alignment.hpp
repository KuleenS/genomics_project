#include <iterator>
#include <tuple>

class BlockBasedIterator {
    private:
        int m_length_x;
        int m_length_y;
        int m_micro_strategy;
        int m_macro_strategy;
        int m_box_size;

        int num_box_x;
        int num_box_y;

        int box_x_pos;
        int box_y_pos;
        int box_k;
        
        int micro_x_pos;
        int micro_y_pos;
        int micro_k;

        bool can_iterate;

        void macro_step();

    public:
        BlockBasedIterator(int length_x, int length_y, int box_size, int micro_strategy, int macro_strategy);
        BlockBasedIterator &operator++();

        bool has_more();
        std::tuple<int, int> get_location();
};