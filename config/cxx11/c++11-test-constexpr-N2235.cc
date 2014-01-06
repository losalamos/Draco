#include <limits>

constexpr int square(int x)
{
    return x*x;
}

constexpr int the_answer()
{
    return 42;
}

template <typename T>
constexpr T maxTvalue()
{
    return std::numeric_limits<T>::max();
}

int main()
{
    constexpr float factor = 1.0/(1.0 + maxTvalue<int>());
    int test_arr[square(3)];
    bool ret = (
        (square(the_answer()) == 1764) &&
        (sizeof(test_arr)/sizeof(test_arr[0]) == 9) &&
        (factor < 1.0)
        );
    return ret ? 0 : 1;
}
