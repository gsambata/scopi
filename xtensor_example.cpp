#include <iostream>
#include <xtensor/xarray.hpp>

int main() {
    // Creating a 2D array (tensor) with integer values
    xt::xarray<int> arr{{1, 2, 3}, {4, 5, 6}};
    
    // Print the array
    std::cout << "Array: \n" << arr << std::endl;
    
    // Element-wise addition of two arrays
    xt::xarray<int> arr2{{7, 8, 9}, {10, 11, 12}};
    auto sum = arr + arr2;
    
    std::cout << "Sum of arrays: \n" << sum << std::endl;
    
    // Slicing the array (selecting a sub-array)
    auto slice = arr({0, 1}, {0, 2});  // Select first two rows and first two columns
    std::cout << "Slice of array: \n" << slice << std::endl;

    // Element-wise multiplication (Hadamard product)
    auto product = arr * arr2;
    std::cout << "Element-wise product: \n" << product << std::endl;

    return 0;
}


