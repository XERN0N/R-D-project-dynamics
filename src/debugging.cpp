#include <iostream>
#include <format>


void print_array(double *array, int m, int n) {
    int number_of_elements = m*n;

    int string_lengths[number_of_elements];
    int max_length = 0;
    for (int i = 0; i < number_of_elements; i++) {
        string_lengths[i] = std::format("{}", array[i]).length();
        if (string_lengths[i] > max_length) {
            max_length = string_lengths[i];
        }
    }
    
    std::cout << '[';
    for (int i = 0; i < m; i++) {
        if (i==0) {
            std::cout << '[';
        } else {
            std::cout << " [";        
        }
        for (int j = 0; j < n; j++) {
            std::string leading_whitespaces;
            if (max_length != string_lengths[i*m+j]) {
                leading_whitespaces.append(std::string(max_length-string_lengths[i*m+j], ' '));
            }
            std::cout << leading_whitespaces << array[i*m+j];
            if (j < n-1) {
                std::cout << ", ";
            }
        }
        if (i < m-1) {
            std::cout << "],\n";
        } else {
            std::cout << "]]\n";
        }
    }
}