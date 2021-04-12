#include <iostream>
#include <cmath>
#include <limits>
#include <boost/multiprecision/cpp_bin_float.hpp>

template <typename T>
T f_1(T x)
{
	return x * x - 1;
}

template <typename T>
T derivative(T (*f)(T), const T &x)
{
	const T dx = sqrt(std::numeric_limits<T>::epsilon()) * x;
	return ((*f)(x + dx) - (*f)(x)) / dx;
}

template <typename T>
T newtons_iteration(T (*f)(T), const T &x)
{
	return x - f(x) / derivative(f, x);
}

int main()
{

	std::cout.precision(100);
	std::cout << std::fixed;

	std::cout << "Task #2: Solving x^2 == 1 using Newton's iteration \nUsing initial approximation x = 1.1...\n\n";

	boost::multiprecision::cpp_bin_float_100 x("1.1"), x_prime = 0;

	for (int i = 0; x != x_prime; ++i)
	{
		std::cout << "Iteration #" << i << ": x = " << x;

		x_prime = x;
		x = newtons_iteration(f_1, x);

		std::cout << " (press ENTER to continue)";
		std::cin.get();
	}

	return 0;
}