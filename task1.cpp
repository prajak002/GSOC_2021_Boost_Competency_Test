#include <iostream>
#include <cmath>
#include <limits>
#include <boost/multiprecision/cpp_bin_float.hpp>

int main()
{

	std::cout.precision(100);
	std::cout << std::fixed;

	std::cout << "\nTask #1: Gauss-Lengendre algorithm for Pi\n\n";

	using mp_float = boost::multiprecision::cpp_bin_float_100;

	mp_float a(1), b = 1 / sqrt(mp_float(2)), t("0.25"), pi = 0;
	int correct_digits = 0;
	for (int i = 1; correct_digits < 99; ++i)
	{
		std::cout << "Iteration #" << i << ": ";

		mp_float a_prime = (a + b) / 2;
		mp_float b_prime = sqrt(a * b);
		mp_float t_prime = t - exp2(i - 1) * (a - a_prime) * (a - a_prime);

		a = a_prime;
		b = b_prime;
		t = t_prime;

		pi = (a + b) * (a + b) / (4 * t);

		correct_digits = (int)-log10(abs(pi - boost::math::constants::pi<mp_float>()));

		std::cout << "pi = " << pi;
		std::cout << "\n | " << correct_digits << " digits correct | \n";
		std::cout << " (press ENTER to continue)";
		std::cin.get();
	}

	return 0;
}