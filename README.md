<h2>Task 1</h2>
Implemented Gauss-Lengendre algorithm for Pi. Using correct_digits (integer) and mp_float (boost::multiprecision::cpp_bin_float_100).
<br>The implementation with float literal can have problems of rounding error.
<br>I have iterated several times to get correct digits up to 100th precision.
<br><br>
The implementation using cpp_bin_float_100 is written in ‘<a href="https://github.com/NewCyberGypsy/GSOC_2021_Boost_Competency_Test/blob/master/task1.cpp">task1.cpp</a>’.

<h2>Task 2</h2>
Implemented Newton's iteration using initial approximation x_prime (integer) and x("1.1") (boost::multiprecision::cpp_bin_float_100).
<br><br>
Made three templates,<br>
1. T f_1(T x) for x * x - 1, <br>
2. T derivative(T (*f)(T), const T &x) for finding derivative, <br>
3. T newtons_iteration(T (*f)(T), const T &x)  for finding newton’s iteration. <br><br>
The implementation having accuracy upto 100th precision is written in ‘<a href="https://github.com/NewCyberGypsy/GSOC_2021_Boost_Competency_Test/blob/master/task2.cpp">task2.cpp</a>’.

<h2>Task 3</h2>
Provided rough idea of quad double implementation. This is not the final product. A lot of improvements are required.
<br><br>
The rough idea is written in ‘<a href="https://github.com/NewCyberGypsy/GSOC_2021_Boost_Competency_Test/blob/master/task3.hpp">task3.hpp</a>’.
