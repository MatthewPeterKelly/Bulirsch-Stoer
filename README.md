# Bulirsch-Stoer Method 
The Bulirsch-Stoer method is an algorithm that solves initial value problems for smooth functions.

## Method Details:
The algorithm implemented here is largely taken from a [report by Sujit Kirpekar](http://web.mit.edu/ehliu/Public/Spring2006/18.304/implementation_bulirsch_stoer.pdf
), and is also covered in Numerical Recipes in C, section 16.4.

The basic idea is to compute a sequence of approximations to the solution using the modified mid-point method. Then Richardson extrapolation is used to estimate the limit of that sequence as the step size goes to zero.



