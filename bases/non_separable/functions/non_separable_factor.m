function factor = non_separable_factor(n, a, b, i, j, domain_length)
    factor = domain_length ^ (a + b) * factorial(n) * factorial(i + a) * factorial(j + b) / (factorial(i) * factorial(j) * factorial(n + a + b));
end