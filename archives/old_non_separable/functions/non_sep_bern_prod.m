function non_sep_bern_prod_out = non_sep_bern_prod(n, i, j, a_1, a_2)
a = 1;
for k = 1 : a_1
    a = a * (i + k);
end
b = 1;
for k = 1 : a_2
    b = b * (j + k);
end
c = 1;
for k = 1 : a_1 + a_2
    c = c * (n + k);
end
non_sep_bern_prod_out = a * b / c;
end