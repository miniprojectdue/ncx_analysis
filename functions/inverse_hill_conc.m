function C = inverse_hill_conc(beta, IC50, h)
%INVERSE_HILL_CONC Convert fractional block beta to concentration.
    beta = min(max(beta, 1e-9), 1 - 1e-9);
    C = IC50 * (beta / (1 - beta))^(1 / h);
end