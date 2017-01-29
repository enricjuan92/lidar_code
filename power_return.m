% Return power
function p_r = power_return(R, r_pbl, alpha_aer, alpha_mol, beta_aer, beta_mol, K)

if R <= r_pbl
    p_r = K / R^2 * (beta_aer + beta_mol) * exp(-2 * (alpha_aer + alpha_mol) * R);
elseif R > r_pbl
    p_r = K / R^2 * (beta_mol) * exp(-2 * (alpha_aer + alpha_mol) * r_pbl) * exp(-2 * alpha_mol * (R - r_pbl));
end

end

