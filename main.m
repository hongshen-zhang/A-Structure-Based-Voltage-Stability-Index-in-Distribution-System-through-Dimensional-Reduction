mpc = case533mt_hi;
rand_num = 100;

[g_index,p_index,posi_p_index] = Get_positive_load(mpc);

perturb_set = [mpc.bus(posi_p_index,3) + ...
    mpc.bus(posi_p_index,4)*j] * ones(1,rand_num);
mu = 1;
sigma =  sqrt(0.3);
mm = mu + sigma * randn(numel(posi_p_index),rand_num);

perturb_set = perturb_set .*  abs(mm);


[lam_truth,vol_truth] = Get_cpf_result...
    (mpc,posi_p_index,perturb_set);


critical_value = Predict_lambda_dis_sys(mpc,g_index,p_index,posi_p_index,...
    perturb_set / mpc.baseMVA);


error = average_median_error(lam_truth,critical_value)