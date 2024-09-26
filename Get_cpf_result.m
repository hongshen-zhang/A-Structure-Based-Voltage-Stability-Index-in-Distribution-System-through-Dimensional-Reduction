function [lam_truth,vol_truth] = Get_cpf_result(mpc,posi_load,perturb_set)
mpopt = mpoption('OUT_ALL',0,'VERBOSE',0);
Test_num = numel(perturb_set(1,:));
lam_truth = zeros(1,Test_num);
base_step = 1000;
parfor loop = 1 : Test_num
    m = 1;
    mpc_temp = mpc;
    mpc_temp.bus(posi_load,3) = mpc_temp.bus(posi_load,3) + real(perturb_set(:,loop)) * base_step;
    mpc_temp.bus(posi_load,4) = mpc_temp.bus(posi_load,4) + imag(perturb_set(:,loop)) * base_step;
    cpf_temp = runcpf(mpc,mpc_temp,mpopt);
    while cpf_temp.cpf.max_lam > 20
        m = m * cpf_temp.cpf.max_lam;
        mpc_temp = mpc;
        mpc_temp.bus(posi_load,3) = mpc_temp.bus(posi_load,3) + real(perturb_set(:,loop)) ...
            * base_step * m;
        mpc_temp.bus(posi_load,4) = mpc_temp.bus(posi_load,4) + imag(perturb_set(:,loop)) ...
            * base_step * m;
        cpf_temp = runcpf(mpc,mpc_temp,mpopt);
    end
    lam_truth(loop) = cpf_temp.cpf.max_lam * base_step * m;
    vol_truth(:,loop) = cpf_temp.cpf.V(:,end);
end
end