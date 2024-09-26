function [critical_value,lam_value,lam_posi,V_step] = Predict_lambda_dis_sys(mpc,gen,load,posi_load,perturb_set)

%% Initialization
set_num = numel(perturb_set(1,:));
nonzero_load = intersect(find(mpc.bus(:,4) ~= 0),load);

Y = makeYbus(mpc);
pf_vol = get_pf_vol(mpc);
gen_vol = pf_vol(gen);
[weight,Y_eq]= Father_matrix(Y,gen,nonzero_load);

nonzero_S = -(mpc.bus(nonzero_load,3) + mpc.bus(nonzero_load,4)*j) / ...
    mpc.baseMVA;


[~,posi_loc] = ismember(posi_load,nonzero_load);
add_S = zeros(numel(nonzero_load),set_num);
add_S(posi_loc,:) = -perturb_set;

%% All operations are based on load
Y_load_load_verse = full(Y(load,load))^(-1);
[~,nonzero_loc] = ismember(nonzero_load,load);
Y_nonzero_nonzero_verse = Y_load_load_verse(nonzero_loc,nonzero_loc);
diag_Y_nonzero_nonzero_verse = diag(Y_nonzero_nonzero_verse);

%% Form the V_critical_mat based on one node perturbation 
V_g = (weight * pf_vol(gen)) ./ Y_eq;


V_critical_mat = V_g - (Y_nonzero_nonzero_verse) .* ...
    (V_g.' ./ diag_Y_nonzero_nonzero_verse.') / 2;
V_critical_load = diag(V_critical_mat);

% lambda * S0 + S1
%% Choose the most vulurable node

S0_set = (conj(Y_nonzero_nonzero_verse) ./...
    (V_critical_mat.') * ...
    add_S) ...
    .*  V_critical_load ./ ...
    conj(diag_Y_nonzero_nonzero_verse);


S1_set = (conj(Y_nonzero_nonzero_verse) ./...
    (V_critical_mat.') * ...
    nonzero_S) ...
    .*  V_critical_load ./ ...
    conj(diag_Y_nonzero_nonzero_verse) * ...
    ones(1,set_num);

[lam_set] = predict_initial(weight,Y_eq,gen_vol,S0_set,S1_set);
[lam_value,lam_posi] = min(abs(lam_set(posi_loc,:)),[],1);


%% Step analysis for the vulurable node
V_step = V_critical_mat(:,posi_loc(lam_posi));
I_step_1 = -conj(nonzero_S) ./ conj(V_step);
I_step_load_1 = zeros(numel(load),set_num);
I_step_load_1(nonzero_loc,:) =  I_step_1;

I_step_2 = -conj(add_S) ./ conj(V_step);
I_step_load_2 = zeros(numel(load),set_num);
I_step_load_2(nonzero_loc,:) =  I_step_2;

%% Get I
I_step = I_step_load_1 + I_step_load_2 .* lam_value;
V_step = Y_load_load_verse(nonzero_loc,:) ...
    * (-I_step - full(Y(load,gen)) * pf_vol(gen));

S0_set = (conj(Y_nonzero_nonzero_verse)  * ...
    (add_S ./ (V_step))) ...
    .*  V_step ./ ...
    conj(diag_Y_nonzero_nonzero_verse);


S1_set = (conj(Y_nonzero_nonzero_verse)  * ...
    (nonzero_S ./ (V_step))) ...
    .*  V_step ./ ...
    conj(diag_Y_nonzero_nonzero_verse);

[lam_set] = predict_initial(weight,Y_eq,gen_vol,S0_set,S1_set);
critical_value = abs(diag(lam_set(posi_loc(lam_posi),:))');



end



function [pf_vol] = get_pf_vol(mpc)
%mpopt = mpoption('OUT_ALL',0,'VERBOSE',0);
%pf = runpf(mpc,mpopt);
%if ~pf.success
%    error('This mpc can not runpf directly.');
%end
pf = runpf(mpc);
rad  = pf.bus(:,9) * pi  / 180;
pf_vol = pf.bus(:,8) .* cos(rad) + pf.bus(:,8) .* sin(rad) * j;
end

function [weight,Y_eq]= Father_matrix(Y,g_index,p_index)
N = size(Y,1);
load_index = setdiff(1:N,g_index);
[~,p_loc] = ismember(p_index,load_index);
B = (full(Y(load_index,load_index)))^(-1);
C = full(Y(load_index,g_index));
weight = -B(p_loc,:) * C ./ (diag(B(p_loc,p_loc)));
Y_eq = diag(eye(numel(p_loc)) ./ diag(B(p_loc,p_loc)));
end


function [lam] = predict_initial(weight,Y_eq,gen_vol,S,S_initial)
A = -4 * imag(conj(S) ./ Y_eq).^2;
B = -8 * imag(conj(S)  ./ Y_eq) .* imag(conj(S_initial) ./ Y_eq) + ...
    4 * real(conj(S) ./ Y_eq) .* abs(weight * gen_vol ./ Y_eq) .^ 2;
C = -4 * imag(conj(S_initial) ./ Y_eq) .^2 + 4 * real(conj(S_initial) ./ Y_eq) ...
    .* abs(weight * gen_vol ./ Y_eq) .^ 2 + abs(weight * gen_vol ./ Y_eq) .^ 4;
lam = (-B - sqrt(B.^2 - 4 * A .* C)) ./ (2 .* A);
lam(find(abs(real(S)) < 1e-9)) = -C(find(abs(real(S)) < 1e-9))  ./ B(find(abs(real(S)) < 1e-9));
end