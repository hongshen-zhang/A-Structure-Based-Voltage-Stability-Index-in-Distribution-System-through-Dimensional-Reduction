function [gen,load,posi_load] = Get_positive_load(mpc)
g_bus = find(mpc.bus(:,2) ~=1);
gen = intersect(unique(mpc.gen(:,1)),mpc.bus(g_bus,1));
load = setdiff(mpc.bus(:,1),gen);
Positive_bus = intersect(find(mpc.bus(:,3)>0), find(mpc.bus(:,4)>0));
posi_load = intersect(find(mpc.bus(:,2)==1), Positive_bus);
end