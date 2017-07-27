function [pd,ps,alpha_org,alpha]=p_tot_div_to_ps_and_pd(p_tot,v_max,f_max,beta_ps,beta_pd,ratio_v,dvfs_stage)
%this function is used to calculate pd and ps, when given the max allowable
%dynamic power value pd_max, current total power budget value p_tot_current
%and the threshold temperature value t_core.
%syms x
%f=p_tot-beta_pd*x^2*(f_max/(1-ratio_v))*(x/v_max-ratio_v)-x*beta_ps;

%f=p_tot_current-x^3*pd_max_per_core-x^(index)*(0.00268279*ps_org_para*t_core^2.*exp(45.1565/t_core));
%v=solve (f);
%v=eval(v);
%A=v==real(v);
%v=v(A);



v=(0.0037*p_tot + ((0.0037*p_tot - 0.0193)^2 + 0.0016)^(1/2) - 0.0193)^(1/3) - 0.1167/(0.0037*p_tot + ((0.0037*p_tot -0.0193)^2 + 0.0016)^(1/2) - 0.0193)^(1/3) + 0.1067;

%pd=alpha^3*pd_max_per_core;
%ps=alpha^(index)*(0.00268279*ps_org_para*t_core^2.*exp(45.1565/t_core));
f=(f_max/(1-ratio_v))*(v/v_max-ratio_v);
alpha_org=f/f_max;
for j=1:size(dvfs_stage,2)
    if alpha_org>=dvfs_stage(j)
        alpha=dvfs_stage(j);
        break
    end
end  
f=f_max*alpha;
v=v_max*(ratio_v+(1-ratio_v)*alpha);
pd=beta_pd*v^2*f;
ps=p_tot-pd;
