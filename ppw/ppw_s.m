function [active_power_map,active_core_map]=ppw_s(A)
v_max=0.8;
f_max=3;
ratio_v=0.4;
core_num=16;
T_active_core=85*ones(core_num,1);
T_idle_core=30*ones(core_num,1);
dvfs_stage=1:-0.05:0;
T_max=200*ones(core_num,1);
p_tot_max=A\T_max;

par=1:-0.2:0.2;
a=0.0188279*0.2;
b=30.1565;
pd_max=0.5*p_tot_max;
beta_pd=max(pd_max./(v_max^2*f_max));
ps_max=p_tot_max-pd_max;
beta_ps=max(ps_max./v_max);
all_core=1:core_num;
active_core_map=zeros(core_num,core_num);
active_power_map=zeros(core_num,1);
for light_core_num=1:core_num
    light_core_com=nchoosek(all_core,light_core_num);
    ppw_second_max=0;
    ppw_org_second_max=0;
    active_core_final=0;
    for k=1:size(light_core_com,1)
        light_core=light_core_com(k,:);
        ppw_first_max=0;
        ppw_org_first_max=0;
        for i=1:light_core_num
            active_core=light_core(i);
            p_tot=zeros(core_num,1);
            p_tot(light_core)=A(light_core,light_core)\T_idle_core(light_core);
            p_tot(active_core)=A(active_core,active_core)\T_active_core(active_core);
            t_core=A*p_tot;
            ps_light=zeros(core_num,1);
            pd_active=0;
            alpha_active=0;
            alpha_active_org=0;
            for z=1:core_num
                if p_tot(z)~=0
                    ps_light(z)=a*t_core(z)*t_core(z)*exp(b/t_core(z));
                end
                [pd,ps,alpha_org,alpha]=p_tot_div_to_ps_and_pd(p_tot(active_core),v_max,f_max,beta_ps,beta_pd,ratio_v,dvfs_stage);
                pd_active=pd;
                ps_light(active_core)=ps;
                alpha_active_org=alpha_org;
                alpha_active=alpha;
            end
            sum_ps=sum(ps_light);
            sum_pd=sum(pd_active);
            sum_power=sum_pd+sum_ps;
            ppw=alpha_active/sum_power;
            ppw_org=alpha_active_org/sum_power;
            if ppw>ppw_first_max || (ppw==ppw_first_max && ppw_org>ppw_org_first_max)
                ppw_first_max=ppw;
                ppw_org_first_max=ppw_org;
                active_core_mid=active_core;
                sum_power_mid=sum_power;
            end
        end
        if ppw_first_max>ppw_second_max || (ppw_first_max==ppw_second_max && ppw_org_first_max>ppw_org_second_max)
            ppw_second_max=ppw_first_max;
            ppw_org_second_max=ppw_org_first_max;
            light_core_current=light_core;
            active_core_final=active_core_mid;
            sum_power_final=sum_power_mid;
        end
    end
    active_core_map(light_core_num,light_core_current)=1;
    active_core_map(light_core_num,active_core_final)=2;
    active_power_map(light_core_num)=sum_power_final;
    
end
% sum_power=sum_pd_all+sum_ps_all;
% for i=1:5
%     w=((1-par(i))*sum_power(1)+par(i)./all_core.*sum_power)./(1-par(i)+par(i)./all_core);
%     ppw(i,:)=sum_alpha_all./w;
% end
% ppw_norm=ppw/max(max(ppw));
% ppw_norm
% figure
% for i=1:5
% plot(all_core,ppw_norm(i,:))
% hold on
% axis([1 16 0 1]);
% end
