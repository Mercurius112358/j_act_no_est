G=importdata('data/16Gmatrix');
C=importdata('data/16Cmatrix');
v_max=0.8;
f_max=3;
ratio_v=0.4;
core_num=16;
pkg_num=3*core_num+12;
tot_num=4*core_num+12;
B=[eye(core_num,core_num);zeros(pkg_num,core_num)];
L=[eye(core_num,core_num),zeros(core_num,pkg_num)];
A=L*(G\B);
T_core=85*ones(core_num,1);
dvfs_stage=1:-0.05:0;
T_max=200*ones(core_num,1);
p_tot_max=A\T_max;
pd_max=0.5*p_tot_max;
ps_max=p_tot_max-pd_max;
beta_pd=max(pd_max./(v_max^2*f_max));
beta_ps=max(ps_max./v_max);
light_core_former=[];
all_core=1:core_num;
sum_alpha_all=zeros(5,core_num);
sum_pd_all=zeros(5,core_num);
sum_ps_all=zeros(5,core_num);
light_core_map=[];
sum_power_seq=0;
par=1:-0.2:0.2;

[active_power_map,active_core_map]=ppw_s(A);
for light_core_num=1:core_num
    light_core_com=nchoosek(all_core,light_core_num);
    sum_alpha_max=[0,0,0,0,0];
    sum_alpha_max_org=[0,0,0,0,0];
    sum_pd_max=[0,0,0,0,0];
    sum_ps_max=[0,0,0,0,0];
    ppw_max=[0,0,0,0,0];
    ppw_org_max=[0,0,0,0,0];
    light_core_current=zeros(5,core_num);
    for k=1:size(light_core_com,1)
        light_core=light_core_com(k,:);
        all_core_now=all_core;
        all_core_now(light_core)=0;
        dark_core=all_core_now(all_core_now~=0);
        p_tot_dark=zeros(core_num,1);
        p_tot_dark(light_core)=A(light_core,light_core)\T_core(light_core);
        t_core_dark=A*p_tot_dark;
        ps_dark=zeros(core_num,1);
        pd_dark=zeros(core_num,1);
        alpha_dark_org=zeros(core_num,1);
        alpha_dark=zeros(core_num,1);
        for i=1:core_num
            if p_tot_dark(i)~=0
                t_core=t_core_dark(i);
                p_tot_current=p_tot_dark(i);
                [pd,ps,alpha_org,alpha]=p_tot_div_to_ps_and_pd(p_tot_current,v_max,f_max,beta_ps,beta_pd,ratio_v,dvfs_stage);
                pd_dark(i)=pd;
                ps_dark(i)=ps;
                alpha_dark_org(i)=alpha_org;
                alpha_dark(i)=alpha;
            end
        end
        sum_alpha_org=sum(alpha_dark_org);
        sum_alpha=sum(alpha_dark);
        sum_pd=sum(pd_dark);
        sum_ps=sum(ps_dark);
        sum_power=sum_pd+sum_ps;
        for j=1:5
%             if light_core_num==1 && sum_power>sum_power_seq
%                 sum_power_seq=sum_power;
%             end
%             w=((1-par(j))*sum_power_seq+par(j)/light_core_num*sum_power)/(1-par(j)+par(j)/light_core_num);
            w=((1-par(j))*active_power_map(light_core_num)+par(j)/light_core_num*sum_power)/(1-par(j)+par(j)/light_core_num);
            ppw=sum_alpha/w;
            ppw_org=sum_alpha_org/w;
            if ppw>ppw_max(j) || (ppw==ppw_max(j) && ppw_org>ppw_org_max(j))
                light_core_current(j,:)=0;
                ppw_org_max(j)=ppw_org;
                ppw_max(j)=ppw;
                sum_alpha_max(j)=sum_alpha;
                sum_pd_max(j)=sum_pd;
                sum_ps_max(j)=sum_ps;
                light_core_current(j,light_core)=1;
            end
        end   
    end
    for i=1:5
        sum_alpha_all(i,light_core_num)=sum_alpha_max(i);
        sum_pd_all(i,light_core_num)=sum_pd_max(i);
        sum_ps_all(i,light_core_num)=sum_ps_max(i);
        light_core_map(:,:,light_core_num)=light_core_current;
        ppw_all(i,light_core_num)=ppw_max(i);
    end
end

% sum_power_all=sum_pd_all+sum_ps_all;
% sum_power_all
% for i=1:5
%     sum_power_all(i,:)=sum_power_all(i,:)./all_core;
% end
% figure
% for i=1:5
%     plot(all_core,sum_power_all(i,:))
%     hold on
% end
% previous section is implemented to prove sum_power_1 > sum_power_n/n


ppw_all_norm=ppw_all/max(max(ppw_all));
ppw_all_norm
figure
for i=1:5
plot(all_core,ppw_all_norm(i,:))
hold on
axis([1 16 0 1]);
end