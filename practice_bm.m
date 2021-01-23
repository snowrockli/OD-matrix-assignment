%===============铁路线网配流模型================
clc
clear
close all
load data_OD;%OD间需求
distance=xlsread('OD.xlsx','dis');%两地距离
TT_c=xlsread('OD.xlsx','time_c');%普速站点间运行时间矩阵
TT_h=xlsread('OD.xlsx','time_h');%高铁站点间运行时间矩阵
t=1;
t_T=25;
M_punish=13000;
K_punish=100000;%惩罚因子

theita=0.2;%logit效用感知系数
learn=0.9;%学习强度
theita_l_so=1.5;%社会交互水平
lmemory=t_T;%记忆长度
kesei_p=0.03;%票价系数
kesei_t=0.48;%时间系数
kesei_crowd=0.1;%拥挤系数
air_compete=900;%旅行长度阈值，既超过此值，飞机加入竞争
N_agent_length=8;%OD内元胞空间边长
%===========基本参数===========
[num_station,~]=size(OD_demand);%站点数量
h_speed=300;
c_speed=120;
a_speed=800;
%===========线路=====================
%=======高速============
route(1,:)=[1,3,5,8,7];
route(2,:)=[1,3,6,9,0];
route(3,:)=[1,3,5,4,0];
route(4,:)=[1,2,0,0,0];
route(5,:)=[2,5,8,7,0];
route(6,:)=[2,3,6,9,0];
route(7,:)=[4,5,6,9,0];
route(8,:)=[4,5,8,0,0];
route(9,:)=[6,9,8,7,0];
route(10,:)=[6,5,8,7,0];
%=======普速===========
route(11,:)=[1,3,5,4,7];
route(12,:)=[1,3,5,8,7];
route(13,:)=[1,3,6,9,0];
route(14,:)=[1,3,5,4,0];
route(15,:)=[1,2,0,0,0];
route(16,:)=[2,5,4,7,0];
route(17,:)=[2,3,6,9,0];
route(18,:)=[4,5,6,9,0];
route(19,:)=[6,9,8,7,0];
%========差异化高铁单位里程票价=======
% p=[0.925146961036997,0.584317495831686;
%     0.753277451772915,0.594264797961231;
%     0.762977120200327,0.371241739349683;
%     0.674698757073415,0.540307352901685;
%     0.949416062596097,0.540996841860156;
%     0.741498555807319,0.426233598022061;
%     0.979280571240929,0.447458821311605;
%     0.840530106693025,0.321365995092424;
%     0.740342994073744,0.520546032036620;
%     0.818320606580949,0.599360679397732];
p=[0.75,0.25;0.75,0.25;0.75,0.25;0.75,0.25;0.75,0.25;0.75,0.25;0.75,0.25;0.75,0.25;0.75,0.25;0.75,0.25];
% p=[0.702,0.43;
%    0.592,0.35;
%    0.600,0.37;
%    0.480,0.30;
%    0.745,0.46;
%    0.817,0.48;
%    0.660,0.40;
%    0.974,0.50;
%    0.920,0.50;
%    0.920,0.50];
%========高铁超员======================
capacity(1)=1200;capacity(2)=1200;capacity(3)=1200;capacity(4)=1200;capacity(5)=1200;
capacity(6)=1200;capacity(7)=1200;capacity(8)=1200;capacity(9)=1200;capacity(10)=1200;
capacity_r(1:10)=[1.1382 1.0725 1.0793 1.1977 1.0516 1.0340 1.1044 1.1043 1.1117 1.0368];
%=======普铁超员=======================
capacity(11)=1800;capacity(12)=1800;capacity(13)=1800;capacity(14)=1800;capacity(15)=1800;capacity(16)=1800;
capacity(17)=1800;capacity(18)=1800;capacity(19)=1800;
capacity_r(11)=1;capacity_r(12)=1;capacity_r(13)=1;capacity_r(14)=1;capacity_r(15)=1;capacity_r(16)=1;
capacity_r(17)=1;capacity_r(18)=1;capacity_r(19)=1;
%=======普速铁路单位里程票价==========
p_yz=0.1;p_yw=0.2;
%=======民航单位里程票价==============
p_air=1.2;
p_air=p_game(p,p_yz,p_yw,p_air,theita,h_speed,c_speed,a_speed,OD_demand);%民航票价博弈调整
%==============开始迭代===================
while t<=t_T
    for i=1:num_station
        for j=1:num_station
            cell{i,j}.route_num=0;%线路数量
            cell{i,j}.route=[];
            cell{i,j}.route_length=[];
            cell{i,j}.train_travel_time=[];
            cell{i,j}.air_travel_time=[];
            cell{i,j}.train_fare_1=[];
            cell{i,j}.train_fare_2=[];
            cell{i,j}.air_fare=[];
            cell{i,j}.real_route=[];
            cell{i,j}.direction=[];
            cell{i,j}.train_fee_1=[];
            cell{i,j}.train_fee_2=[];
            cell{i,j}.air_fee=[];
            
            %=========确定每个OD间线路数量与名称============
            for k=1:19
                if ismember(i,route(k,:))&&ismember(j,route(k,:))&&(i~=j)
                    cell{i,j}.route_num=cell{i,j}.route_num+1;
                    cell{i,j}.route=[cell{i,j}.route,k];%添加线路名称
                    cell{i,j}.i_locate(cell{i,j}.route_num,1)=find(ismember(route(k,:),i));%定位起点在线路中的位置
                    cell{i,j}.j_locate(cell{i,j}.route_num,1)=find(ismember(route(k,:),j));%定位终点在线路中的位置
                    if cell{i,j}.i_locate(cell{i,j}.route_num,1)<cell{i,j}.j_locate(cell{i,j}.route_num,1)%方向
                        cell{i,j}.direction(cell{i,j}.route_num,1)=1;%上行
                    else
                        cell{i,j}.direction(cell{i,j}.route_num,1)=-1;%下行
                    end
                end
            end
            %==============计算每个OD间的实际线路==============
            cell{i,j}.real_route=zeros(cell{i,j}.route_num,num_station);
            for k=1:cell{i,j}.route_num
                if cell{i,j}.direction(k)==1%上行实际的线路
                    cell{i,j}.route_length(k)=size(route(cell{i,j}.route(k),cell{i,j}.i_locate(k,1):cell{i,j}.j_locate(k,1)),2);
                    cell{i,j}.real_route(k,1:cell{i,j}.route_length(k))=route(cell{i,j}.route(k),cell{i,j}.i_locate(k,1):cell{i,j}.j_locate(k,1));
                elseif cell{i,j}.direction(k)==-1%下行实际的线路
                    cell{i,j}.route_length(k)=size(route(cell{i,j}.route(k),cell{i,j}.j_locate(k,1):cell{i,j}.i_locate(k,1)),2);
                    cell{i,j}.real_route(k,1:cell{i,j}.route_length(k))=route(cell{i,j}.route(k),cell{i,j}.j_locate(k,1):cell{i,j}.i_locate(k,1));
                end
            end
            %============列车运行时间、票价==============
            for k=1:cell{i,j}.route_num
                cell{i,j}.train_travel_time(k)=0;
                non_zero=cell{i,j}.real_route(k,(find(cell{i,j}.real_route(k,:)~=0)));%提取实际线路
                for kk=1:length(non_zero)-1
                    if cell{i,j}.route(k)<=10%高铁
                        cell{i,j}.train_travel_time(k)=cell{i,j}.train_travel_time(k)+TT_h(non_zero(kk),non_zero(kk+1));%高铁运行时间
                    elseif cell{i,j}.route(k)>10%普速
                        cell{i,j}.train_travel_time(k)=cell{i,j}.train_travel_time(k)+TT_c(non_zero(kk),non_zero(kk+1));%普速运行时间
                    end
                end
                if cell{i,j}.route(k)<=10%高铁
                    cell{i,j}.train_fare_1(k)=p(cell{i,j}.route(k),1)*h_speed*cell{i,j}.train_travel_time(k);%一等座票价
                    cell{i,j}.train_fare_2(k)=p(cell{i,j}.route(k),2)*h_speed*cell{i,j}.train_travel_time(k);%二等座票价
                elseif cell{i,j}.route(k)>10%普速
                    cell{i,j}.train_fare_1(k)=p_yw*c_speed*cell{i,j}.train_travel_time(k);%硬卧票价
                    cell{i,j}.train_fare_2(k)=p_yz*c_speed*cell{i,j}.train_travel_time(k);%硬座票价
                end
            end
            cell{i,j}.air_travel_time=distance(i,j)/a_speed;%民航运行时间
            cell{i,j}.air_fare=distance(i,j)*p_air;%民航票价
            %===========计算各种出行方式人数=====================
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if t==1&&i~=j
                %=============民航广义费用==============
                cell{i,j}.air_fee=kesei_t*cell{i,j}.air_travel_time+kesei_p*cell{i,j}.air_fare;
                cell{i,j}.air_fee_b=cell{i,j}.air_fee;
                %=============铁路广义费用===============
                for k=1:cell{i,j}.route_num
                    cell{i,j}.train_fee_1(k)=kesei_t*cell{i,j}.train_travel_time(k)+kesei_p*cell{i,j}.train_fare_1(k);
                    cell{i,j}.train_fee_2(k)=kesei_t*cell{i,j}.train_travel_time(k)+kesei_p*cell{i,j}.train_fare_2(k);
                    cell{i,j}.train_fee_b_1(k)=cell{i,j}.train_fee_1(k);
                    cell{i,j}.train_fee_b_2(k)=cell{i,j}.train_fee_2(k);
                end
            elseif t>1&&i~=j
                %=============民航广义费用==============
                cell{i,j}.air_fee=kesei_t*cell{i,j}.air_travel_time+kesei_p*cell{i,j}.air_fare+theita_l_so*(sum(cell{i,j}.air_q(:))-sum(cell{i,j}.train_q_1(:))-sum(cell{i,j}.train_q_2(:)))/mean(OD_demand(:));
                cell{i,j}.air_fee_b=kesei_t*cell{i,j}.air_travel_time+kesei_p*cell{i,j}.air_fare;
                %=============铁路广义费用===============
                for k=1:cell{i,j}.route_num
                    cell{i,j}.train_fee_1(k)=kesei_t*cell{i,j}.train_travel_time(k)+kesei_p*cell{i,j}.train_fare_1(k)+kesei_crowd*cell{i,j}.train_crowd(k)/capacity(cell{i,j}.route(k))+theita_l_so*(sum(cell{i,j}.train_q_1(:))+sum(cell{i,j}.train_q_2(:))-sum(cell{i,j}.air_q(:)))/mean(OD_demand(:));
                    cell{i,j}.train_fee_2(k)=kesei_t*cell{i,j}.train_travel_time(k)+kesei_p*cell{i,j}.train_fare_2(k)+kesei_crowd*cell{i,j}.train_crowd(k)*capacity_r(cell{i,j}.route(k))/capacity(cell{i,j}.route(k))+theita_l_so*(sum(cell{i,j}.train_q_1(:))+sum(cell{i,j}.train_q_2(:))-sum(cell{i,j}.air_q(:)))/mean(OD_demand(:));
                    cell{i,j}.train_fee_b_1(k)=kesei_t*cell{i,j}.train_travel_time(k)+kesei_p*cell{i,j}.train_fare_1(k)+kesei_crowd*cell{i,j}.train_crowd(k)/capacity(cell{i,j}.route(k));
                    cell{i,j}.train_fee_b_2(k)=kesei_t*cell{i,j}.train_travel_time(k)+kesei_p*cell{i,j}.train_fare_2(k)+kesei_crowd*cell{i,j}.train_crowd(k)*capacity_r(cell{i,j}.route(k))/capacity(cell{i,j}.route(k));
                end
            end

            %==============选择不同方式的人数(期望效用+bm模型)=============================
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cell=multi_agent_bm(t,i,j,cell,air_compete,learn,lmemory,distance,OD_demand,N_agent_length);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if i~=j
                air_q(i,j)=cell{i,j}.train_q_1(1);
            else
                air_q(i,j)=0;
            end 
        end
    end
    %==============记录一些结果(以OD{7,1}为例)===========
    flow_1(t)=cell{1,8}.train_q_1(1);
    flow_2(t)=cell{1,8}.train_q_1(2);
    flow_3(t)=cell{1,8}.train_q_2(1);
    flow_4(t)=cell{1,8}.train_q_2(2);
    flow_5(t)=cell{7,6}.air_q;
    
    %===========检查收敛性=======================
    air_qq(t)=mean(mean(air_q(:,:)));%检查收敛性
    if t==1||t==2
        q_var(t)=0.95;
    else
        q_var(t)=abs(air_qq(t)-air_qq(t-1))/air_qq(t-1);%检查收敛性
    end
    %=============计算每条线路的使用情况===========
    q=zeros(num_station,num_station,19);
    for ii=1:num_station
        for jj=1:num_station
            for k=1:19
                if ii==jj
                    q(ii,jj,k)=0;
                else
                    if isempty(find(cell{ii,jj}.route==k))==1
                        q(ii,jj,k)=0;
                    else
                        q(ii,jj,k)=cell{ii,jj}.train_q_1(find(cell{ii,jj}.route==k))+cell{ii,jj}.train_q_2(find(cell{ii,jj}.route==k));
                    end
                end
            end
        end
    end
    %===========计算每个路段的客流拥挤程度（上行）================
    for k=1:size(route,1)
        for ii=1:size(route,2)-1
            if route(k,ii)*route(k,ii+1)==0
                crowd(k,ii,1)=0;%上行
            else
                if ii==1
                    crowd(k,ii,1)=0;
                    for jj=ii+1:size(route,2)
                        if route(k,jj)~=0
                            crowd(k,ii,1)=crowd(k,ii,1)+q(route(k,ii),route(k,jj),k); 
                        end
                    end
                else
                   q_hou=0;
                   q_qian=0;
                   for jj=ii+1:size(route,2)
                       if route(k,jj)~=0
                           q_hou=q_hou+q(route(k,ii),route(k,jj),k);
                       end
                   end
                   for jj=1:ii
                       if route(k,jj)~=0
                           q_qian=q_qian+q(route(k,jj),route(k,ii),k);
                       end
                   end
                   crowd(k,ii,1)=crowd(k,ii-1,1)+q_hou-q_qian;
                end
            end   
        end
    end
    %===========计算每个路段的拥挤程度（下行）================
    for k=1:size(route,1)
        for ii=size(route,2):-1:2
            if route(k,ii)*route(k,ii-1)==0
                crowd(k,ii-1,2)=0;%下行
            else
                if ii==size(route,2)
                    crowd(k,ii-1,2)=0;
                    for jj=ii-1:-1:1
                        if route(k,jj)~=0
                            crowd(k,ii-1,2)=crowd(k,ii-1,2)+q(route(k,ii),route(k,jj),k);
                        else
                            crowd(k,ii-1,2)=0;
                        end
                    end
                else
                    q_hou=0;
                    q_qian=0;
                    for jj=1:ii-1
                        if route(k,jj)~=0
                            q_hou=q_hou+q(route(k,ii),route(k,jj),k);
                        end
                    end
                    for jj=ii:size(route,2)
                        if route(k,jj)~=0
                            q_qian=q_qian+q(route(k,jj),route(k,ii),k);
                        end
                    end
                    crowd(k,ii-1,2)=crowd(k,ii,2)+q_hou-q_qian;
                end
            end
        end
    end
    %==================计算每个OD间的拥挤度===================
    for i=1:num_station
        for j=1:num_station
            for k=1:cell{i,j}.route_num
                if cell{i,j}.direction(k)==1
                    cell{i,j}.train_crowd(k)=sum(crowd(cell{i,j}.route(k),cell{i,j}.i_locate(k):cell{i,j}.j_locate(k)-1,1));
                elseif cell{i,j}.direction(k)==-1
                    cell{i,j}.train_crowd(k)=sum(crowd(cell{i,j}.route(k),cell{i,j}.j_locate(k):cell{i,j}.i_locate(k)-1,2));
                end
            end
        end
    end
    
    
    t=t+1;
end
%===============计算目标函数=============
%目标1：铁路运输成本最小化
%目标2：广义出行费用最小化
for i=1:num_station
    for j=1:num_station
        if i~=j
            for k=1:cell{i,j}.route_num
                for t=1:10
                    q1_t(i,j,k,t)=OD_demand(i,j)*length(find(cell{i,j}.path(:,:,t_T-t)==k))/(N_agent_length^2);
                    q2_t(i,j,k,t)=OD_demand(i,j)*length(find(cell{i,j}.path(:,:,t_T-t)==cell{i,j}.route_num+k))/(N_agent_length^2);
                end
%                 profit_train(k)=cell{i,j}.train_fare_1(k)*cell{i,j}.train_q_1(k)+cell{i,j}.train_fare_2(k)*cell{i,j}.train_q_2(k);
%                 g_fee_train(k)=cell{i,j}.train_fee_1(k)*cell{i,j}.train_q_1(k)+cell{i,j}.train_fee_2(k)*cell{i,j}.train_q_2(k);
                profit_train(k)=cell{i,j}.train_fare_1(k)*mean(q1_t(i,j,k,:))+cell{i,j}.train_fare_2(k)*mean(q2_t(i,j,k,:));
                g_fee_train(k)=cell{i,j}.train_fee_1(k)*mean(q1_t(i,j,k,:))+cell{i,j}.train_fee_2(k)*mean(q2_t(i,j,k,:));
                g_fee_train_b(k)=cell{i,j}.train_fee_b_1(k)*mean(q1_t(i,j,k,:))+cell{i,j}.train_fee_b_2(k)*mean(q2_t(i,j,k,:));
            end
            profit(i,j)=sum(profit_train(:));%铁路运输成本最小化
            g_fee(i,j)=sum(g_fee_train(:));%出行广义费用最小化
            g_fee_b(i,j)=sum(g_fee_train_b(:));
        else
            profit(i,j)=0;
            g_fee(i,j)=0;
            g_fee_b(i,j)=0;
        end
    end
end
%=======目标函数，总成本最小==============
if mean(mean(profit))<=M_punish
    f=-mean(mean(profit))+mean(mean(g_fee));
elseif mean(mean(profit))>M_punish
    f=-mean(mean(profit))+mean(mean(g_fee))+K_punish*(mean(mean(profit))-M_punish);
end

%plot(2:t_T,q_var(2:t_T));
plot(1:t_T,flow_1(1:t_T),1:t_T,flow_2(1:t_T),1:t_T,flow_3(1:t_T),1:t_T,flow_4(1:t_T),1:t_T,flow_5(1:t_T));
%==========计算最大断面流量==========
max_q_sum=max(max(max(q(:,:,:))));