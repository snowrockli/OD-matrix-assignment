function cell= multi_agent_bm(t,i,j,cell,air_compete,learn,lmemory,distance,OD_demand,N_agent_length)
%%%%%%%%%%%%%%%%%%%%%%%%小于一定里程不坐飞机%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i~=j&&distance(i,j)<air_compete%小于一定里程不坐飞机
    if t==1
        for ii=1:N_agent_length
            for jj=1:N_agent_length
                for k=1:cell{i,j}.route_num*2
                    randr(k)=rand;
                end
                randxigema=sum(randr(:));
                for k=1:cell{i,j}.route_num*2
                    cell{i,j}.w(ii,jj,k,1)=randr(k)/randxigema;
                    cell{i,j}.ws(ii,jj,k,1)=randr(k)/randxigema;
                    cell{i,j}.wss(ii,jj,k,1)=randr(k)/randxigema;
                    cell{i,j}.PT(ii,jj,k,1)=0;
                end
                cell{i,j}.path(ii,jj,t)=randsrc(1,1,[1:cell{i,j}.route_num*2]);%运输方式
            end
        end
    end
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            if cell{i,j}.path(ii,jj,t)<=cell{i,j}.route_num
                cell{i,j}.ET(ii,jj,t)=cell{i,j}.train_fee_1(cell{i,j}.path(ii,jj,t));%一等座、硬卧经验
                cell{i,j}.PT(ii,jj,cell{i,j}.path(ii,jj,t),t)=cell{i,j}.train_fee_1(cell{i,j}.path(ii,jj,t));%一等座、硬卧感知
            else
                cell{i,j}.ET(ii,jj,t)=cell{i,j}.train_fee_2(cell{i,j}.path(ii,jj,t)-cell{i,j}.route_num);%二等座、硬座经验
                cell{i,j}.PT(ii,jj,cell{i,j}.path(ii,jj,t),t)=cell{i,j}.train_fee_2(cell{i,j}.path(ii,jj,t)-cell{i,j}.route_num);%二等座、硬座感知
            end
        end
    end
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            if t==1
                cell{i,j}.EET(ii,jj,t)=0;%experience
            else
                if t<=lmemory
                    cell{i,j}.EET(ii,jj,t)=sum(cell{i,j}.ET(ii,jj,1:t-1))/(t-1);%所有经验的均值
                else
                    cell{i,j}.EET(ii,jj,t)=sum(cell{i,j}.ET(ii,jj,t-lmemory:t-1))/(lmemory-1);%所有经验的均值
                end
            end
            for k=1:cell{i,j}.route_num*2
                if t==1
                    cell{i,j}.PPT(ii,jj,k,t)=cell{i,j}.PT(ii,jj,k,t);%perception
                else
                    if t<=lmemory
                        cell{i,j}.PPT(ii,jj,k,t)=sum(cell{i,j}.PT(ii,jj,k,1:t))/t;%某一种经验的均值
                    else
                        cell{i,j}.PPT(ii,jj,k,t)=sum(cell{i,j}.PT(ii,jj,k,t-lmemory:t))/lmemory;%某一种经验的均值
                    end
                end
            end
            if t==1
                cell{i,j}.SC(ii,jj,t)=0;
            else
                if cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t)>=0
                    cell{i,j}.SC(ii,jj,t)=(cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t))/abs(max(cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,find(cell{i,j}.PPT(ii,jj,:,t)~=cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t)),t)));%stimulus
                else
                    cell{i,j}.SC(ii,jj,t)=(cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t))/abs(min(cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,find(cell{i,j}.PPT(ii,jj,:,t)~=cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t)),t)));
                end
            end 
        end
    end
    %===========================更新出行方式选择概率=======================
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            for k=1:cell{i,j}.route_num*2
                if k==cell{i,j}.path(ii,jj,t)
                    if cell{i,j}.SC(ii,jj,t)>=0
                        cell{i,j}.w(ii,jj,k,t+1)=cell{i,j}.w(ii,jj,k,t)+(1-cell{i,j}.w(ii,jj,k,t))*learn*cell{i,j}.SC(ii,jj,t);
                    else
                        cell{i,j}.w(ii,jj,k,t+1)=cell{i,j}.w(ii,jj,k,t)+cell{i,j}.w(ii,jj,k,t)*learn*cell{i,j}.SC(ii,jj,t);
                    end
                else
                    if cell{i,j}.SC(ii,jj,t)>=0
                        cell{i,j}.w(ii,jj,k,t+1)=cell{i,j}.w(ii,jj,k,t)-cell{i,j}.w(ii,jj,k,t)*learn*cell{i,j}.SC(ii,jj,t);
                    else
                        cell{i,j}.w(ii,jj,k,t+1)=cell{i,j}.w(ii,jj,k,t)-cell{i,j}.w(ii,jj,k,t)*cell{i,j}.w(ii,jj,cell{i,j}.path(ii,jj,t),t)*learn*cell{i,j}.SC(ii,jj,t)/(1-cell{i,j}.w(ii,jj,cell{i,j}.path(ii,jj,t),t));
                    end
                end
            end
        end
    end
    %==================================运输方式选择概率标准化===============
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            for k=1:cell{i,j}.route_num*2
                cell{i,j}.ws(ii,jj,k,t+1)=(cell{i,j}.w(ii,jj,k,t)-min(cell{i,j}.w(ii,jj,:,t)))/(max(cell{i,j}.w(ii,jj,:,t))-min(cell{i,j}.w(ii,jj,:,t)));
            end
        end
    end
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            for k=1:cell{i,j}.route_num*2
                cell{i,j}.wss(ii,jj,k,t+1)=cell{i,j}.ws(ii,jj,k,t+1)/sum(cell{i,j}.ws(ii,jj,:,t+1));
            end
        end
    end
    %===============铁路客流分配=======================
    for k=1:cell{i,j}.route_num
        %=======一等座、硬卧========
        cell{i,j}.train_q_1(k)=OD_demand(i,j)*mean(mean(cell{i,j}.wss(:,:,k,t+1)));
        %=======二等座、硬座=========
        cell{i,j}.train_q_2(k)=OD_demand(i,j)*mean(mean(cell{i,j}.wss(:,:,cell{i,j}.route_num+k,t+1)));
    end
    cell{i,j}.air_q=0;
    %=================更新出行方式选择==================
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            randroute=rand;
            for k=1:cell{i,j}.route_num*2
                cell{i,j}.cwp(ii,jj,k)=sum(cell{i,j}.wss(ii,jj,1:k,t+1));
            end
            for k=1:cell{i,j}.route_num*2-1
                if randroute<=cell{i,j}.cwp(ii,jj,1)
                    cell{i,j}.path(ii,jj,t+1)=1;
                elseif randroute>cell{i,j}.cwp(ii,jj,k)&&randroute<=cell{i,j}.cwp(ii,jj,k+1)
                    cell{i,j}.path(ii,jj,t+1)=k+1;
                end
            end
            if isnan(cell{i,j}.path(ii,jj,t+1))
                cell{i,j}.path(ii,jj,t+1)=randsrc(1,1,[1:cell{i,j}.route_num*2]);
            end
        end
    end 
%%%%%%%%%%%%%%%%%%%%%%%%%大于一定里程引入民航%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif i~=j&&distance(i,j)>=air_compete%大于一定里程引入民航
    if t==1
        for ii=1:N_agent_length
            for jj=1:N_agent_length
                for k=1:cell{i,j}.route_num*2+1
                    randr(k)=rand;
                end
                randxigema=sum(randr(:));
                for k=1:cell{i,j}.route_num*2+1
                    cell{i,j}.w(ii,jj,k,1)=randr(k)/randxigema;
                    cell{i,j}.ws(ii,jj,k,1)=randr(k)/randxigema;
                    cell{i,j}.wss(ii,jj,k,1)=randr(k)/randxigema;
                    cell{i,j}.PT(ii,jj,k,1)=0;
                end
                cell{i,j}.path(ii,jj,t)=randsrc(1,1,[1:cell{i,j}.route_num*2+1]);%运输方式
            end
        end
    end
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            if cell{i,j}.path(ii,jj,t)<=cell{i,j}.route_num
                cell{i,j}.ET(ii,jj,t)=cell{i,j}.train_fee_1(cell{i,j}.path(ii,jj,t));%一等座、硬卧经验
                cell{i,j}.PT(ii,jj,cell{i,j}.path(ii,jj,t),t)=cell{i,j}.train_fee_1(cell{i,j}.path(ii,jj,t));%一等座、硬卧感知
            elseif cell{i,j}.path(ii,jj,t)>cell{i,j}.route_num&&cell{i,j}.path(ii,jj,t)<=cell{i,j}.route_num*2
                cell{i,j}.ET(ii,jj,t)=cell{i,j}.train_fee_2(cell{i,j}.path(ii,jj,t)-cell{i,j}.route_num);%二等座、硬座经验
                cell{i,j}.PT(ii,jj,cell{i,j}.path(ii,jj,t),t)=cell{i,j}.train_fee_2(cell{i,j}.path(ii,jj,t)-cell{i,j}.route_num);%二等座、硬座感知
            elseif cell{i,j}.path(ii,jj,t)==cell{i,j}.route_num*2+1
                cell{i,j}.ET(ii,jj,t)=cell{i,j}.air_fee;
                cell{i,j}.PT(ii,jj,cell{i,j}.path(ii,jj,t),t)=cell{i,j}.air_fee;
            end
        end
    end
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            if t==1
                cell{i,j}.EET(ii,jj,t)=0;%experience
            else
                if t<=lmemory
                    cell{i,j}.EET(ii,jj,t)=sum(cell{i,j}.ET(ii,jj,1:t-1))/(t-1);%所有经验的均值
                else
                    cell{i,j}.EET(ii,jj,t)=sum(cell{i,j}.ET(ii,jj,t-lmemory:t-1))/(lmemory-1);%所有经验的均值
                end
            end
            for k=1:cell{i,j}.route_num*2+1
                if t==1
                    cell{i,j}.PPT(ii,jj,k,t)=cell{i,j}.PT(ii,jj,k,t);%perception
                else
                    if t<=lmemory
                        cell{i,j}.PPT(ii,jj,k,t)=sum(cell{i,j}.PT(ii,jj,k,1:t))/t;%某一种经验的均值
                    else
                        cell{i,j}.PPT(ii,jj,k,t)=sum(cell{i,j}.PT(ii,jj,k,t-lmemory:t))/lmemory;%某一种经验的均值
                    end
                end
            end
            if t==1
                cell{i,j}.SC(ii,jj,t)=0;
            else
                if cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t)>=0
                    cell{i,j}.SC(ii,jj,t)=(cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t))/abs(max(cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,find(cell{i,j}.PPT(ii,jj,:,t)~=cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t)),t)));%stimulus
                else
                    cell{i,j}.SC(ii,jj,t)=(cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t))/abs(min(cell{i,j}.EET(ii,jj,t)-cell{i,j}.PPT(ii,jj,find(cell{i,j}.PPT(ii,jj,:,t)~=cell{i,j}.PPT(ii,jj,cell{i,j}.path(ii,jj,t),t)),t)));
                end
            end
        end
    end
    %===========================更新出行方式选择概率=======================
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            for k=1:cell{i,j}.route_num*2+1
                if k==cell{i,j}.path(ii,jj,t)
                    if cell{i,j}.SC(ii,jj,t)>=0
                        cell{i,j}.w(ii,jj,k,t+1)=cell{i,j}.w(ii,jj,k,t)+(1-cell{i,j}.w(ii,jj,k,t))*learn*cell{i,j}.SC(ii,jj,t);
                    else
                        cell{i,j}.w(ii,jj,k,t+1)=cell{i,j}.w(ii,jj,k,t)+cell{i,j}.w(ii,jj,k,t)*learn*cell{i,j}.SC(ii,jj,t);
                    end
                else
                    if cell{i,j}.SC(ii,jj,t)>=0
                        cell{i,j}.w(ii,jj,k,t+1)=cell{i,j}.w(ii,jj,k,t)-cell{i,j}.w(ii,jj,k,t)*learn*cell{i,j}.SC(ii,jj,t);
                    else
                        cell{i,j}.w(ii,jj,k,t+1)=cell{i,j}.w(ii,jj,k,t)-cell{i,j}.w(ii,jj,k,t)*cell{i,j}.w(ii,jj,cell{i,j}.path(ii,jj,t),t)*learn*cell{i,j}.SC(ii,jj,t)/(1-cell{i,j}.w(ii,jj,cell{i,j}.path(ii,jj,t),t));
                    end
                end
            end
        end
    end
    %==================================运输方式选择概率标准化===============
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            for k=1:cell{i,j}.route_num*2+1
                cell{i,j}.ws(ii,jj,k,t+1)=(cell{i,j}.w(ii,jj,k,t)-min(cell{i,j}.w(ii,jj,:,t)))/(max(cell{i,j}.w(ii,jj,:,t))-min(cell{i,j}.w(ii,jj,:,t)));
            end
        end
    end
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            for k=1:cell{i,j}.route_num*2+1
                cell{i,j}.wss(ii,jj,k,t+1)=cell{i,j}.ws(ii,jj,k,t+1)/sum(cell{i,j}.ws(ii,jj,:,t+1));
            end
        end
    end
    %===============铁路客流分配========================
    for k=1:cell{i,j}.route_num
        %=======一等座、硬卧========
        cell{i,j}.train_q_1(k)=OD_demand(i,j)*mean(mean(cell{i,j}.wss(:,:,k,t+1)));
        %=======二等座、硬座=========
        cell{i,j}.train_q_2(k)=OD_demand(i,j)*mean(mean(cell{i,j}.wss(:,:,cell{i,j}.route_num+k,t+1)));
    end
    cell{i,j}.air_q=OD_demand(i,j)*mean(mean(cell{i,j}.wss(:,:,2*cell{i,j}.route_num+1,t+1)));
    %=================更新出行方式选择==================
    for ii=1:N_agent_length
        for jj=1:N_agent_length
            randroute=rand;
            for k=1:cell{i,j}.route_num*2+1
                cell{i,j}.cwp(ii,jj,k)=sum(cell{i,j}.wss(ii,jj,1:k,t+1));
            end
            for k=1:cell{i,j}.route_num*2
                if randroute<=cell{i,j}.cwp(ii,jj,1)
                    cell{i,j}.path(ii,jj,t+1)=1;
                elseif randroute>cell{i,j}.cwp(ii,jj,k)&&randroute<=cell{i,j}.cwp(ii,jj,k+1)
                    cell{i,j}.path(ii,jj,t+1)=k+1;
                end
            end
            if isnan(cell{i,j}.path(ii,jj,t+1))
                cell{i,j}.path(ii,jj,t+1)=randsrc(1,1,[1:cell{i,j}.route_num*2+1]);
            end
        end
    end
    
   
     
end



end