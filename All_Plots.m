close all
clear all
warning off
clc

%% Customer Behavior Comparison for All Pricing Methodologies
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_avcons_S.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_avcons_C.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_avcons_U.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_avcons_W.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    % Titles
    title1 = "Average Consumption per Capita: Private Information - Pricing: " + num2str(mod(i-1,3));
    title2 = "Average Consumption per Capita: Broadcast - Pricing: " + num2str(mod(i-1,3));
    if i / 4 < 1
        title1 = title1 + " (Graph Nodes)";
        title2 = title2 + " (Graph Nodes)";
    else
        title1 = title1 + " (Neighborhoods)";
        title2 = title2 + " (Neighborhoods)";
    end
    % Broadcast Check
    average_consumption_S(2,1) = average_consumption_S(1,1);
    average_consumption_U(2,1) = average_consumption_U(1,1);
    average_consumption_W(2,1) = average_consumption_W(1,1);
    average_consumption_C(2,1) = average_consumption_C(1,1);
    % Figures
    figure(2*i-1);
    hold on
    plot(1:1:Hours,average_consumption_S(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_U(1,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_W(1,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_C(1,:)/N,'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Consumption per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $S, $\Omega = $P','$\Gamma = $U, $\Omega = $P',...
        '$\Gamma = $W, $\Omega = $P','$\Gamma = $CF, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(2*i);
    hold on
    plot(1:1:Hours,average_consumption_S(2,:)/N,'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_U(2,:)/N,'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_W(2,:)/N,'s-','Color',ccc(6,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_C(2,:)/N,'s-','Color',ccc(8,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $S, $\Omega = $B','$\Gamma = $U, $\Omega = $B',...
        '$\Gamma = $W, $\Omega = $B','$\Gamma = $CF, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
end
toc;

%% Pricing Policy Comparison for All Customer Behaviors
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_avcons_S.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_avcons_C.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_avcons_U.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_avcons_W.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    % Broadcast Check
    average_consumption_S(2,1) = average_consumption_S(1,1);
    average_consumption_U(2,1) = average_consumption_U(1,1);
    average_consumption_W(2,1) = average_consumption_W(1,1);
    average_consumption_C(2,1) = average_consumption_C(1,1);
    %
    for j = 1:2
        % Titles
        if j == 1
             title4 = "Average Consumption per Capita: Private Information - ";
             no = 1;
        else
             title4 = "Average Consumption per Capita: Broadcast - ";
             no = 2;
        end
        title1 = title4 + "Selfish";
        title2 = title4 + "Altruistic";
        title3 = title4 + "Welfare Maximizer";
        title4 = title4 + "Community Friendly";
        if i / 4 < 1
            title1 = title1 + " (Graph Nodes)";
            title2 = title2 + " (Graph Nodes)";
            title3 = title3 + " (Graph Nodes)";
            title4 = title4 + " (Graph Nodes)";
            base_no = 12;
        else
            title1 = title1 + " (Neighborhoods)";
            title2 = title2 + " (Neighborhoods)";
            title3 = title3 + " (Neighborhoods)";
            title4 = title4 + " (Neighborhoods)";
            base_no = 20;
        end
        % Figures
        figure(base_no+no);
        hold on
        plot(1:1:Hours,average_consumption_S(no,:)/N,'s-','Color',ccc(i,:),'LineWidth',2)
        grid on
        title(title1)
        ylabel('Average Consumption per Capita (kWh)')
        xlabel('Hours')
        lll = legend('$p^0$','$p^1$','$p^2$');
        set(lll,'Interpreter','Latex');
        set(lll,'FontSize',14);
        hold off
        
        figure(base_no+no+2);
        hold on
        plot(1:1:Hours,average_consumption_U(no,:)/N,'s-','Color',ccc(i,:),'LineWidth',2)
        grid on
        title(title2)
        ylabel('Average Consumption per Capita (kWh)')
        xlabel('Hours')
        lll = legend('$p^0$','$p^1$','$p^2$');
        set(lll,'Interpreter','Latex');
        set(lll,'FontSize',14);
        hold off
        
        figure(base_no+no+4);
        hold on
        plot(1:1:Hours,average_consumption_W(no,:)/N,'s-','Color',ccc(i,:),'LineWidth',2)
        grid on
        title(title3)
        ylabel('Average Consumption per Capita (kWh)')
        xlabel('Hours')
        lll = legend('$p^0$','$p^1$','$p^2$');
        set(lll,'Interpreter','Latex');
        set(lll,'FontSize',14);
        hold off
        
        figure(base_no+no+6);
        hold on
        plot(1:1:Hours,average_consumption_C(no,:)/N,'s-','Color',ccc(i,:),'LineWidth',2)
        grid on
        title(title4)
        ylabel('Average Consumption per Capita (kWh)')
        xlabel('Hours')
        lll = legend('$p^0$','$p^1$','$p^2$');
        set(lll,'Interpreter','Latex');
        set(lll,'FontSize',14);
        hold off
    end
end
toc;

%% Pricing Policy Comparison According to Revenue and Cost
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_revenue_container_S.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_revenue_container_C.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_revenue_container_U.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_revenue_container_W.mat";
    path5 = "Results_Graph/sim" + num2str(i) + "_cost_container_S.mat";
    path6 = "Results_Graph/sim" + num2str(i) + "_cost_container_C.mat";
    path7 = "Results_Graph/sim" + num2str(i) + "_cost_container_U.mat";
    path8 = "Results_Graph/sim" + num2str(i) + "_cost_container_W.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    load(path5)
    load(path6)
    load(path7)
    load(path8)
    % Data Preparation
    average_revenue_S = reshape(mean(revenue_container_S),[3,6]);
    average_revenue_U = reshape(mean(revenue_container_U),[3,6]);
    average_revenue_W = reshape(mean(revenue_container_W),[3,6]);
    average_revenue_C = reshape(mean(revenue_container_C),[3,6]);
    average_pay_S = mean(cost_container_S);
    average_pay_S = reshape(average_pay_S(:,:,:,1),[3,6]);
    average_pay_S = average_revenue_S - average_pay_S;
    average_pay_U = mean(cost_container_U);
    average_pay_U = reshape(average_pay_U(:,:,:,1),[3,6]);
    average_pay_U = average_revenue_U - average_pay_U;
    average_pay_W = mean(cost_container_W);
    average_pay_W = reshape(average_pay_W(:,:,:,1),[3,6]);
    average_pay_W = average_revenue_W - average_pay_W;
    average_pay_C = mean(cost_container_C);
    average_pay_C = reshape(average_pay_C(:,:,:,1),[3,6]);
    average_pay_C = average_revenue_C - average_pay_C;
    % Broadcast Check
    average_revenue_S(2,1) = average_revenue_S(1,1);
    average_revenue_U(2,1) = average_revenue_U(1,1);
    average_revenue_W(2,1) = average_revenue_W(1,1);
    average_revenue_C(2,1) = average_revenue_C(1,1);
    average_pay_S(2,1) = average_pay_S(1,1);
    average_pay_U(2,1) = average_pay_U(1,1);
    average_pay_W(2,1) = average_pay_W(1,1);
    average_pay_C(2,1) = average_pay_C(1,1);
    %
    for j = 1:2
        % Titles
        if j == 1
             title4 = "Revenue and Profit: Private Information - ";
             no = 1;
        else
             title4 = "Revenue and Profit: Broadcast - ";
             no = 2;
        end
        title1 = title4 + "Selfish";
        title2 = title4 + "Altruistic";
        title3 = title4 + "Welfare Maximizer";
        title4 = title4 + "Community Friendly";
        if i / 4 < 1
            title1 = title1 + " (Graph Nodes)";
            title2 = title2 + " (Graph Nodes)";
            title3 = title3 + " (Graph Nodes)";
            title4 = title4 + " (Graph Nodes)";
            base_no = 28;
        else
            title1 = title1 + " (Neighborhoods)";
            title2 = title2 + " (Neighborhoods)";
            title3 = title3 + " (Neighborhoods)";
            title4 = title4 + " (Neighborhoods)";
            base_no = 36;
        end
        % Figures
        figure(base_no+no);
        hold on
        plot(1:1:Hours,average_revenue_S(no,:)/100,'s-','Color',ccc(mod(i,3)+1,:),'LineWidth',2)
        plot(1:1:Hours,average_pay_S(no,:)/100,'s-','Color',ccc(mod(i,3)+4,:),'LineWidth',2)
        grid on
        title(title1)
        ylabel('Revenue and Profit (Dollar)')
        xlabel('Hours')
        lll = legend('$p^0$, Revenue','$p^0$, Profit',...
            '$p^1$, Revenue','$p^1$, Profit',...
            '$p^2$, Revenue','$p^2$, Profit');
        set(lll,'Interpreter','Latex');
        set(lll,'FontSize',14);
        hold off
        
        figure(base_no+no+2);
        hold on
        plot(1:1:Hours,average_revenue_U(no,:)/100,'s-','Color',ccc(mod(i,3)+1,:),'LineWidth',2)
        plot(1:1:Hours,average_pay_U(no,:)/100,'s-','Color',ccc(mod(i,3)+4,:),'LineWidth',2)
        grid on
        title(title2)
        ylabel('Revenue and Profit (Dollar)')
        xlabel('Hours')
        lll = legend('$p^0$, Revenue','$p^0$, Profit',...
            '$p^1$, Revenue','$p^1$, Profit',...
            '$p^2$, Revenue','$p^2$, Profit');
        set(lll,'Interpreter','Latex');
        set(lll,'FontSize',14);
        hold off
        
        figure(base_no+no+4);
        hold on
        plot(1:1:Hours,average_revenue_W(no,:)/100,'s-','Color',ccc(mod(i,3)+1,:),'LineWidth',2)
        plot(1:1:Hours,average_pay_W(no,:)/100,'s-','Color',ccc(mod(i,3)+4,:),'LineWidth',2)
        grid on
        title(title3)
        ylabel('Revenue and Profit (Dollar)')
        xlabel('Hours')
        lll = legend('$p^0$, Revenue','$p^0$, Profit',...
            '$p^1$, Revenue','$p^1$, Profit',...
            '$p^2$, Revenue','$p^2$, Profit');
        set(lll,'Interpreter','Latex');
        set(lll,'FontSize',14);
        hold off
        
        figure(base_no+no+6);
        hold on
        plot(1:1:Hours,average_revenue_C(no,:)/100,'s-','Color',ccc(mod(i,3)+1,:),'LineWidth',2)
        plot(1:1:Hours,average_pay_C(no,:)/100,'s-','Color',ccc(mod(i,3)+4,:),'LineWidth',2)
        grid on
        title(title4)
        ylabel('Revenue and Profit (Dollar)')
        xlabel('Hours')
        lll = legend('$p^0$, Revenue','$p^0$, Profit',...
            '$p^1$, Revenue','$p^1$, Profit',...
            '$p^2$, Revenue','$p^2$, Profit');
        set(lll,'Interpreter','Latex');
        set(lll,'FontSize',14);
        hold off
    end
end
toc;

%% Hybrid Model Comparison (1) for All Pricing Methodologies
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_avcons_H1.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_avcons_H2.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_avcons_H3.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_avcons_H4.mat";
    path5 = "Results_Graph/sim" + num2str(i) + "_avcons_H5.mat";
    path6 = "Results_Graph/sim" + num2str(i) + "_avcons_H6.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    load(path5)
    load(path6)
    % Titles
    title1 = "Average Consumption per Capita - Private Information - Pricing: " + num2str(mod(i-1,3));
    title2 = "Average Consumption per Capita - Broadcast - Pricing: " + num2str(mod(i-1,3));
    if i / 4 < 1
        title1 = title1 + " (Graph Nodes)";
        title2 = title2 + " (Graph Nodes)";
    else
        title1 = title1 + " (Neighborhoods)";
        title2 = title2 + " (Neighborhoods)";
    end
    % Broadcast Check
    average_consumption_H1(2,1) = average_consumption_H1(1,1);
    average_consumption_H2(2,1) = average_consumption_H2(1,1);
    average_consumption_H3(2,1) = average_consumption_H3(1,1);
    average_consumption_H4(2,1) = average_consumption_H4(1,1);
    average_consumption_H5(2,1) = average_consumption_H5(1,1);
    average_consumption_H6(2,1) = average_consumption_H6(1,1);
    % Figures
    figure(4*i+41);
    hold on
    plot(1:1:Hours,average_consumption_H1(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H2(1,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H3(1,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Consumption per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H1, $\Omega = $P','$\Gamma = $H2, $\Omega = $P',...
        '$\Gamma = $H3, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(4*i+42);
    hold on
    plot(1:1:Hours,average_consumption_H1(2,:)/N,'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H2(2,:)/N,'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H3(2,:)/N,'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H1, $\Omega = $B','$\Gamma = $H2, $\Omega = $B',...
        '$\Gamma = $H3, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(4*i+43);
    hold on
    plot(1:1:Hours,average_consumption_H4(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H5(1,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H6(1,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Consumption per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H4, $\Omega = $P','$\Gamma = $H5, $\Omega = $P',...
        '$\Gamma = $H6, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(4*i+44);
    hold on
    plot(1:1:Hours,average_consumption_H4(2,:)/N,'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H5(2,:)/N,'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H6(2,:)/N,'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H4, $\Omega = $B','$\Gamma = $H5, $\Omega = $B',...
        '$\Gamma = $H6, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
end
toc;

%% Hybrid Model Comparison (2) for All Pricing Methodologies
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_avcons_H1.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_avcons_H2.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_avcons_H3.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_avcons_H4.mat";
    path5 = "Results_Graph/sim" + num2str(i) + "_avcons_H5.mat";
    path6 = "Results_Graph/sim" + num2str(i) + "_avcons_H6.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    load(path5)
    load(path6)
    % Titles
    title1 = "Average Consumption per Capita - Pricing: " + num2str(mod(i-1,3));
    if i / 4 < 1
        title1 = title1 + " (Graph Nodes)";
    else
        title1 = title1 + " (Neighborhoods)";
    end
    % Broadcast Check
    average_consumption_H1(2,1) = average_consumption_H1(1,1);
    average_consumption_H2(2,1) = average_consumption_H2(1,1);
    average_consumption_H3(2,1) = average_consumption_H3(1,1);
    average_consumption_H4(2,1) = average_consumption_H4(1,1);
    average_consumption_H5(2,1) = average_consumption_H5(1,1);
    average_consumption_H6(2,1) = average_consumption_H6(1,1);
    % Figures
    figure(3*i+66);
    hold on
    plot(1:1:Hours,average_consumption_H1(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H4(1,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H1(2,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H4(2,:)/N,'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H1, $\Omega = $P','$\Gamma = $H4, $\Omega = $P',...
        '$\Gamma = $H1, $\Omega = $B','$\Gamma = $H4, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(3*i+67);
    hold on
    plot(1:1:Hours,average_consumption_H2(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H5(1,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H2(2,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H5(2,:)/N,'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H2, $\Omega = $P','$\Gamma = $H5, $\Omega = $P',...
        '$\Gamma = $H2, $\Omega = $B','$\Gamma = $H5, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(3*i+68);
    hold on
    plot(1:1:Hours,average_consumption_H3(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H6(1,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H3(2,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H6(2,:)/N,'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H3, $\Omega = $P','$\Gamma = $H6, $\Omega = $P',...
        '$\Gamma = $H3, $\Omega = $B','$\Gamma = $H6, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
end
toc;

%% Hybrid Model Comparison (3) for All Pricing Methodologies
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_avcons_H7.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_avcons_H8.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_avcons_H9.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_avcons_H10.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    % Titles
    title1 = "Average Consumption per Capita: Private Information - Pricing: " + num2str(mod(i-1,3));
    title2 = "Average Consumption per Capita: Broadcast - Pricing: " + num2str(mod(i-1,3));
    if i / 4 < 1
        title1 = title1 + " (Graph Nodes)";
        title2 = title2 + " (Graph Nodes)";
    else
        title1 = title1 + " (Neighborhoods)";
        title2 = title2 + " (Neighborhoods)";
    end
    % Broadcast Check
    average_consumption_H7(2,1) = average_consumption_H7(1,1);
    average_consumption_H8(2,1) = average_consumption_H8(1,1);
    average_consumption_H9(2,1) = average_consumption_H9(1,1);
    average_consumption_H10(2,1) = average_consumption_H10(1,1);
    % Figures
    figure(2*i+85);
    hold on
    plot(1:1:Hours,average_consumption_H7(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H8(1,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H9(1,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H10(1,:)/N,'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H7, $\Omega = $P','$\Gamma = $H8, $\Omega = $P',...
        '$\Gamma = $H9, $\Omega = $P','$\Gamma = $H10, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(2*i+86);
    hold on
    plot(1:1:Hours,average_consumption_H7(2,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H8(2,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H9(2,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_H10(2,:)/N,'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H7, $\Omega = $B','$\Gamma = $H8, $\Omega = $B',...
        '$\Gamma = $H9, $\Omega = $B','$\Gamma = $H10, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
end
toc;

%% Hybrid Model Comparison (4) for All Pricing Methodologies
clear
ccc = colormap(hsv(8));
Hours = 24;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_avcons_24.mat";
    % Load Results
    load(path1)
    % Titles
    title1 = "Average Consumption per Capita: Private Information - Pricing: " + num2str(mod(i-1,3));
    title2 = "Average Consumption per Capita: Broadcast - Pricing: " + num2str(mod(i-1,3));
    if i / 4 < 1
        title1 = title1 + " (Graph Nodes)";
        title2 = title2 + " (Graph Nodes)";
    else
        title1 = title1 + " (Neighborhoods)";
        title2 = title2 + " (Neighborhoods)";
    end
    % Broadcast Check
    average_consumption_24(2,1) = average_consumption_24(1,1);
    average_consumption_24(4,1) = average_consumption_24(3,1);
    average_consumption_24(6,1) = average_consumption_24(5,1);
    average_consumption_24(8,1) = average_consumption_24(7,1);
    % Figures
    figure(2*i+97);
    hold on
    plot(1:1:Hours,average_consumption_24(1,:)/N,'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_24(3,:)/N,'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_24(5,:)/N,'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_24(7,:)/N,'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H7, $\Omega = $P','$\Gamma = $H8, $\Omega = $P',...
        '$\Gamma = $H9, $\Omega = $P','$\Gamma = $H10, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(2*i+98);
    hold on
    plot(1:1:Hours,average_consumption_24(2,:)/N,'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_24(4,:)/N,'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_24(6,:)/N,'s-','Color',ccc(6,:),'LineWidth',2)
    plot(1:1:Hours,average_consumption_24(8,:)/N,'s-','Color',ccc(8,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Consumption per Capita (kWh)')
    xlabel('Hours')
    lll = legend('$\Gamma = $H7, $\Omega = $B','$\Gamma = $H8, $\Omega = $B',...
        '$\Gamma = $H9, $\Omega = $B','$\Gamma = $H10, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
end
toc;

%% Utility Comparison (1) for Hybrid Models
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_utility_H1.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_utility_H2.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_utility_H3.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_utility_H4.mat";
    path5 = "Results_Graph/sim" + num2str(i) + "_utility_H5.mat";
    path6 = "Results_Graph/sim" + num2str(i) + "_utility_H6.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    load(path5)
    load(path6)
    % Data Preparation
    average_utility_H1 = mean(utility_H1);
    average_utility_H1 = reshape(average_utility_H1,[3,100,6]);
    average_utility_H1 = reshape(mean(average_utility_H1,2),[3,6]);
    average_utility_H2 = mean(utility_H2);
    average_utility_H2 = reshape(average_utility_H2,[3,100,6]);
    average_utility_H2 = reshape(mean(average_utility_H2,2),[3,6]);
    average_utility_H3 = mean(utility_H3);
    average_utility_H3 = reshape(average_utility_H3,[3,100,6]);
    average_utility_H3 = reshape(mean(average_utility_H3,2),[3,6]);
    average_utility_H4 = mean(utility_H4);
    average_utility_H4 = reshape(average_utility_H4,[3,100,6]);
    average_utility_H4 = reshape(mean(average_utility_H4,2),[3,6]);
    average_utility_H5 = mean(utility_H5);
    average_utility_H5 = reshape(average_utility_H5,[3,100,6]);
    average_utility_H5 = reshape(mean(average_utility_H5,2),[3,6]);
    average_utility_H6 = mean(utility_H6);
    average_utility_H6 = reshape(average_utility_H6,[3,100,6]);
    average_utility_H6 = reshape(mean(average_utility_H6,2),[3,6]);
    % Broadcast Check
    average_utility_H1(2,1) = average_utility_H1(1,1);
    average_utility_H2(2,1) = average_utility_H2(1,1);
    average_utility_H3(2,1) = average_utility_H3(1,1);
    average_utility_H4(2,1) = average_utility_H4(1,1);
    average_utility_H5(2,1) = average_utility_H5(1,1);
    average_utility_H6(2,1) = average_utility_H6(1,1);
    % Titles
    title1 = "Average Utility per Capita: Private Information - Pricing: " + num2str(mod(i-1,3));
    title2 = "Average Utility per Capita: Broadcast - Pricing: " + num2str(mod(i-1,3));
    if i / 4 < 1
        title1 = title1 + " (Graph Nodes)";
        title2 = title2 + " (Graph Nodes)";
    else
        title1 = title1 + " (Neighborhoods)";
        title2 = title2 + " (Neighborhoods)";
    end
    % Figures
    figure(4*i+107);
    hold on
    plot(1:1:Hours,average_utility_H1(1,:),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H2(1,:),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H3(1,:),'s-','Color',ccc(5,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H1, $\Omega = $P','$\Gamma = $H2, $\Omega = $P',...
        '$\Gamma = $H3, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(4*i+108);
    hold on
    plot(1:1:Hours,average_utility_H1(2,:),'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H2(2,:),'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H3(2,:),'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H1, $\Omega = $B','$\Gamma = $H2, $\Omega = $B',...
        '$\Gamma = $H3, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(4*i+109);
    hold on
    plot(1:1:Hours,average_utility_H4(1,:),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H5(1,:),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H6(1,:),'s-','Color',ccc(5,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H4, $\Omega = $P','$\Gamma = $H5, $\Omega = $P',...
        '$\Gamma = $H6, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(4*i+110);
    hold on
    plot(1:1:Hours,average_utility_H4(2,:),'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H5(2,:),'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H6(2,:),'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H4, $\Omega = $B','$\Gamma = $H5, $\Omega = $B',...
        '$\Gamma = $H6, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
end
toc;

%% Utility Comparison (2) for Hybrid Models
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_utility_H1.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_utility_H2.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_utility_H3.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_utility_H4.mat";
    path5 = "Results_Graph/sim" + num2str(i) + "_utility_H5.mat";
    path6 = "Results_Graph/sim" + num2str(i) + "_utility_H6.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    load(path5)
    load(path6)
    % Data Preparation
    average_utility_H1 = mean(utility_H1);
    average_utility_H1 = reshape(average_utility_H1,[3,100,6]);
    average_utility_H1 = reshape(mean(average_utility_H1,2),[3,6]);
    average_utility_H2 = mean(utility_H2);
    average_utility_H2 = reshape(average_utility_H2,[3,100,6]);
    average_utility_H2 = reshape(mean(average_utility_H2,2),[3,6]);
    average_utility_H3 = mean(utility_H3);
    average_utility_H3 = reshape(average_utility_H3,[3,100,6]);
    average_utility_H3 = reshape(mean(average_utility_H3,2),[3,6]);
    average_utility_H4 = mean(utility_H4);
    average_utility_H4 = reshape(average_utility_H4,[3,100,6]);
    average_utility_H4 = reshape(mean(average_utility_H4,2),[3,6]);
    average_utility_H5 = mean(utility_H5);
    average_utility_H5 = reshape(average_utility_H5,[3,100,6]);
    average_utility_H5 = reshape(mean(average_utility_H5,2),[3,6]);
    average_utility_H6 = mean(utility_H6);
    average_utility_H6 = reshape(average_utility_H6,[3,100,6]);
    average_utility_H6 = reshape(mean(average_utility_H6,2),[3,6]);
    % Broadcast Check
    average_utility_H1(2,1) = average_utility_H1(1,1);
    average_utility_H2(2,1) = average_utility_H2(1,1);
    average_utility_H3(2,1) = average_utility_H3(1,1);
    average_utility_H4(2,1) = average_utility_H4(1,1);
    average_utility_H5(2,1) = average_utility_H5(1,1);
    average_utility_H6(2,1) = average_utility_H6(1,1);
    % Titles
    title1 = "Average Utility per Capita - Pricing: " + num2str(mod(i-1,3));
    if i / 4 < 1
        title1 = title1 + " (Graph Nodes)";
    else
        title1 = title1 + " (Neighborhoods)";
    end
    % Figures
    figure(3*i+132);
    hold on
    plot(1:1:Hours,average_utility_H1(1,:),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H4(1,:),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H1(2,:),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H4(2,:),'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H1, $\Omega = $P','$\Gamma = $H4, $\Omega = $P',...
        '$\Gamma = $H1, $\Omega = $B','$\Gamma = $H4, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(3*i+133);
    hold on
    plot(1:1:Hours,average_utility_H2(1,:),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H5(1,:),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H2(2,:),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H5(2,:),'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H2, $\Omega = $P','$\Gamma = $H5, $\Omega = $P',...
        '$\Gamma = $H2, $\Omega = $B','$\Gamma = $H5, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(3*i+134);
    hold on
    plot(1:1:Hours,average_utility_H3(1,:),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H6(1,:),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H3(2,:),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,average_utility_H6(2,:),'s-','Color',ccc(7,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('$\Gamma = $H3, $\Omega = $P','$\Gamma = $H6, $\Omega = $P',...
        '$\Gamma = $H3, $\Omega = $B','$\Gamma = $H6, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
end
toc;

%% Utility Comparison (3) for Hybrid Models
clear
ccc = colormap(hsv(8));
Hours = 6;
N = 100;
%
tic;
%
run_time = 50;
hybrid = [0.8, 0.1, 0.1, 0; 0.1, 0.8, 0.1, 0; 0.1, 0.1, 0.8, 0; 0, 0.1, 0.1, 0.8;...
    0, 0.8, 0.1, 0.1; 0, 0.1, 0.8, 0.1; 0.20, 0.30, 0.15, 0.35; 0.15, 0.35, 0.20, 0.30;...
    0.25, 0.35, 0.15, 0.35; 0.20, 0.30, 0.20, 0.30];
hybrid = hybrid * N;
all_customers = zeros(6,50,3,100,6);
for hyb = 1:6
    customers = zeros(1,N);
    sp = 1;
    for nums = 1:4
        if sp > 100
            break
        end
        customers(sp:sp+hybrid(hyb,nums)-1) = nums;
        sp = sp+hybrid(hyb,nums);
    end
    for run = 1:run_time
        rng(run);
        customers = customers(randperm(N));
        all_customers(hyb,run,1,:,1) = customers;
    end
end
for i = 1:3
    for j = 1:6
        all_customers(:,:,i,:,j) = all_customers(:,:,1,:,1);
    end
end
%
for i = 1:6
    % Paths
    path1 = "Results_Graph/sim" + num2str(i) + "_utility_H1.mat";
    path2 = "Results_Graph/sim" + num2str(i) + "_utility_H2.mat";
    path3 = "Results_Graph/sim" + num2str(i) + "_utility_H3.mat";
    path4 = "Results_Graph/sim" + num2str(i) + "_utility_H4.mat";
    path5 = "Results_Graph/sim" + num2str(i) + "_utility_H5.mat";
    path6 = "Results_Graph/sim" + num2str(i) + "_utility_H6.mat";
    % Load Results
    load(path1)
    load(path2)
    load(path3)
    load(path4)
    load(path5)
    load(path6)
    % Data Preparation
    average_utility_H1 = zeros(3,3,6);
    raw_utility_H1 = utility_H1 .* reshape(all_customers(1,:,:,:,:) == 1,[50,3,100,6]);
    average_utility_H1(1,:,:) = sum(raw_utility_H1,[1,3]); % S
    average_utility_H1(1,:,:) = average_utility_H1(1,:,:) / run_time / 80;
    raw_utility_H1 = utility_H1 .* reshape(all_customers(1,:,:,:,:) == 2,[50,3,100,6]);
    average_utility_H1(2,:,:) = sum(raw_utility_H1,[1,3]); % U
    average_utility_H1(2,:,:) = average_utility_H1(2,:,:) / run_time / 10;
    raw_utility_H1 = utility_H1 .* reshape(all_customers(1,:,:,:,:) == 3,[50,3,100,6]);
    average_utility_H1(3,:,:) = sum(raw_utility_H1,[1,3]); % W
    average_utility_H1(3,:,:) = average_utility_H1(3,:,:) / run_time / 10;
    
    average_utility_H2 = zeros(3,3,6);
    raw_utility_H2 = utility_H2 .* reshape(all_customers(2,:,:,:,:) == 1,[50,3,100,6]);
    average_utility_H2(1,:,:) = sum(raw_utility_H2,[1,3]); % S
    average_utility_H2(1,:,:) = average_utility_H2(1,:,:) / run_time / 10;
    raw_utility_H2 = utility_H2 .* reshape(all_customers(2,:,:,:,:) == 2,[50,3,100,6]);
    average_utility_H2(2,:,:) = sum(raw_utility_H2,[1,3]); % U
    average_utility_H2(2,:,:) = average_utility_H2(2,:,:) / run_time / 80;
    raw_utility_H2 = utility_H2 .* reshape(all_customers(2,:,:,:,:) == 3,[50,3,100,6]);
    average_utility_H2(3,:,:) = sum(raw_utility_H2,[1,3]); % W
    average_utility_H2(3,:,:) = average_utility_H2(3,:,:) / run_time / 10;
    
    average_utility_H3 = zeros(3,3,6);
    raw_utility_H3 = utility_H3 .* reshape(all_customers(3,:,:,:,:) == 1,[50,3,100,6]);
    average_utility_H3(1,:,:) = sum(raw_utility_H3,[1,3]); % S
    average_utility_H3(1,:,:) = average_utility_H3(1,:,:) / run_time / 10;
    raw_utility_H3 = utility_H3 .* reshape(all_customers(3,:,:,:,:) == 2,[50,3,100,6]);
    average_utility_H3(2,:,:) = sum(raw_utility_H3,[1,3]); % U
    average_utility_H3(2,:,:) = average_utility_H3(2,:,:) / run_time / 10;
    raw_utility_H3 = utility_H3 .* reshape(all_customers(3,:,:,:,:) == 3,[50,3,100,6]);
    average_utility_H3(3,:,:) = sum(raw_utility_H3,[1,3]); % W
    average_utility_H3(3,:,:) = average_utility_H3(3,:,:) / run_time / 80;
    
    average_utility_H4 = zeros(3,3,6);
    raw_utility_H4 = utility_H4 .* reshape(all_customers(4,:,:,:,:) == 4,[50,3,100,6]);
    average_utility_H4(1,:,:) = sum(raw_utility_H4,[1,3]); % C
    average_utility_H4(1,:,:) = average_utility_H4(1,:,:) / run_time / 80;
    raw_utility_H4 = utility_H4 .* reshape(all_customers(4,:,:,:,:) == 2,[50,3,100,6]);
    average_utility_H4(2,:,:) = sum(raw_utility_H4,[1,3]); % U
    average_utility_H4(2,:,:) = average_utility_H4(2,:,:) / run_time / 10;
    raw_utility_H4 = utility_H4 .* reshape(all_customers(4,:,:,:,:) == 3,[50,3,100,6]);
    average_utility_H4(3,:,:) = sum(raw_utility_H4,[1,3]); % W
    average_utility_H4(3,:,:) = average_utility_H4(3,:,:) / run_time / 10;
    
    average_utility_H5 = zeros(3,3,6);
    raw_utility_H5 = utility_H5 .* reshape(all_customers(5,:,:,:,:) == 4,[50,3,100,6]);
    average_utility_H5(1,:,:) = sum(raw_utility_H5,[1,3]); % C
    average_utility_H5(1,:,:) = average_utility_H5(1,:,:) / run_time / 10;
    raw_utility_H5 = utility_H5 .* reshape(all_customers(5,:,:,:,:) == 2,[50,3,100,6]);
    average_utility_H5(2,:,:) = sum(raw_utility_H5,[1,3]); % U
    average_utility_H5(2,:,:) = average_utility_H5(2,:,:) / run_time / 80;
    raw_utility_H5 = utility_H5 .* reshape(all_customers(5,:,:,:,:) == 3,[50,3,100,6]);
    average_utility_H5(3,:,:) = sum(raw_utility_H5,[1,3]); % W
    average_utility_H5(3,:,:) = average_utility_H5(3,:,:) / run_time / 10;
    
    average_utility_H6 = zeros(3,3,6);
    raw_utility_H6 = utility_H6 .* reshape(all_customers(6,:,:,:,:) == 4,[50,3,100,6]);
    average_utility_H6(1,:,:) = sum(raw_utility_H6,[1,3]); % C
    average_utility_H6(1,:,:) = average_utility_H6(1,:,:) / run_time / 10;
    raw_utility_H6 = utility_H6 .* reshape(all_customers(6,:,:,:,:) == 2,[50,3,100,6]);
    average_utility_H6(2,:,:) = sum(raw_utility_H6,[1,3]); % U
    average_utility_H6(2,:,:) = average_utility_H6(2,:,:) / run_time / 10;
    raw_utility_H6 = utility_H6 .* reshape(all_customers(6,:,:,:,:) == 3,[50,3,100,6]);
    average_utility_H6(3,:,:) = sum(raw_utility_H6,[1,3]); % W
    average_utility_H6(3,:,:) = average_utility_H6(3,:,:) / run_time / 80;
    
    % Broadcast Check
    average_utility_H1(:,2,1) = average_utility_H1(:,1,1);
    average_utility_H2(:,2,1) = average_utility_H2(:,1,1);
    average_utility_H3(:,2,1) = average_utility_H3(:,1,1);
    average_utility_H4(:,2,1) = average_utility_H4(:,1,1);
    average_utility_H5(:,2,1) = average_utility_H5(:,1,1);
    average_utility_H6(:,2,1) = average_utility_H6(:,1,1);
    % Titles
    title1 = "Average Utility per Capita: Private Information - Pricing: " + num2str(mod(i-1,3));
    title2 = "Average Utility per Capita: Broadcast - Pricing: " + num2str(mod(i-1,3));
    if i / 4 < 1
        title1 = title1 + " (Graph Nodes)";
        title2 = title2 + " (Graph Nodes)";
    else
        title1 = title1 + " (Neighborhoods)";
        title2 = title2 + " (Neighborhoods)";
    end
    
    % Figures
    figure(6*i+147)
    hold on
    plot(1:1:Hours,reshape(average_utility_H1(1,1,:),[1,6]),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H1(2,1,:),[1,6]),'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H1(3,1,:),[1,6]),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H4(1,1,:),[1,6]),'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H4(2,1,:),[1,6]),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H4(3,1,:),[1,6]),'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('H1: $\Gamma = $S, $\Omega = $P','H1: $\Gamma = $U, $\Omega = $P',...
            'H1: $\Gamma = $W, $\Omega = $P','H4: $\Gamma = $C, $\Omega = $P',...
            'H4: $\Gamma = $U, $\Omega = $P','H4: $\Gamma = $W, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
   
    figure(6*i+148)
    hold on
    plot(1:1:Hours,reshape(average_utility_H1(1,2,:),[1,6]),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H1(2,2,:),[1,6]),'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H1(3,2,:),[1,6]),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H4(1,2,:),[1,6]),'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H4(2,2,:),[1,6]),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H4(3,2,:),[1,6]),'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('H1: $\Gamma = $S, $\Omega = $B','H1: $\Gamma = $U, $\Omega = $B',...
            'H1: $\Gamma = $W, $\Omega = $B','H4: $\Gamma = $C, $\Omega = $B',...
            'H4: $\Gamma = $U, $\Omega = $B','H4: $\Gamma = $W, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(6*i+149)
    hold on
    plot(1:1:Hours,reshape(average_utility_H2(1,1,:),[1,6]),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H2(2,1,:),[1,6]),'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H2(3,1,:),[1,6]),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H5(1,1,:),[1,6]),'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H5(2,1,:),[1,6]),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H5(3,1,:),[1,6]),'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('H2: $\Gamma = $S, $\Omega = $P','H2: $\Gamma = $U, $\Omega = $P',...
            'H2: $\Gamma = $W, $\Omega = $P','H5: $\Gamma = $C, $\Omega = $P',...
            'H5: $\Gamma = $U, $\Omega = $P','H5: $\Gamma = $W, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
   
    figure(6*i+150)
    hold on
    plot(1:1:Hours,reshape(average_utility_H2(1,2,:),[1,6]),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H2(2,2,:),[1,6]),'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H2(3,2,:),[1,6]),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H5(1,2,:),[1,6]),'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H5(2,2,:),[1,6]),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H5(3,2,:),[1,6]),'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('H2: $\Gamma = $S, $\Omega = $B','H2: $\Gamma = $U, $\Omega = $B',...
            'H2: $\Gamma = $W, $\Omega = $B','H5: $\Gamma = $C, $\Omega = $B',...
            'H5: $\Gamma = $U, $\Omega = $B','H5: $\Gamma = $W, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
    
    figure(6*i+151)
    hold on
    plot(1:1:Hours,reshape(average_utility_H3(1,1,:),[1,6]),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H3(2,1,:),[1,6]),'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H3(3,1,:),[1,6]),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H6(1,1,:),[1,6]),'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H6(2,1,:),[1,6]),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H6(3,1,:),[1,6]),'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title1)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('H3: $\Gamma = $S, $\Omega = $P','H3: $\Gamma = $U, $\Omega = $P',...
            'H3: $\Gamma = $W, $\Omega = $P','H6: $\Gamma = $C, $\Omega = $P',...
            'H6: $\Gamma = $U, $\Omega = $P','H6: $\Gamma = $W, $\Omega = $P');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
   
    figure(6*i+152)
    hold on
    plot(1:1:Hours,reshape(average_utility_H3(1,2,:),[1,6]),'s-','Color',ccc(1,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H3(2,2,:),[1,6]),'s-','Color',ccc(2,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H3(3,2,:),[1,6]),'s-','Color',ccc(3,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H6(1,2,:),[1,6]),'s-','Color',ccc(4,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H6(2,2,:),[1,6]),'s-','Color',ccc(5,:),'LineWidth',2)
    plot(1:1:Hours,reshape(average_utility_H6(3,2,:),[1,6]),'s-','Color',ccc(6,:),'LineWidth',2)
    grid on
    title(title2)
    ylabel('Average Utility per Capita')
    xlabel('Hours')
    lll = legend('H3: $\Gamma = $S, $\Omega = $B','H3: $\Gamma = $U, $\Omega = $B',...
            'H3: $\Gamma = $W, $\Omega = $B','H6: $\Gamma = $C, $\Omega = $B',...
            'H6: $\Gamma = $U, $\Omega = $B','H6: $\Gamma = $W, $\Omega = $B');
    set(lll,'Interpreter','Latex');
    set(lll,'FontSize',14);
    hold off
end
toc;

clear
