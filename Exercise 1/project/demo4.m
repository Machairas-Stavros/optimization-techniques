clc
close all

syms f1(x) f2(x) f3(x)
f1 = (x-1)^3+(x-4)^2*cos(x);
f2 = exp(-2*x)+(x-2)^2;
f3 = x^2*log(0.5*x)+sin(0.2*x)^2;
df1 = diff(f1,x);
df2 = diff(f2,x);
df3 = diff(f3,x);

handlefunctions = {matlabFunction(f1,'Vars',{x}),matlabFunction(f2,'Vars',{x}),matlabFunction(f3,'Vars',{x})};
dfhandlefunctions = {matlabFunction(df1,'Vars',{x}),matlabFunction(df2,'Vars',{x}),matlabFunction(df3,'Vars',{x})};
legendfunctions = {'f_1(x)={(x-1)}^3+{(x-4)}^2*cos(x)','f_2(x)=e^{-2x}+{(x-2)}^2','f_3(x)=x^2*ln(0.5x)+{sin}^2(0.2x)'};

%matlabFunction(f1,'Vars',{x})

search_interval = [0,3];
lamda = 0.002:0.001:0.006;

plot_interval = 0:0.001:3;

figure
hold on
plot(plot_interval,handlefunctions{1}(plot_interval))
plot(plot_interval,handlefunctions{2}(plot_interval))
plot(plot_interval,handlefunctions{3}(plot_interval))
hold off
title('Initial Functions Plot')
ylabel('f(x)')
xlabel('x')
legend(legendfunctions)

results = cell(1,numel(handlefunctions));

for i=1:numel(handlefunctions)
    for j=1:numel(lamda)
        results{i}{end+1} = BisectDiffAlgorithm(search_interval,lamda(j),dfhandlefunctions{i});
    end
end

for i=1:numel(handlefunctions)
    legendLabels = {legendfunctions{i}};
    ak = [];
    bk = [];

    figure()
    subplot(2,2,[1 3])
    hold on
    plot(plot_interval,handlefunctions{i}(plot_interval))
    for j=1:numel(lamda)
        plot(results{i}{j}(end,:),handlefunctions{i}(results{i}{j}(end,:)),'x','MarkerSize',10)
        legendLabels{end+1} = 'lamda = '+string(lamda(j));
        ak(end+1) = results{i}{j}(end,1);
        bk(end+1) = results{i}{j}(end,2);
    end
    hold off
    title('Interval of Minimum Value')
    subtitle('Applied in f_'+string(i)+'(x)')
    ylabel('f_'+string(i)+'(x)')
    xlabel('x')    
    legend(legendLabels)

    subplot(2,2,2)
    hold on
    plot(lamda,ak,'Marker','x','Color','Blue')
    plot(lamda,bk,'Marker','x','Color','Red')
    hold off    
    title('Interval of Minimum Values')
    subtitle('x Values')
    ylabel('x')
    xlabel('lamda')
    legend('a_k','b_k')
    
    subplot(2,2,4)
    hold on
    plot(lamda,handlefunctions{i}(ak),'Marker','x','Color','Blue')
    plot(lamda,handlefunctions{i}(bk),'Marker','x','Color','Red')
    hold off    
    title('Interval of Minimum Values')
    subtitle('f(x) Values')
    ylabel('f_'+string(i)+'(x)')
    xlabel('lamda')
    legend('a_k','b_k')
    sgtitle('Bisect Derivatives Algorithm with Varying lamda - '+string(legendfunctions{i}))
end

for i=1:numel(handlefunctions)
    figure()
    for j=1:numel(lamda)
        nexttile()
        hold on
        plot(1:1:size(results{i}{j},1),results{i}{j}(:,1),'Marker','x','MarkerSize',10,'Color','Blue')
        plot(1:1:size(results{i}{j},1),results{i}{j}(:,2),'Marker','x','MarkerSize',10,'Color','Red')
        hold off
        title('Lower & Upper Limit Behavior in f_'+string(i)+'(x) - lamda = '+string(lamda(j)))
        ylabel('a_k - b_k')
        xlabel('No. Iterations (k)')
        legend('a_k','b_k')
        sgtitle('Bisect Derivatives Algorithm with Varying lamda - '+string(legendfunctions{i}))
    end
end

clear ak bk i j legendfunctions legendLabels plot_interval results_temp