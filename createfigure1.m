function createfigure1(YMatrix1, YMatrix2, Y1)

%CREATEFIGURE1(YMATRIX1,YMATRIX2,Y1)
%  YMATRIX1:  matrix of y data
%  YMATRIX2:  matrix of y data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 05-Jun-2011 14:32:30

% Create figure
figure1 = figure(1); clf;

xdata=[0:30:1600];

months1=['Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';'Jan';'Feb'];
months=[months1;months1;months1];
% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',months,...
    'XTick',xdata,...
    'Position',[0.13 0.6409 0.8037 0.2891],...
    'FontSize',11);
box('on');
hold on;

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix1);
set(plot1(1),'DisplayName','Egg','Color',[0 0 1]);
set(plot1(2),'DisplayName','Larva','Color',[0 1 0]);
set(plot1(3),'DisplayName','Pupa','Color',[1 0 0]);
set(plot1(4),'DisplayName','Nurse Bee','Color',[1 0 1]);
set(plot1(5),'DisplayName','House Bee','Color',[0 1 1]);
set(plot1(6),'DisplayName','Forager','Color',[0 0 0]);

% Create ylabel
ylabel('Number of Bees');
ylim([0,2.2e4]);

% Create title
title('A Healthy Honey Bee Colony Population Dynamics --four-year-simulation');

% Create legend
legend(axes1,'show');


% Create axes
axes2 = axes('Parent',figure1,...
    'XTickLabel',months,...
    'XTick',xdata,...
    'Position',[0.1294 0.302 0.8056 0.2597]);
box('on');
hold on;

% Create multiple lines using matrix input to plot
plot2 = plot(YMatrix2);
set(plot2(1),'DisplayName','Pollen','Marker','o','LineStyle','none',...
    'Color',[0 1 0]);
set(plot2(2),'DisplayName','Honey','Marker','o','LineStyle','none',...
    'Color',[0 0 1]);

% Create ylabel
ylabel('Number of Cells');
ylim([0,2e5])

% Create legend
legend(axes2,'show');

% Create axes
axes3 = axes('Parent',figure1,...
    'XTickLabel',months,...
    'XTick',xdata,...
    'Position',[0.1294 0.03867 0.805 0.2081]);
box('on');
%hold('all');

% Create plot for weight

Weight = (.047*YMatrix1(:, 2)+.158*YMatrix1(:, 3) + .133*(YMatrix1(:, 4)+YMatrix1(:, 5)+YMatrix1(:, 6)) + .23*YMatrix2(:, 1)+ .5*YMatrix2(:, 2))/1000;
disp(Weight)
plot(Weight,'Parent',axes3,'Marker','o','Color',[0 0 0]);

% Create ylabel
ylabel('Weight (kg)');

print -depsc output.eps

end 
