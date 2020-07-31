%% This fits a Rayleigh distribution to the measured coalescence times to sample data
% **Replace sample data with the required data set before use ****

% User Data:
SampleNames = {'data/G3.mat','data/2019Vermant.mat', 'data/Sample3630_250PPM.mat','data/1umFilteredOil.mat'}; % Data sets with coalescence times.
K=[1,1,1,2];

%Plot color code
SampleMarkerColor={'b',[1 0.5 0],'g','k','b'};
SampleAlphaMarkerColor={[0.73 .83 0.96],[0.95 0.87 0.73],[0.75,0.95,0.75],[0.6 0.6 0.6],[0.73 .83 0.96]};
SampleMarkerSymbol={'o','<','s','>','p'};
LineStyleList={'-','-','-','-','-.'};


for i = 1:length(SampleNames)

% Load Sample data: 
Data = load(SampleNames{i});


% Compute fractions
if ~isfield(Data,'DrainFracTrials')
Data.DrainFracTrials=(1:length(Data.DrainTimeForTrials))/length(Data.DrainTimeForTrials);
end

h(i) = plot(sort(Data.DrainTimeForTrials),Data.DrainFracTrials,SampleMarkerSymbol{i},'Color',SampleAlphaMarkerColor{i},'MarkerFaceColor',SampleAlphaMarkerColor{i},'Linewidth',1.4);
hold on


%Compute the Rayleigh Distribution:
N = length(Data.DrainTimeForTrials);
mu=mean(Data.DrainTimeForTrials);
FracTime=0:0.01:max(Data.DrainTimeForTrials)+5*mu;

[muRayleigh,ratios]=EMRayleigh(Data.DrainTimeForTrials,K(i),10000);
DrainCDF = 0;
%Handle low mu:
if muRayleigh(1)<0.1
    ratios(2:end) = ratios(2:end) + ratios(1)/(K(i)-1);
    ratios(1)=0;
end

for k=1:K(i)
DrainCDF=DrainCDF+ ratios(k)*(1-exp(-FracTime.^2./(2*muRayleigh(k))));
end 
hNormal(i) = plot(FracTime,DrainCDF,'Color',SampleMarkerColor{i},'LineStyle',LineStyleList{i},'LineWidth',1.4);
end


%Plot legend
fontname = 'Helvetica';
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
set(0,'defaultTextInterpreter','latex');
  

set(gca,'FontName','Helvetica','FontSize',15,'Linewidth',1.1,'xScale','log');
 
Lhandle=legend(h,{['Bubbles', 10, 'Suja et.al (2018)'], ['Antibubbles' ,10, 'Vermant et.al (2019)'], [' Silicone oil drops',10, 'Milad et.al (2020)'],[' Bubbles' ,10, 'Suja et.al (2020)']},'Fontsize',14);
set(Lhandle,'box','off','Position',[0.3226 0.6303 0.1441 0.2703],'FontSize',14);
M = findobj(Lhandle,'type','line');
set(M,'linewidth',1.4) 
  
xlabelhandle= xlabel('Time (s)','FontName','Helvetica','FontSize',19);
ylabelhandle= ylabel('Fraction of bubbles that ruptured before time $\tau$','FontName','Helvetica','FontSize',19);
xlim([0.1 10000]);
ylim([0 1]);
axis square
titlehandle = title('Rayleigh Distribution','FontName','Helvetica','FontSize',17);