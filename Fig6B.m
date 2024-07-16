% BA concentration along intestine, input of 1 with recirculation,
% Comparison between ASBT and ΔASBT (with passive reabsorption)

p1=[0.1 0    0     0    0.0015    0.0050];    % ΔASBT
p2=[0.1 0.11    0     0    0.0015    0.0050]; % ASBT

[int1 int2]=simulation(p2);
b1=bar([1:15]+0.2,[int1;int2],0.7,'stacked','LineWidth',1);
b1(2).FaceColor=[28    117    188]/256;
b1(1).FaceColor=[225    223    35]/256;

hold on
[int1 int2]=simulation(p1);
b2=bar([1:15],[int1;int2],0.7,'stacked','LineWidth',1);
b2(2).FaceColor=0.7*[1 1 1]+0.3*[28    117    188]/256;
b2(1).FaceColor=0.7*[1 1 1]+0.3*[225    223    35]/256;
axis([0 16 0 22])
set(gca,'FontSize',15,'YTick',[0 2 10 20],'XTick',[],'box','off')
xlabel('Length of the intestinal tract','FontSize',20)
ylabel('BA (norm.)','FontSize',20)