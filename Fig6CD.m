% BA concentrations in WT: SI=215, 88% primary; COL=9, 32% primary
% BA concentrations in KO: SI=44, 60% primary COL=43, 14% primary

% Optimized parameters
p=[0.475    0.005     0.067    0.0012    0.0010    0.0089    1.0000];

% simulations
T=1000000;
opts = odeset('Events',@(t,x) events(t,x,p));

% simulate WT 
p=abs([p 0]);
[tt1,xx1]=ode113(@(t,x) F(t,x,p), [0 T], zeros(1,30),opts);

% simulate KO
p(end)=1;
[tt2,xx2]=ode113(@(t,x) F(t,x,p), [0 T], zeros(1,30),opts);

x1=xx1(end,:);
x2=xx2(end,:);

figure
subplot(2,1,1)
int1=x1(1:15);
int2=x1(16:30);
iint=int1+int2;

hold on
b=bar([1:8],[int1(1:8);int2(1:8)]',1,'stacked');
b(2).FaceColor=[28    117    188]/256;
b(1).FaceColor=[225    223    35]/256;
b=bar([9:10],[int1(9:10);int2(9:10)]',1,'stacked');
b(2).FaceColor=[28    117    188]/256;
b(1).FaceColor=[225    223    35]/256;
b=bar([11:15],[int1(11:15);int2(11:15)]',1,'stacked');
b(2).FaceColor=[28    117    188]/256;
b(1).FaceColor=[225    223    35]/256;
axis([0 16 0 320])
set(gca,'XTick',[],'YTick',[0 100 200 300],'FontSize',15)
ylabel('BA (\mumol)','FontSize',20)
xlabel('Length of the intestinal tract','FontSize',20)
legend({'Primary BA','Secondary BA'},'FontSize',20)


subplot(2,1,2)
int1=x2(1:15);
int2=x2(16:30);
iint=int1+int2;

hold on
b=bar([1:8],[int1(1:8);int2(1:8)]',1,'stacked');
b(2).FaceColor=[28    117    188]/256;
b(1).FaceColor=[225    223    35]/256;
b=bar([9:10],[int1(9:10);int2(9:10)]',1,'stacked');
b(2).FaceColor=[28    117    188]/256;
b(1).FaceColor=[225    223    35]/256;
b=bar([11:15],[int1(11:15);int2(11:15)]',1,'stacked');
b(2).FaceColor=[28    117    188]/256;
b(1).FaceColor=[225    223    35]/256;
axis([0 16 0 320])
set(gca,'XTick',[],'YTick',[0 100 200 300],'FontSize',15)
ylabel('BA (\mumol)','FontSize',20)
xlabel('Length of the intestinal tract','FontSize',20)
legend({'Primary BA','Secondary BA'},'FontSize',20)


% ODEs for simulations
function dx=F(t,x,p)

    % Transit parameters
    trsi1=0.1;  % small intestine
    trsi2=0.05; % ileum
    trco=0.01;  % colon
    
    a1=p(1);
    psi1=p(2);
    pil1=p(3);
    pco1=p(4);
    a2=p(1);
    psi2=p(2);
    pil2=p(3);
    pco2=p(4);
    cwt=p(5);
    cko=p(6);
    inwt=0.09;
    inko=0.43;

    tr=[trsi1*ones(1,8) trsi2*ones(1,2) trco*ones(1,5)];
    if p(end)==0
        conv=[zeros(1,10) cwt*ones(1,5)];
        p1=[psi1*ones(1,8) (pil1+a1)*ones(1,2) pco1*ones(1,5)];
        p2=[psi2*ones(1,8) (pil2+a2)*ones(1,2) pco2*ones(1,5)];
        in=inwt;
    else
        conv=[zeros(1,10) cko*ones(1,5)];
        p1=[psi1*ones(1,8) (pil1)*ones(1,2) pco1*ones(1,5)];
        p2=[psi2*ones(1,8) (pil2)*ones(1,2) pco2*ones(1,5)];
        in=inko;
    end

    dx(1)=in+sum(p1.*x(1:15)')-(tr(1)+p1(1)+conv(1)).*x(1);
    dx(2:15)=tr(1:14).*x(1:14)'-(tr(2:15)+p1(2:15)+conv(2:15)).*x(2:15)';
    dx(16)=sum(p2.*x(16:30)')-(tr(1)+p2(1)).*x(16)+conv(1).*x(1);
    dx(17:30)=tr(1:14).*x(16:29)'-(tr(2:15)+p2(2:15)).*x(17:30)'+conv(2:15).*x(2:15)';
    dx=dx';
end

function [value,isterminal,direction]=events(t,x,p)
    dx=F(t,x,p);
    value=norm(dx./x)-1e-6;
    isterminal=1;
    direction=0;
end