function [x1 x2]=simulation(p)

% PARAMETERS
% p = [in a psi pil pco c]
% in  -> synthesis
% a   -> active absorption in the ileum
% psi -> passive absorption in small intestine
% pil -> passive absorption in the ileum
% pco -> passive absorption in the colon
% c   -> conversion of primary BA to secondary BA

% simulations
T=1000000;
opts = odeset('Events',@(t,x) events(t,x,p));
[tt,xx]=ode113(@(t,x) F(t,x,p), [0 T], zeros(1,30),opts);
x1=xx(end,1:15);
x2=xx(end,16:30);

% ODEs for simulations
function dx=F(t,x,p)

    % Transit parameters
    trsi1=0.1;  % small intestine
    trsi2=0.05; % ileum
    trco=0.01;  % colon
    
    in=p(1);
    a1=p(2);
    psi1=p(3);
    pil1=p(4);
    pco1=p(5);
    a2=p(2);
    psi2=p(3);
    pil2=p(4);
    pco2=p(5);
    c=p(6);

    tr=[trsi1*ones(1,8) trsi2*ones(1,2) trco*ones(1,5)];
    conv=[zeros(1,10) c*ones(1,5)];
    p1=[psi1*ones(1,8) (pil1+a1)*ones(1,2) pco1*ones(1,5)];
    p2=[psi2*ones(1,8) (pil2+a2)*ones(1,2) pco2*ones(1,5)];

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

end