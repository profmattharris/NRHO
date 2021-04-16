function [dv2b,dv3b,tflight] = TwoImpOptimizer(rev,flag)

%--------------------------------------------------------------------------
% The function has inputs of NRHO rev fraction and guess type. It computes
% locally optimal transfers from LLO to NRHO. It outputs a two-body delta-v
% estimate, the three-body locally optimal delta-v, and flight time.

% Inputs:
%   1) rev (NRHO rev fraction between 1 and 2)
%   2) flag (1 for Hohmann guess, 2 for two-body optimal guess)
%
% Outputs:
%   1) dv2b (two-body delta-v estimate)
%   2) dv3b (three-body locally optimal delta-v)
%   3) tflight (flight time from LLO to NRHO)
%
%--------------------------------------------------------------------------

%------%
% DATA %
%------%
Rm  = 1736;             % Lunar radius (km)
RM  = 384400;           % Distance to moon from earth (km)
muM = 4902.800076;      % Moon gravitational parameter (km,s)
muE = 398600.436233;    % Earth gravitational parameter (km,s)
f   = 1/sqrt(RM^3/muE); % Moon orbital frequency
tf  = rev*6.993*24*3600; % NRHO flight time past periapse (s)
options = odeset('AbsTol',1e-12,'RelTol',1e-12); % Integration tolerances

% Initial Moon states in ECI...
sm0     = [RM;0;0;0;sqrt(muE/RM);0];

% Initial Gateway states in LCI...
rg0     = 1.0e+05 * [3.841440114932523; 0; 0.045503057686036];
vg0     = [0; 2.452499891637010; 0];
sg0     = [rg0; vg0] - sm0;

% Calculate the final Gateway state in LCI...
[~,s] = ode45(@three_body,[0,tf],sg0,options,muE,muM,RM,f,1);
rgf   = s(end,1:3).'; 	% Final Gateway position in LCI
vgf   = s(end,4:6).'; 	% Final Gateway velocity in LCI

% Needed transfer orbit plane...
khat = cross(rgf,vgf)  / norm( cross(rgf,vgf)  );

% Calculate guess for departure from circular orbit...
rdep = -(Rm+100)*rgf/norm(rgf);
vdep = null([khat.'; rdep.'])*sqrt(muM/(Rm+100));

% Hohmann transfer guess...
if flag == 1
rp = norm(rdep);
ra = norm(rgf);
a  = 1/2*(rp+ra);
e  = 1 - rp/a;
tp = 2*pi*sqrt(a^3/muM);

vdir = sign(1.65-rev)*null([rdep';khat']);
Vc   = sqrt( muM/rp )*vdir/norm(vdir);
vdep = sqrt( 2*(-muM/(2*a) + muM/rp) )*vdir/norm(vdir);
varr = sqrt( 2*(-muM/(2*a) + muM/ra) )*vdir/norm(vdir);
dv2b = norm(vdep-Vc) + norm(vgf-varr);
tflight = tp/2;
end

 
% Two-body optimization guess...
if flag == 2
dth = 160*pi/180;
ops = optimoptions(@fmincon,'Display','iter','MaxFunEvals',5000);
[unk,dv2b] = fmincon(@obj,[rdep/norm(rdep);vdep;vdep;dth],[],[],[],[],[],[],...
                  @con,ops,muM,Rm,rgf,vgf);
%
rdep = unk(1:3)*(Rm+100);
vdep = unk(4:6);
Vc = unk(7:9);
dth = unk(10);
tflight = time_calc(muM,rdep,vdep,dth);
end

% Now back up and find the Gateway position at the departure time.
if (tf-tflight) <= 0
    disp('TIME ERROR')
    return;
else
    [~,s] = ode45(@three_body,[0,tf-tflight],sg0,options,muE,muM,RM,f,1);
end

rg0     = s(end,1:3).';       % Initial Gateway position in LCI
vg0     = s(end,4:6).';       % Initial Gateway velocity in LCI

tgrid = linspace(tf-tflight,tf);
[~,sg] = ode45(@three_body,tgrid,[rg0;vg0],options,muE,muM,RM,f,1);
[~,sv] = ode45(@three_body,tgrid,[rdep;vdep],options,muE,muM,RM,f,1);
[~,s2] = ode45(@two_body,tgrid,[rdep;vdep],options,muM);

figure, plot3(sg(:,1),sg(:,2),sg(:,3),'linewidth',2), grid on
hold on, plot3(sv(:,1),sv(:,2),sv(:,3),'linewidth',2)
hold on, plot3(s2(:,1),s2(:,2),s2(:,3),'--','linewidth',2)
xlabel('X (km)'), ylabel('Y (km)'), zlabel('Z (km)')
legend('Gateway (3BP)','Vehicle (3BP)','Vehicle (2BP)')
title('Two-body Initial Guess in LCI')


%---------------------%
% THREE-BODY SHOOTING %
%---------------------%
% ops = optimoptions('fsolve','Display','iter','MaxFunEvals',5e3);
% pshooting = @(v0)shooting(v0,muE,muM,RM,f,rdep,rgf,tf,tflight,1);
% vvd = fsolve(pshooting,vdep,ops); % Velocity of the vehicle at departure in ECI
% 
% [~,sv] = ode45(@three_body,tgrid,[rdep;vvd],options,muE,muM,RM,f,1);
% figure, plot3(sg(:,1),sg(:,2),sg(:,3),'linewidth',2), grid on
% hold on, plot3(sv(:,1),sv(:,2),sv(:,3),'linewidth',2)
% legend('Gateway (3BP)','Vehicle (3BP)')
% title('3BP 2IMP Shooting in LCI')
% 
% ops = optimoptions(@fmincon,'Display','iter','MaxFunEvals',5000);
% v0 = fmincon(@obj2,vvd,[],[],[],[],[],[],...
%              @con2,ops,rdep,vvd,muM);
% dv1 = vvd - v0;
% dv2 = sg(end,4:6).'-sv(end,4:6).';
% dv3b = norm(dv1)+norm(dv2);


%-------------------------%
% THREE-BODY OPTIMIZATION %
%-------------------------%
ops = optimoptions(@fmincon,'Display','iter','MaxFunEvals',5000,...
                            'OptimalityTolerance',1e-3);
unk = [Vc;rdep/(Rm+100);vdep;tflight];
for i = 0:.05:1
[unk,dv3b] = fmincon(@obj3,unk,[],[],[],[],[],[],...
                     @con3,ops,muE,muM,RM,f,Rm,rgf,vgf,tf,i);
end
rdep = unk(4:6)*(Rm+100);
vdep = unk(7:9);
tflight = unk(10);

[~,sv] = ode45(@three_body,[tf-tflight,tf],[rdep;vdep],options,muE,muM,RM,f,1);
figure, plot3(sg(:,1),sg(:,2),sg(:,3),'linewidth',2), grid on
hold on, plot3(sv(:,1),sv(:,2),sv(:,3),'linewidth',2)
xlabel('X (km)'), ylabel('Y (km)'), zlabel('Z (km)')
legend('Gateway (3BP)','Vehicle (3BP)')
title('Three-body Optimized Solution in LCI')

end % the main function

%------------------------%
% THREE-BODY INTEGRATION %
%------------------------%
function xdot = three_body(t,x,muE,muM,RM,f,scale)
rm = RM*[cos(f*t);sin(f*t);0*t];
vm = RM*f*[-sin(f*t);cos(f*t);0*t];
del = x(1:3);    
nu = x(4:6);
r13 = -rm;
r23 = -(rm+del);

deldot = nu;
nudot  = -muM*del/norm(del)^3 - scale*muE*(r13/norm(r13)^3 - r23/norm(r23)^3);
xdot = [deldot; nudot];
end

%----------------------%
% TWO-BODY INTEGRATION %
%----------------------%
function xdot = two_body(~,x,mu)
r = x(1:3);
v = x(4:6);
rdot = v;
vdot = -mu*r/norm(r)^3;
xdot = [rdot; vdot];
end

%------------------------------%
% THREE-BODY SHOOTING FUNCTION %
%------------------------------%
function F = shooting(v0,muE,muM,RM,f,r0,rgf,tf,tflight,scale)
options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[~,s] = ode45(@three_body,[tf-tflight,tf],[r0;v0],options,muE,muM,RM,f,scale);
F = s(end,1:3).'-rgf;
end

%-----------------------------%
% TWO-BODY F & G OPTIMIZATION %
%-----------------------------%
function J = obj(x,mu,Rm,Rf,Vf)

R0 = x(1:3)*(Rm+100);
V0 = x(4:6);
Vc = x(7:9);
dth  = x(10);
H  = cross(R0,V0);
h  = norm( H );
r0 = norm(R0);
vr0 = dot(V0,R0/r0);
rf = h^2/mu*1/( 1 + (h^2/(mu*r0)-1)*cos(dth) - h*vr0/mu*sin(dth) );

f  = 1 - mu*rf/h^2*(1-cos(dth));
g  = rf*r0/h*sin(dth);
fd = mu/h*(1-cos(dth))/sin(dth)*(mu/h^2*(1-cos(dth))-1/r0-1/rf);
gd = 1-mu*r0/h^2*(1-cos(dth));

R2 = f*R0  + g*V0;
V2 = fd*R0 + gd*V0;

J = norm(V0-Vc)+norm(Vf-V2);
end

function [cin,ceq] = con(x,mu,Rm,Rf,Vf)

R0 = x(1:3)*(Rm+100);
V0 = x(4:6);
Vc = x(7:9);
dth  = x(10);
H  = cross(R0,V0);
h  = norm( H );
r0 = norm(R0);
vr0 = dot(V0,R0/r0);
rf = h^2/mu*1/( 1 + (h^2/(mu*r0)-1)*cos(dth) - h*vr0/mu*sin(dth) );

f  = 1 - mu*rf/h^2*(1-cos(dth));
g  = rf*r0/h*sin(dth);
fd = mu/h*(1-cos(dth))/sin(dth)*(mu/h^2*(1-cos(dth))-1/r0-1/rf);
gd = 1-mu*r0/h^2*(1-cos(dth));

R2 = f*R0  + g*V0;
V2 = fd*R0 + gd*V0;

cin = [];
eq1 = norm(x(1:3)) - 1;
eq2 = norm(Vc) - sqrt(mu/(Rm+100));
eq3 = dot(R0,Vc);
eq4 = (R2-Rf)/rf;
eq5 = dth - 68*pi/180;
ceq = [eq1;eq2;eq3;eq4;0*eq5];

end

function tflight = time_calc(muM,R0,V0,dth)

coe0 = coe_from_sv(R0,V0,muM);
e = coe0(2);
ta0 = coe0(6);
sma = coe0(7); T = 2*pi*sqrt(sma^3/muM);
E0 = 2*atan( sqrt( (1-e)/(1+e) ) * tan(ta0/2) );
M0 = E0 - e*sin(E0);
t0 = T*M0/(2*pi);
taf = ta0+dth;
Ef = 2*atan( sqrt( (1-e)/(1+e) ) * tan(taf/2) );
if Ef < 0
    Ef = 2*pi+Ef;
end
Mf = Ef - e*sin(Ef);
tf = T*Mf/(2*pi);
tflight = (tf-t0);

end


%--------------------------------%
% CIRCULAR VELOCITY OPTIMIZATION %
%--------------------------------%
function J = obj2(v0,rdep,vdep,muM)

J = norm(vdep-v0);

end

function [cin,ceq] = con2(v0,rdep,vdep,muM)

cin = [];
eq1 = norm(v0) - sqrt(muM/norm(rdep));
eq2 = dot(v0,rdep);
ceq = [eq1,eq2];

end


%-------------------------%
% THREE-BODY OPTIMIZATION %
%-------------------------%
function J = obj3(unk,muE,muM,RM,f,Rm,rgf,vgf,tf,scale)
v0   = unk(1:3);
rdep = unk(4:6)*(Rm+100);
vdep = unk(7:9);
tflight = unk(10);

options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[~,s] = ode45(@three_body,[tf-tflight,tf],[rdep;vdep],options,muE,muM,RM,f,scale);

varr = s(end,4:6).';
J = norm(vdep-v0) + norm(vgf-varr);

end

function [cin,ceq] = con3(unk,muE,muM,RM,f,Rm,rgf,vgf,tf,scale)
v0   = unk(1:3);
rdep = unk(4:6)*(Rm+100);
vdep = unk(7:9);
tflight = unk(10);

options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[~,s] = ode45(@three_body,[tf-tflight,tf],[rdep;vdep],options,muE,muM,RM,f,scale);

eq1 = norm(v0) - sqrt(muM/(Rm+100));
eq2 = norm(rdep)/(Rm+100) - 1;
eq3 = dot(rdep/(Rm+100),v0);
eq4 = ( s(end,1:3).'-rgf )/(Rm+100);
cin = [];
ceq = [eq1; eq2; eq3; eq4];
end