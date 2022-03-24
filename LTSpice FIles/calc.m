rt = @(str)(readtable(str,"PreserveVariableNames",true));

%% Capacitors
% ape=zeros(2,8);
% cs=4e-12;rs=4e10;
% rf = 1e6;cf = 5e-12;
% i=1;
% table=rt("1.76pF-100kHz-1MOhm-5pF.txt");zt=1/(2*pi*1e5*1.76e-12);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("10pF-100kHz-1MOhm-5pF.txt");zt=1/(2*pi*1e5*1e-11);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("100pF-100kHz-1MOhm-5pF.txt");zt=1/(2*pi*1e5*1e-10);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("1nF-100kHz-1MOhm-5pF.txt");zt=1/(2*pi*1e5*1e-9);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("10nF-100kHz-1MOhm-5pF.txt");zt=1/(2*pi*1e5*1e-8);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("100nF-100kHz-1MOhm-5pF.txt");zt=1/(2*pi*1e5*1e-7);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% rf = 1e3;cf = 5e-9;
% table=rt("10pF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1e-11);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("100pF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1e-10);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("1nF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1e-9);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("10nF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1e-8);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("100nF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1e-7);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("1uF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1e-6);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("10uF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1e-5);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("100uF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1e-4);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% table=rt("1.59mF-100kHz-1kOhm-5nF.txt");zt=1/(2*pi*1e5*1.59e-3);
% [ape(1,i),ape(2,i)]=calct(table,rf,cf,1e5,zt,-90,cs,rs);i=i+1;
% ape(:,1:6) %1e-12 ~ 1e-7
% ape(:,7:end) % 1e-11 ~ 1e-3


%% Resistors
ape=zeros(2,13);
cs=4e-12;rs=4e10;
%rf = 1e6;cf = 5e-12;
% table=rt("10MOhm-100kHz-1MOhm-5pF.txt");zt=1e7;
% [ape(1,1),ape(2,1)]=calct(table,rf,cf,zt,0,cs,rs);
% table=rt("1MOhm-100kHz-1MOhm-5pF.txt");zt=1e6;
% [ape(1,2),ape(2,2)]=calct(table,rf,cf,zt,0,cs,rs);
% table=rt("100kOhm-100kHz-1MOhm-5pF.txt");zt=1e5;
% [ape(1,3),ape(2,3)]=calct(table,rf,cf,zt,0,cs,rs);
% table=rt("10kOhm-100kHz-1MOhm-5pF.txt");zt=1e4;
% [ape(1,4),ape(2,4)]=calct(table,rf,cf,zt,0,cs,rs);
rf = 1e3;cf = 5e-9;
i=1;
% table=rt("1kOhm-100kHz-1kOhm-5nF.txt");zt=1e3;
% [ape(1,5),ape(2,5)]=calct(table,rf,cf,zt,0,cs,rs);
% table=rt("100Ohm-100kHz-1kOhm-5nF.txt");zt=1e2;
% [ape(1,6),ape(2,6)]=calct(table,rf,cf,zt,0,cs,rs);
% table=rt("10Ohm-100kHz-1kOhm-5nF.txt");zt=1e1;
% [ape(1,7),ape(2,7)]=calct(table,rf,cf,zt,0,cs,rs);
table=rt("R0-5-3.txt");zt=1e0;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,0,cs,rs);i=i+1;
table=rt("R3-5-3.txt");zt=1e3;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,0,cs,rs);i=i+1;
table=rt("R6-5-3.txt");zt=1e6;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,0,cs,rs);i=i+1;
table=rt("R6-5-6.txt");zt=1e6;rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,0,cs,rs);i=i+1;
table=rt("R7-5-6.txt");zt=1e7;rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,0,cs,rs);i=i+1;

table=rt("C-11-5-6.txt");zt=1/(2*pi*1e5*1e-11);rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,-90,cs,rs);i=i+1;
table=rt("C-10-5-6.txt");zt=1/(2*pi*1e5*1e-10);rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,-90,cs,rs);i=i+1;

table=rt("L2-5-6.txt");zt=2*pi*1e5*1e2;rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
table=rt("L1-5-6.txt");zt=2*pi*1e5*1e1;rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
table=rt("L0-5-6.txt");zt=2*pi*1e5*1e0;rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
table=rt("L-1-5-6.txt");zt=2*pi*1e5*1e-1;rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
table=rt("L0-5-3.txt");zt=2*pi*1e5*1e0;rf=1e3;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
table=rt("L-1-5-3.txt");zt=2*pi*1e5*1e-1;rf=1e3;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
table=rt("L-2-5-3.txt");zt=2*pi*1e5*1e-2;rf=1e3;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
table=rt("L3-5-6.txt");zt=2*pi*1e5*1432;rf=1e6;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
table=rt("L3-5-8.txt");zt=2*pi*1e5*1432;rf=1e8;cf=5e-6/rf;
[ape(1,i),ape(2,i)]=calct(table,rf,cf,zt,90,cs,rs);i=i+1;
ape'

%% Inductors

%zt_actu=2*pi*f*h;

%% Calculation

function [ape,pe] = calct(table,rf,cf,zt_actu,phase_actu,cs,rs)
    g.table=table;

    g.v0_t0 = find(g.table.("V(v_0)")>2,1); % starting index for high V0
    g.v0_t1 = find(g.table.("V(v_0)")(g.v0_t0:end)<2,1)+g.v0_t0-1; % ending index for high V0
    g.v0_t2 = find(g.table.("V(v_0)")(g.v0_t1+1:end)>2,1)+g.v0_t1; % starting index for high V0
    g.v90_t0 = find(g.table.("V(v_90)")>2,1); % starting index for high V0
    g.v90_t1 = find(g.table.("V(v_90)")(g.v90_t0:end)<2,1)+g.v90_t0-1; % ending index for high V0
    g.v90_t2 = find(g.table.("V(v_90)")(g.v90_t1+1:end)>2,1)+g.v90_t1; % starting index for high V0
    g.func = @(var,time,t0,t1)(sum(var(t0:t1-1).*(time(t0+1:t1)-time(t0:t1-1))/(time(t1)-time(t0))));% inner product average
    f=1/(g.table.time(g.v0_t2)-g.table.time(g.v0_t0)); % calculate frequency from square wave, offsetting error from oscillator
    dth=(g.table.time(g.v90_t0)-g.table.time(g.v0_t0))*f*2*pi; % calculate phase shift from square wave, offsetting error from oscillator
    
    %f=1e5;
    %dth=pi/2;
    
    zs=rs/(1+2*pi*1i*f*cs*rs); % calculate parasitic input impedance Zs
    zf=rf/(1+2*pi*1i*f*cf*rf); % calculate feedback impedance Zf
    
    vyvx = calcz(g,g.table.("V(v_y)")-g.table.("V(v_x)"),dth);
    vxvout = calcz(g,g.table.("V(v_x)")-g.table.("V(v_out)"),dth);
    vx = calcz(g,g.table.("V(v_x)"),dth);

    zt.val=1/(vx.val/vyvx.val/zs+vxvout.val/vyvx.val/zf);
    zt.amp=abs(zt.val);
    zt.rad=angle(zt.val);
    zt.deg=zt.rad/pi*180;

    pe=zt.deg-phase_actu;
    ape=(zt.amp/zt_actu-1)*100;
end

function z = calcz(g,var,dth)
    z.dc = g.func(var,g.table.time,g.v0_t0,g.v0_t2);% calculate DC offset
    z.v0 = g.func(var,g.table.time,g.v0_t0,g.v0_t1)/2-z.dc/2; % calculate 0 phase product
    z.v90 = g.func(var,g.table.time,g.v90_t0,g.v90_t1)/2-z.dc/2; % calculate dth phase product
    z.v90 = (z.v90-z.v0*cos(dth))/sin(dth);
    z.amp = sqrt((z.v0^2+z.v90^2))*pi; % calculate amplitude
    z.rad = atan2(-z.v90,z.v0);
    z.deg = z.rad/pi*180;
    z.val = z.amp*exp(1j*z.rad);
end

