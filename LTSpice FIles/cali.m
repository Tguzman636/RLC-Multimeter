g.table = readtable("cali-100kHz-100kOhm-LT1055.txt");
g.v0_t0 = find(g.table.V_v_0_>2,1); % starting index for high V0
g.v0_t1 = find(g.table.V_v_0_(g.v0_t0:end)<2,1)+g.v0_t0-1; % ending index for high V0
g.v0_t2 = find(g.table.V_v_0_(g.v0_t1+1:end)>2,1)+g.v0_t1; % starting index for high V0
g.v90_t0 = find(g.table.V_v_90_>2,1); % starting index for high V0
g.v90_t1 = find(g.table.V_v_90_(g.v90_t0:end)<2,1)+g.v90_t0-1; % ending index for high V0
g.v90_t2 = find(g.table.V_v_90_(g.v90_t1+1:end)>2,1)+g.v90_t1; % starting index for high V0
g.func = @(var,time,t0,t1)(sum(var(t0:t1-1).*(time(t0+1:t1)-time(t0:t1-1))/(time(t1)-time(t0))));


vx = calcz(g,g.table.V_v_x_);
vinvx = calcz(g,g.table.V_v_in_v_x_);
vxvout = calcz(g,g.table.V_v_x_v_out_);

zp=vx.val/(vinvx.val-vxvout.val)*1e5
rt=1/real(1/zp)
ct=imag(1/zp)*1e-5

function z = calcz(g,var)
    z.dc = g.func(var,g.table.time,g.v0_t0,g.v0_t2);
    z.v0 = g.func(var,g.table.time,g.v0_t0,g.v0_t1)/2-z.dc/2;
    z.v90 = g.func(var,g.table.time,g.v90_t0,g.v90_t1)/2-z.dc/2;
    z.amp = sqrt((z.v0^2+z.v90^2))*pi;
    z.rad = atan2(-z.v90,z.v0);
    z.deg = z.rad/pi*180;
    z.val = z.amp*exp(1j*z.rad);
end

