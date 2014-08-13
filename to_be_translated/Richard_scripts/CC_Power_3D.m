function P=CC_power_3D(u,depths,turbine_top, turbine_bottom,A)


depths_in_turbine=depths>=turbine_bottom & depths<=turbine_top;

nem = 0.0016*u.^4-0.0324*u.^3+ 0.1369*u.^2-0.1534*u + 0.8396;

rho=1025;

TSR=4.3;

Cp = -0.0242*TSR^2 + 0.1963*TSR - 0.0049;

P = Cp*nem.*(0.5*rho*A*u.^3);

end