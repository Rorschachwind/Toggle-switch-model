function dydt = repress_fun(Time,y)
global d kmax1 K2 Ot3 K5 kmax2 K3 k0 d1 
%y(1)=NRI y(2)=LacI
%This function is the ODE for toggle switch model
dydt = [-1.0 * (d * y(1)) + (k0 + (kmax1 * K2 * Ot3 * y(1) ^ 2.0) / ((1.0 + K2 * y(1) ^ 2.0) * (1.0 + K5 * y(2) ^ 2.0)));
        -1.0 * (d1 * y(2)) + (kmax2 * K3 * y(1) ^ 2.0) / (1.0 + K3 * y(1) ^ 2.0);
        ];
end
