%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transfer Matrix Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------- SOA condition -----------------
T_for               = zeros(1,div_n - 1);
T_back              = zeros(1,div_n - 1);
R               	= zeros(1,div_n - 1);
T_for(1,1)      	= sqrt(1 - Ref);
T_for(1,div_n - 1)	= sqrt(1 - Ref);
T_back(1,1)         = sqrt(1 - Ref);
T_back(1,div_n - 1)	= sqrt(1 - Ref);
R(1,1)              = sqrt(Ref);
R(1,div_n - 1)		= sqrt(Ref);
carrier_density		= zeros(dat_len,div_n);
tmp_index   		= ones(dat_len,div_n);
index               = tmp_index * index_ini;
total_power         = zeros(dat_len,div_n);
alpha               = zeros(dat_len,div_n);
del_g_sig           = zeros(dat_len,div_n);
del_carrier_density = zeros(dat_len,div_n);
P_ase_to_sig        = zeros(dat_len,div_n);
del_index           = zeros(dat_len,div_n);
differential_gain   = zeros(dat_len,div_n);
effective_carrier_lifetime  = zeros(dat_len,div_n);
phase_noise_change_carrier  = zeros(dat_len,div_n);

%--------------- input signal for SOA --------------
P_for_sig           = zeros(dat_len,div_n);
P_back_sig          = zeros(dat_len,div_n);
phi_for_sig         = zeros(dat_len,div_n);
phi_back_sig   		= zeros(dat_len,div_n);
photon_density_sig 	= zeros(dat_len,div_n);
A_sig               = zeros(dat_len,div_n);
B_sig               = zeros(dat_len,div_n);
gain_sig            = zeros(dat_len,div_n);
%alph_sig   		= zeros(dat_len,div_n);
C_sig               = zeros(dat_len,div_n);
D_sig               = zeros(dat_len,div_n);
F_sig               = zeros(dat_len,div_n);
spm_sig             = zeros(dat_len,div_n);
xpm_sig             = zeros(dat_len,div_n);
fwm_sig             = zeros(dat_len,div_n);
P_spm_sig           = zeros(dat_len,div_n);
P_xpm_sig           = zeros(dat_len,div_n);
P_fwm_sig           = zeros(dat_len,div_n);
phi_spm_sig         = zeros(dat_len,div_n);
phi_xpm_sig         = zeros(dat_len,div_n);
phi_fwm_sig         = zeros(dat_len,div_n);
rel_spm_sig         = zeros(dat_len,div_n);
rel_xpm_sig         = zeros(dat_len,div_n);
rel_fwm_sig         = zeros(dat_len,div_n);
central_rel         = zeros(1,div_n);
%phi_sig             = zeros(dat_len,div_n);
g_sig               = zeros(dat_len,div_n);

%---------------pump light 1 for SOA---------------
P_for_pu1   		= zeros(dat_len,div_n);
P_back_pu1  		= zeros(dat_len,div_n);
phi_for_pu1 		= zeros(dat_len,div_n);
phi_back_pu1		= zeros(dat_len,div_n);
photon_density_pu1  = zeros(dat_len,div_n);
A_pu1           	= zeros(dat_len,div_n);
B_pu1           	= zeros(dat_len,div_n);
gain_pu1        	= zeros(dat_len,div_n);
%alph_pu1       	= zeros(dat_len,div_n);
C_pu1           	= zeros(dat_len,div_n);
D_pu1           	= zeros(dat_len,div_n);
F_pu1               = zeros(dat_len,div_n);
spm_pu1             = zeros(dat_len,div_n);
xpm_pu1             = zeros(dat_len,div_n);
fwm_pu1             = zeros(dat_len,div_n);
g_pu1               = zeros(dat_len,div_n);

%---------------pump light 2 for SOA---------------
P_for_pu2   		= zeros(dat_len,div_n);
P_back_pu2  		= zeros(dat_len,div_n);
phi_for_pu2 		= zeros(dat_len,div_n);
phi_back_pu2		= zeros(dat_len,div_n);
photon_density_pu2  = zeros(dat_len,div_n);
A_pu2           	= zeros(dat_len,div_n);
B_pu2           	= zeros(dat_len,div_n);
gain_pu2        	= zeros(dat_len,div_n);
%alph_pu2       	= zeros(dat_len,div_n);
C_pu2           	= zeros(dat_len,div_n);
D_pu2           	= zeros(dat_len,div_n);
F_pu2               = zeros(dat_len,div_n);
spm_pu2             = zeros(dat_len,div_n);
xpm_pu2             = zeros(dat_len,div_n);
fwm_pu2             = zeros(dat_len,div_n);
g_pu2               = zeros(dat_len,div_n);

%-------------------ase in SOA--------------------
P_for_ase           = zeros(dat_len,div_n,30);
P_back_ase          = zeros(dat_len,div_n,30);

photon_density_ase  = zeros(dat_len,div_n);
g_ase               = zeros(dat_len,div_n);
gain_ase            = zeros(dat_len,div_n,30);
n_sp                = zeros(dat_len,div_n);
E_ase               = zeros(30,1);
tmp_S_ase_div       = zeros(30,1);
%P_for_ase_div       = zeros(30,1);
P_back_ase_div      = zeros(30,1);
G_ase               = zeros(30,1);
tmp_s_ase           = zeros(30,1);

P_ase_to_sig        = zeros(dat_len,div_n);

t              	= zeros(dat_len - 1,1);

for i = 1:dat_len-1
	P_for_sig(i,1)          = input_sig(i,1);
	P_for_pu1(i,1)      	= input_pu1(i,1);
	P_for_pu2(i,1)          = input_pu2(i,1);    
    P_for_ase(i,1)          = 0;
    for ii = 2:div_n -1
		carrier_density(i,ii)	= carrier_den_tr;
    end
end

for j = 2:dat_len - 1
	t(j,1) = del_t * j;
	
	for m = 2:div_n -1 
		tmp_N           = carrier_density(j,m);

		% peak-gain
		g1              = a_g1*tmp_N+b_g1;
		g2              = a_g2*tmp_N+b_g2;
		g_peak          = 0.5*(g1+g2) - sqrt(((g1-g2)/2)^2+cg);

		%polynomial coefficient
		p1              = a_p1*tmp_N+b_p1;
		p2              = a_p2*tmp_N+b_p2;
		p               = 0.5*(p1+p2) - sqrt(((p1-p2)/2)^2+cp);

		%peak-gain wavelength
		lam_peak1       = a_l1*tmp_N+b_l1;
		lam_peak2       = a_l2*tmp_N+b_l2;
		lam_peak        = 0.5*(lam_peak1+lam_peak2) + sqrt(((lam_peak1-lam_peak2)/2)^2+cl);

		%energy of spontaneous emission light
        for v = 1 : 30
            E_ase(v,1)          = (h * c) / (1499e-9+v*4e-9);
        end
        
		%photon density
		tmp_S_sig       = (abs(P_for_sig(j,m)) + abs(P_back_sig(j,m))) / (group_velocity * A_cross * E_sig);
		tmp_S_pu1       = (abs(P_for_pu1(j,m)) + abs(P_back_pu1(j,m))) / (group_velocity * A_cross * E_pu1);
		tmp_S_pu2       = (abs(P_for_pu2(j,m)) + abs(P_back_pu2(j,m))) / (group_velocity * A_cross * E_pu2);
        
        tmp_S_ase = 0;
        
        for v=1:30
            tmp_S_ase_div(v,1)  = (abs(P_for_ase(j,m,v)) + abs(P_back_ase(j,m,v))) / (group_velocity * A_cross * E_ase(v,1));
            tmp_S_ase = tmp_S_ase + tmp_S_ase_div(v,1);
        end
        
		photon_density_sig(j,m)     = tmp_S_sig;
		photon_density_pu1(j,m)     = tmp_S_pu1;
		photon_density_pu2(j,m)     = tmp_S_pu2;
        photon_density_ase(j,m)     = tmp_S_ase;
        
		%material gain for signal an P_for_sig & pump light
		g_sig(j,m)           = (g_peak - p*(lam_sig-lam_peak)^2 + gain_const3*(lam_sig-lam_peak)^3)/(1 + gain_compress * (tmp_S_sig + tmp_S_pu1 + tmp_S_pu2));
		g_pu1(j,m)           = (g_peak - p*(lam_pu1-lam_peak)^2 + gain_const3*(lam_pu1-lam_peak)^3)/(1 + gain_compress * (tmp_S_sig + tmp_S_pu1 + tmp_S_pu2));
		g_pu2(j,m)           = (g_peak - p*(lam_pu2-lam_peak)^2 + gain_const3*(lam_pu2-lam_peak)^3)/(1 + gain_compress * (tmp_S_sig + tmp_S_pu1 + tmp_S_pu2));

        del_g_sig(j,m)       = g_sig(j,m) - g_sig(j-1,m);
        del_g_pu1       = g_pu1(j,m) - g_pu1(j-1,m);
        del_g_pu2       = g_pu2(j,m) - g_pu2(j-1,m);
		
        R_ase=0;
        
        for v=1:30
            g_ase(v,1)  =  (g_peak - p*((1499e-9+4e-9*v)-lam_peak)^2 +...
                gain_const3*((1499e-9+4e-9*v)-lam_peak)^3)/(1 + gain_compress * (tmp_S_sig + tmp_S_pu1 + tmp_S_pu2));
            %ASE reconbination rate
            R_ase   = R_ase + gam * g_ase(v,1) * (P_for_ase(j,m,v) + P_back_ase(j,m,v)) / E_ase(v,1) / A_cross;
        end
        
        %solving rate equation
		carrier_density(j+1, m) = del_t * (I_current / (q_unit * volume) - tmp_N * (recom_rate1 + recom_rate2 * tmp_N + recom_rate3 * (tmp_N)^2) - R_ase...
					 - group_velocity * gam * (g_sig(j,m) * tmp_S_sig + g_pu1(j,m) * tmp_S_pu1 + g_pu2(j,m) * tmp_S_pu2)) + carrier_density(j,m);

        del_carrier_density(j,m) = carrier_density(j+1,m)- carrier_density(j,m);
        
        %active layer loss
        loss_in_active  = k_0 + k_1 * carrier_density(j,m);
        
		%modal gain
		G_sig           = gam * (g_sig(j,m) - loss_in_active) - (1 - gam) * loss_in_clad - scatt_loss;
		G_pu1           = gam * (g_pu1(j,m) - loss_in_active) - (1 - gam) * loss_in_clad - scatt_loss;
		G_pu2           = gam * (g_pu2(j,m) - loss_in_active) - (1 - gam) * loss_in_clad - scatt_loss;
		gain_sig(j,m)   = G_sig;
		gain_pu1(j,m)   = G_pu1;
		gain_pu2(j,m)   = G_pu2;
        
        for v=1:30
        	G_ase(v,1)      = abs(gam * (g_ase(v,1) - loss_in_active) - (1 - gam) * loss_in_clad - scatt_loss);
            gain_ase(j,m,v)   = G_ase(v,1);
        end
        
        %spontaneous emission foctor
        n_sp(j,m)       = carrier_density(j,m)/(carrier_density(j,m)-(carrier_den_tr-1e23)); %preventing a divergency 
        
        for v = 1 : 30
            tmp_s_ase(v,1)   = E_ase(v,1) * n_sp(j,m) * del_nyu * (exp((G_ase(v,1))*del_z)-1);
        end

        P_for_sig(j,m)      = P_for_sig(j,m) * exp(G_sig * del_z);
		P_back_sig(j,m)     = P_back_sig(j,m) * exp(G_sig * del_z);
		P_for_pu1(j,m)      = P_for_pu1(j,m) * exp(G_pu1 * del_z);
		P_back_pu1(j,m)     = P_back_pu1(j,m) * exp(G_pu1 * del_z);
		P_for_pu2(j,m)      = P_for_pu2(j,m) * exp(G_pu2 * del_z);
		P_back_pu2(j,m)     = P_back_pu2(j,m) * exp(G_pu2 * del_z);
        
        for v=1:30
            P_for_ase(j,m,v)      = P_for_ase(j,m,v) * exp(G_ase(v,1) * del_z) + tmp_s_ase(v,1);
            P_back_ase(j,m,v)     = P_back_ase(j,m,v) * exp(G_ase(v,1) * del_z) + tmp_s_ase(v,1);
        end
        
        P_ase_to_sig(j,m)       = P_for_ase(j,m,13);
        P_ase_to_sig(j,m)        = P_ase_to_sig(j,m) .* exp(Im * randn());
		
        %index
		index(j,m)  = index_ini + diff_ref_index * carrier_density(j+1,m);
        del_index(j,m) = index(j,m)-index(j-1,m);
        
        %alpha parameter
        alpha(j,m) = -4 * pi() / lam_sig * del_index(j,m) / del_g_sig(j,m);      %usually alpha ÅÅ from 4 to 6
        
        %differential gain
        differential_gain(j,m) = del_g_sig(j,m) / del_carrier_density(j,m);            %as a reference dif. gain = -4.212*10^-4
        
        %different carrier lifetime
        tau_s = 1/(recom_rate1+recom_rate2*tmp_N+recom_rate3*tmp_N^2);
        
        %effective carrier lifetime
        effective_carrier_lifetime(j,m) = 1/(1/tau_s + group_velocity*((g_sig(j,m)*del_z+1)/g_sig(j,m)/del_z*del_g_sig(j,m)/del_carrier_density(j,m)*photon_density_sig(j,m)...
            +(g_pu1(j,m)*del_z+1)/g_pu1(j,m)/del_z*del_g_pu1/del_carrier_density(j,m)*photon_density_pu1(j,m)...
            +(g_pu2(j,m)*del_z+1)/g_pu2(j,m)/del_z*del_g_pu2/del_carrier_density(j,m)*photon_density_pu2(j,m)));
        
        %total input lights power
        total_power(j,m) = abs(P_for_sig(j,m)) + abs(P_for_pu1(j,m)) + abs(P_for_pu2(j,m));
        
        phase_noise_change_carrier(j,m) = alpha(j,m) * gam * group_velocity * differential_gain(j,m) * del_carrier_density(j,m) / 2 * del_t; %íPà ìIÇ…Ç±ÇøÇÁ
	end

	%reflectivity and transmittance
	for k = 2:div_n -2
		R(1,k) 		= (index(j,k+1) - index(j,k)) / (index(j,k+1) + index(j,k));
		T_for(1,k) 	= (2 * index(j,k)) / (index(j,k+1) + index(j,k));
		T_back(1,k)	= (2 * index(j,k+1)) / (index(j,k+1) + index(j,k));
    end

	for s = 1:div_n
        D_sig(j,s) = abs(sqrt(P_for_sig(j,s)));
		D_pu1(j,s) = abs(sqrt(P_for_pu1(j,s)));
  		D_pu2(j,s) = abs(sqrt(P_for_pu2(j,s)));
        C_sig(j,s) = angle(P_for_sig(j,s));
        C_pu1(j,s) = angle(P_for_pu1(j,s));
        C_pu2(j,s) = angle(P_for_pu2(j,s));
    end
    
    for o = 1 : div_n -1
        A_sig(j,o) = D_sig(j,o) * exp(Im * (C_sig(j,o) + phase_noise_change_carrier(j,o)));
        A_pu1(j,o) = D_pu1(j,o) * exp(Im * C_pu1(j,o));
        A_pu2(j,o) = D_pu2(j,o) * exp(Im * C_pu2(j,o));
    end

	%affecting the nonlinearly fenomenon
	for x = 1:div_n-1
        spm_sig(j,x)    = 1*exp(Im*gamma*(abs(A_sig(j,x)))^2*A_sig(j,x));
        spm_pu1(j,x)    = 1*exp(Im*gamma*(abs(A_pu1(j,x)))^2*A_pu1(j,x));
        spm_pu2(j,x)    = 1*exp(Im*gamma*(abs(A_pu2(j,x)))^2*A_pu2(j,x));
        
        xpm_sig(j,x)    = 1*exp(Im*gamma*2*((abs(A_pu1(j,x)))^2 + (abs(A_pu2(j,x)))^2)*A_sig(j,x));
        xpm_pu1(j,x)    = 1*exp(Im*gamma*2*((abs(A_sig(j,x)))^2 + (abs(A_pu2(j,x)))^2)*A_pu1(j,x));
        xpm_pu2(j,x)    = 1*exp(Im*gamma*2*((abs(A_pu1(j,x)))^2 + (abs(A_sig(j,x)))^2)*A_pu2(j,x));
             
        fwm_sig(j,x)    = gamma*2*(A_pu1(j,x))*(A_pu2(j,x))*conj(A_sig(j,x))*exp((-1)*Im*beta);
        fwm_pu1(j,x)    = gamma*(A_sig(j,x))^2*conj(A_pu2(j,x))*exp(Im*beta);
        fwm_pu2(j,x)    = gamma*(A_sig(j,x))^2*conj(A_pu1(j,x))*exp(Im*beta);
        
        F_sig(j,x) = A_sig(j,x);
        F_pu1(j,x) = A_pu1(j,x);
        F_pu2(j,x) = A_pu2(j,x);
        
        A_sig(j,x) = (A_sig(j,x) + fwm_sig(j,x)) * spm_sig(j,x);
        A_pu1(j,x) = (A_pu1(j,x) + fwm_pu1(j,x)) * spm_pu1(j,x) * xpm_pu1(j,x);
        A_pu2(j,x) = (A_pu2(j,x) + fwm_pu2(j,x)) * spm_pu2(j,x) * xpm_pu2(j,x);
    end
   
	%transfer the light to next block 
	for u = 1:div_n -1

        A_sig(j+1,u+1)          = T_for(1,u) * ( A_sig(j,u)) - R(1,u) * ( B_sig(j,u+1));
		B_sig(j+1,u)            = T_back(1,u) * ( B_sig(j,u+1)) + R(1,u) * ( A_sig(j,u));
		A_pu1(j+1,u+1)          = T_for(1,u) * ( A_pu1(j,u)) - R(1,u) * ( B_pu1(j,u+1));
		B_pu1(j+1,u)            = T_back(1,u) * ( B_pu1(j,u+1) ) + R(1,u) * ( A_pu1(j,u));
		A_pu2(j+1,u+1)          = T_for(1,u) * ( A_pu2(j,u)) - R(1,u) * ( B_sig(j,u+1));
		B_pu2(j+1,u)            = T_back(1,u) * ( B_sig(j,u+1)) + R(1,u) * ( A_pu2(j,u));

        for v=1:30
            P_for_ase(j+1,u+1,v)  	= (T_for(1,u))^2 * P_for_ase(j,u,v) + (R(1,u))^2 * P_back_ase(j,u+1,v);
            P_back_ase(j+1,u,v)  	= (T_back(1,u))^2 * P_back_ase(j,u+1,v) + (R(1,u))^2 * P_for_ase(j,u,v);
        end
    end
	 
	%obtaining the optical power and phase from complex amplitude
	for w = 1:div_n -1
        P_for_sig(j+1,w+1)      = (abs(A_sig(j+1,w+1)))^2 * exp(angle(A_sig(j+1,w+1)) * Im);
        P_back_sig(j+1,w)       = ((B_sig(j+1,w)))^2;
		P_for_pu1(j+1,w+1)      = (abs(A_pu1(j+1,w+1)))^2 * exp(angle(A_pu1(j+1,w+1)) * Im);
		P_back_pu1(j+1,w)       = ((B_pu1(j+1,w)))^2;
		P_for_pu2(j+1,w+1)      = (abs(A_pu2(j+1,w+1)))^2 * exp(angle(A_pu2(j+1,w+1)) * Im);
		P_back_pu2(j+1,w)       = ((B_pu2(j+1,w)))^2;
        P_for_sig(j+1,w+1)      = P_for_sig(j+1,w+1) + P_ase_to_sig(j,w);
    end
end

t(dat_len,1)= del_t * dat_len;