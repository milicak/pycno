clear all
% Ehlert and Leverman 2014;
% Mechanism for potential strengthening of Atlantic overturning prior to collapse
% temperature evolution equations
% for git hub
tempevl = true;
FWforce = false;
AADW = 1; % turn on AADW mode

%diffusivity = 'N2 dependent'; % 'constant diff' ; 'constant energy' ; 'N2 dependent'
diffusivity = 'constant diff'; % 'constant diff' ; 'constant energy' ; 'N2 dependent'
%diffusivity = 'constant energy'; % 'constant diff' ; 'constant energy' ; 'N2 dependent'
% Table 1
C1 = 0.1; %nondimensional coeff for MOC transport
m3s2Sv = 1e-6; 
grav = 9.81; %m2/s

Htopo = 5000; %in m ; Average depth of Atlantic Ocean basin
Bwidth = 1e7; % in m ; Average width of Atlantic Ocean
LN = 3.34*1e6; %in m ; Meridional extend of the northen box
LU = 8.9*1e6;  %in m ; Meridional extend of the tropical box
LS = 3.34*1e6; %in m ; Meridional extend of the southern box
LyN = 1.5*1e6; % Meridional extend of the northen outcropping

rho0 = 1027; %kgm-3
S0 = 35; % psu ; Average salinity of Atlantic Ocean
alpha_T = 2.1e-4; % 1/C thermal coefficient
alphaT = alpha_T*rho0;
beta_S = 8e-4; % 1/psu haline coefficient
betaS = beta_S*rho0;
f0 = 7.5e-5; %1/s
fbeta = 2e-11; % 1/s

kappa_GM = 1000; %m2/s GM thickness diffusivity
kappa_cnst = 4.2e-5; %m2/s Background vertical diffusivity
epsilon = 1e-3; % kg s-3 constant dissipation energy in the ocean per area
N0 = 8e-3; % 1/s used in N2 dependent diffusivity 
a0 = 1e-5; % m2/s constant diffusivity used in N2 dependent diffusivity
tau = 0.085; %Nm-2 = kgm-1s-2
Cgm = (1-exp(-tau/0.02)); % it can be 1 for simplivity
FN = 0.2*1e6; % Sv ; Northern meridional freshwater transport
FNd = 0;
FNFW = 0; %initial freshwater into nordic seas
FNFW_cap = 2e6; % cap of freshwater
FS = 0.2*1e6; %Sv ; Southern meridional freshwater transport
gammau = 1/(5*365*86400); % 5 years
gamman = 1/(5*365*86400); % 5 years
gammas = 1/(5*365*86400); % 5 years
TUrelax = 20;
TNrelax = 3;
TSrelax = 7;

% initial conditions 
SN=35; % Northern cell
SU=35; % Upper cell
SS=35; % Southern cell
SD=35; % Deep cell
TN = 5;
TS = 0;
TU = 12.5; %12.5;
TD = 3; 

H_pyc=500; % initial conditions
% initial conditions 
H_deep = Htopo - H_pyc;
if AADW
  SA = 35; % Antarctic DW cell
  TA = 0;
  H_AADW = 1000.0;
  H_deep = Htopo - H_pyc - H_AADW;
  gammaa = 1/(5*365*86400); % 25 years
  TArelax = 2;
  VA = LU*Bwidth*H_AADW;
end
VU = LU*Bwidth*H_pyc;
VN = LN*Bwidth*Htopo;
VS = LS*Bwidth*Htopo;
VD = LU*Bwidth*H_deep;
	
yearinsec = 360*86400; % 1 year in sec
deltat = 15*86400; %30 days to sec
time = 0;
iind = 1;
yearlength = 1200; 
%return

if 0
  load initial
end

%return
% time step
for ind=1:(yearlength*yearinsec/deltat)
    delta_rho = rho0*(beta_S*(SN-SU)-alpha_T*(TN-TU));
    delta_rho_SO = rho0*(beta_S*(SS-SU)-alpha_T*(TS-TU));
    delta_rho_D = rho0*(beta_S*(SN-SD)-alpha_T*(TN-TD));
    if AADW
      delta_rho_SO2 = rho0*(beta_S*(SS-SD)-alpha_T*(TS-TD));
      delta_rho_AADW = rho0*(beta_S*(SA-SD)-alpha_T*(TA-TD));
    end
    %if delta_rho_D<0
    %  display('MI')
    %end
    switch(diffusivity)
    case 'constant diff'
      kappa_v = kappa_cnst;  
      if AADW
        kappa_v2 = epsilon/(grav*delta_rho_AADW);
      end
    case 'N2 dependent'
      kappa_v = a0*(sqrt(grav*delta_rho/(rho0*H_pyc))/N0)^(-1);
    case 'constant energy'
      kappa_v = epsilon/(grav*delta_rho);
    otherwise
      fprintf('not specified diffusivity option')
      return
    end
    % phi_Moc
    phi_N = (C1*grav*delta_rho*H_pyc*H_pyc/(rho0*fbeta*LyN));
    if delta_rho<0
      display('MI'),iind
      phi_N = 0.0;
      return
    end
    % phi_Up
    phi_Up = (LU*Bwidth*kappa_v/H_pyc);
    % phi_Ek
    phi_Ek = (Bwidth*tau/(f0*rho0));
    % phi_GM
    phi_GM1 = Cgm*(Bwidth*kappa_GM*(delta_rho_SO/rho0)*H_pyc/Htopo); 
    % The new version
    %phi_GM1 = 0.01*Bwidth*(grav*delta_rho_SO/rho0)*H_pyc*H_pyc/(f0*LS);

    if AADW
      %phi_GM2 = 0.01*Bwidth*(grav*delta_rho_SO2/rho0)*H_deep*H_deep/(f0*LS);
      %phi_GM2 = 3*Cgm*(Bwidth*kappa_GM*(delta_rho_SO2/rho0)*H_deep/Htopo); 
      % phi_AADW (into the bottom layer AADW)
      %phi_AADW = (phi_GM1+phi_GM2)-phi_Ek;

      % phi Internal Tides Upwelling
      %phi_ITU = (LU*Bwidth*kappa_v2/H_AADW);
      phi_ITU = (LU*Bwidth*kappa_v2/H_deep);
      %phi_ITU = max(0,(LU*Bwidth*kappa_v2/H_deep));
      %
      phi_AADW = phi_ITU;
      phi_GM2 = phi_AADW+phi_Ek-phi_GM1;
      
      if(isinf(phi_ITU)==1); 
          display('phi_ITU')
          phi_ITU = 0.0; 
      end
    end
    %return

    if AADW
      % new AADW depth
      %H_AADW_new = max(0,H_AADW + (deltat/(Bwidth*LU))*(phi_GM1+phi_GM2-phi_Ek-phi_AADW));
      H_AADW_new = max(200.0,H_AADW + (deltat/(Bwidth*LU))*(phi_AADW-phi_ITU));
      if(H_AADW_new < 0)
        display('Mehmet AADW depleted2')
        H_AADW_new = 0.0;
        phi_ITU = 0.0;
        phi_AADW = 0.0;
      end
    end

    % Volume*Salt terms
    VUSU = VU*SU;
    VNSN = VN*SN;
    VSSS = VS*SS;
    VDSD = VD*SD;
    if AADW
      VASA = VA*SA;
    end
    if(tempevl)
      % Volume*Temp terms
      VUTU = VU*TU;
      VNTN = VN*TN;
      VSTS = VS*TS;
      VDTD = VD*TD;
      if AADW
        VATA = VA*TA;
      end
    end
    % New volume salt terms
    VUSUnew = VUSU + deltat*(phi_Up*SD+phi_Ek*SS-SU*(phi_N+phi_GM1)+S0*(FN+FS));
    VNSNnew = VNSN + deltat*(phi_N*(SU-SN)-S0*(FN+FNFW));
    if AADW
      VSSSnew = VSSS + deltat*(phi_Ek*(-SS)+phi_GM1*(SU)+phi_GM2*(SD)-phi_AADW*SS-S0*FS);
      VDSDnew = VDSD + deltat*(phi_N*SN-phi_GM2*SD-SD*(phi_Up)+phi_ITU*SA);
      VASAnew = VASA + deltat*(phi_AADW*SS-phi_ITU*SA);
    else
      VSSSnew = VSSS + deltat*(phi_Ek*(SD-SS)+phi_GM1*(SU-SS)-S0*FS);
      VDSDnew = VDSD + deltat*(phi_N*SN+phi_GM1*SS-SD*(phi_Up+phi_Ek));
    end
    if(tempevl)
      % New volume temp terms
      VUTUnew = VUTU + deltat*(phi_Up*TD+phi_Ek*TS-TU*(phi_N+phi_GM1)+gammau*VU*(TUrelax-TU));
      VNTNnew = VNTN + deltat*(phi_N*(TU-TN)+gamman*VN*(TNrelax-TN));
      if AADW
        VSTSnew = VSTS + deltat*(phi_Ek*(-TS)+phi_GM1*(TU)+phi_GM2*(TD)-phi_AADW*TS+gammas*VS*(TSrelax-TS));
        VDTDnew = VDTD + deltat*(phi_N*TN-phi_GM2*TD-TD*(phi_Up)+phi_ITU*TA);
        VATAnew = VATA + deltat*(phi_AADW*TS-phi_ITU*TA+gammaa*VA*(TArelax-TA));
      else
        VSTSnew = VSTS + deltat*(phi_Ek*(TD-TS)+phi_GM1*(TU-TS)+gammas*VS*(TSrelax-TS));
        VDTDnew = VDTD + deltat*(phi_N*TN+phi_GM1*TS-TD*(phi_Up+phi_Ek));
      end
    end

    % new pycnocline depth
    H_pyc_new = H_pyc + (deltat/(Bwidth*LU))*(phi_Up+phi_Ek-phi_GM1-phi_N);
    %update values
    H_pyc = H_pyc_new;
    if AADW
      % new H_deep
      H_AADW = H_AADW_new;
      H_deep = Htopo - H_pyc_new - H_AADW_new;
    else
      H_deep = Htopo - H_pyc_new;
    end
    VUSU = VUSUnew;
    VNSN = VNSNnew;    
    VSSS = VSSSnew;    
    VDSD = VDSDnew;
    VDSD = VDSDnew;
    if AADW
      VASA = VASAnew;
    end
    if(tempevl)
       VUTU = VUTUnew;
       VNTN = VNTNnew;
       VSTS = VSTSnew;
       VDTD = VDTDnew;
       if AADW
         VATA = VATAnew;
       end
    end
    VU = LU*Bwidth*H_pyc;
    VN = LN*Bwidth*Htopo;
    VS = LS*Bwidth*Htopo;
    VD = LU*Bwidth*H_deep;
    if AADW
      VA = LU*Bwidth*H_AADW;
    end
    SU=VUSU/VU;
    SN=VNSN/VN;
    SS=VSSS/VS;
    SD=VDSD/VD;
    if AADW
      SA = VASA/VA;
    end
    if(tempevl)
       TU = VUTU/VU;
       TN = VNTN/VN;
       TS = VSTS/VS;
       TD = VDTD/VD;
       if AADW
         TA = VATA/VA;
       end
    end
    % time series
    time=time+deltat;        
    if(mod(time,yearinsec)==0)
      if FWforce==true
        if(mod(iind,100)==0)
          FNd=FNFW;
        end
        FNFW=FNd+(1*FNFW_cap)*mod(iind,100)/99;
        if iind>500 & iind < 550
          FNFW = 1.0*FNFW_cap;
        else
          FNFW = 0.0*FNFW_cap;
        end
      end  
      Fwater(iind)=FNFW;
      HUtime(iind) = H_pyc;
      TrN(iind) = phi_N;
      TrW(iind) = phi_Ek;
      TrE(iind) = phi_GM1;
      TrU(iind) = phi_Up;
      SUtime(iind) = SU;
      SNtime(iind) = SN;
      SStime(iind) = SS;
      SDtime(iind) = SD;
      if AADW
        SAtime(iind) = SA;
        TrA(iind) = phi_AADW;
        TrI(iind) = phi_ITU;
      end
      drho(iind) = delta_rho;
      drhoD(iind) = delta_rho_D;
      drhoSO(iind) = delta_rho_SO;
      if(tempevl)
        TUtime(iind) = TU;
        TNtime(iind) = TN;
        TStime(iind) = TS;
        TDtime(iind) = TD;
        if AADW
          TAtime(iind) = TA;
        end
      end
      iind = iind+1;
    end
end   
if (tempevl)
  figure('Position', [100, 100, 1000, 400]);
	subplot(1,3,1)
	plot(SUtime,'k','linewidth',2)
	hold on
	plot(SDtime,'b','linewidth',2)
	plot(SNtime,'r','linewidth',2)
	plot(SStime,'g','linewidth',2)
	subplot(1,3,2)
  plot(TUtime,'k','linewidth',2)
  hold on
  plot(TDtime,'b','linewidth',2)
  plot(TNtime,'r','linewidth',2)
  plot(TStime,'g','linewidth',2)
  legend('Up','Deep','North','South','location','best')
  subplot(1,3,3)
  plot(TrU*m3s2Sv,'k','linewidth',2)
  hold on
  plot(TrW*m3s2Sv,'b','linewidth',2)
  plot(TrN*m3s2Sv,'r','linewidth',2)
  plot(TrE*m3s2Sv,'g','linewidth',2)
  legend('Up','Wind','North','GM','location','northeast')
else
	plot(SUtime,'k','linewidth',2)
	hold on
	plot(SDtime,'b','linewidth',2)
	plot(SNtime,'r','linewidth',2)
	plot(SStime,'g','linewidth',2)
end
%save('initial.mat','SN','SU','SD','SS','TN','TU','TD','TS','H_pyc','VU','VN','VS','VD')
