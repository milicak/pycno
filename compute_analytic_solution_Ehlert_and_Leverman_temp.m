clear all
% Ehlert and Leverman 2014;
% Mechanism for potential strengthening of Atlantic overturning prior to collapse
% temperature evolution equations
% for git hub
tempevl=true;
% Table 1
C1 = 0.1; %nondimensional coeff for MOC transport
m3s2Sv = 1e-6; 
grav = 9.81; %m2/s

H = 4000; %in m ; Average depth of Atlantic Ocean basin
B = 1e7; % in m ; Average width of Atlantic Ocean
LN = 3.34*1e6; %in m ; Meridional extend of the northen box
LU = 8.9*1e6;  %in m ; Meridional extend of the tropical box
LS = 3.34*1e6; %in m ; Meridional extend of the southern box

rho0 = 1027; %kgm-3
S0 = 35; % psu ; Average salinity of Atlantic Ocean
LyN = 1.5*1e6; % Meridional extend of the northen outcropping
kappa_GM = 1000; %m2/s GM thickness diffusivity
kappa = 4e-5; %m2/s Background vertical diffusivity
alpha_T = 2.1e-4; % 1/C thermal coefficient
alphaT = alpha_T*rho0;
beta_S = 8e-4; % 1/psu haline coefficient
betaS = beta_S*rho0;

f0 = 7.5e-5; %1/s
fbeta = 2e-11; % 1/s
tau = 0.1; %Nm-2 = kgm-1s-2
Cgm = (1-exp(-tau/0.02)); % it can be 1 for simplivity
FN = 0.1*1e6; % Sv ; Northern meridional freshwater transport
FNd = 0;
FNFW = 0; %initial freshwater into nordic seas
FS = 0.1*1e6; %Sv ; Southern meridional freshwater transport
TN = 5;
TS = 7;
TU = 12.5; %12.5;
TD = 2; 
gammau = 1/(5*365*86400); % 5 years
gamman = 1/(5*365*86400); % 5 years
gammas = 1/(5*365*86400); % 5 years
TUrelax = 20;
TNrelax = 5;
TSrelax = 7;

SN=35; % initial conditions
SU=35; % initial conditions
SS=35; % initial conditions
SD=35; % initial conditions

H_pyc=500; % initial conditions
VU = LU*B*H_pyc;
VN = LN*B*H;
VS = LS*B*H;
VD = LU*B*(H-H_pyc);
	
yearinsec = 360*86400; % 1 year in sec
deltat = 15*86400; %30 days to sec
time = 0;
iind = 1;

if 1
  load initial
end

% time step
for ind=1:24000    
    delta_rho = rho0*(beta_S*(SN-SU)-alpha_T*(TN-TU));
    delta_rho_SO = rho0*(beta_S*(SS-SU)-alpha_T*(TS-TU));
    delta_rho_D = rho0*(beta_S*(SN-SD)-alpha_T*(TN-TD));
    % phi_Moc
    phi_N = (C1*grav*delta_rho*H_pyc*H_pyc/(rho0*fbeta*LyN));
    % phi_Up
    phi_Up = (LU*B*kappa/H_pyc);
    % phi_Ek
    phi_Ek = (B*tau/(f0*rho0));
    % phi_GM
    phi_GM = Cgm*(B*kappa_GM*(delta_rho_SO/rho0)*H_pyc/H); 
    % new pycnocline depth
    H_pyc_new = H_pyc + (deltat/(B*LU))*(phi_Up+phi_Ek-phi_GM-phi_N);
    % Volume*Salt terms
    VUSU=VU*SU;
    VNSN=VN*SN;
    VSSS=VS*SS;
    VDSD=VD*SD;
    if(tempevl)
      % Volume*Temp terms
      VUTU=VU*TU;
      VNTN=VN*TN;
      VSTS=VS*TS;
      VDTD=VD*TD;
    end
    % New volume salt terms
    VUSUnew = VUSU + deltat*(phi_Up*SD+phi_Ek*SS-SU*(phi_N+phi_GM)+S0*(FN+FS));
    VNSNnew = VNSN + deltat*(phi_N*(SU-SN)-S0*(FN+FNFW));
    VSSSnew = VSSS + deltat*(phi_Ek*(SD-SS)+phi_GM*(SU-SS)-S0*FS);
    VDSDnew = VDSD + deltat*(phi_N*SN+phi_GM*SS-SD*(phi_Up+phi_Ek));
    if(tempevl)
      % New volume temp terms
      VUTUnew = VUTU + deltat*(phi_Up*TD+phi_Ek*TS-TU*(phi_N+phi_GM)+gammau*VU*(TUrelax-TU));
      VNTNnew = VNTN + deltat*(phi_N*(TU-TN)+gamman*VN*(TNrelax-TN));
      VSTSnew = VSTS + deltat*(phi_Ek*(TD-TS)+phi_GM*(TU-TS)+gammas*VS*(TSrelax-TS));
      VDTDnew = VDTD + deltat*(phi_N*TN+phi_GM*TS-TD*(phi_Up+phi_Ek));
    end
    %update values
    H_pyc = H_pyc_new;
    VUSU = VUSUnew;
    VNSN = VNSNnew;    
    VSSS = VSSSnew;    
    VDSD = VDSDnew;
    VDSD = VDSDnew;
    if(tempevl)
       VUTU = VUTUnew;
       VNTN = VNTNnew;
       VSTS = VSTSnew;
       VDTD = VDTDnew;
    end
    VU = LU*B*H_pyc;
    VN = LN*B*H;
    VS = LS*B*H;
    VD = LU*B*(H-H_pyc);
    SU=VUSU/VU;
    SN=VNSN/VN;
    SS=VSSS/VS;
    SD=VDSD/VD;
    if(tempevl)
       TU = VUTU/VU;
       TN = VNTN/VN;
       TS = VSTS/VS;
       TD = VDTD/VD;
    end
    % time series
    time=time+deltat;        
    if(mod(time,yearinsec)==0)
      if(mod(iind,100)==0)
        FNd=FNFW;
      end
      FNFW=FNd+(1*1e6)*mod(iind,100)/99;
      if iind>200
        FNFW = 2e6;
      end
      Fwater(iind)=FNFW;
      HUtime(iind) = H_pyc;
      TrN(iind) = phi_N;
      TrW(iind) = phi_Ek;
      TrE(iind) = phi_GM;
      TrU(iind) = phi_Up;
      SUtime(iind) = SU;
      SNtime(iind) = SN;
      SStime(iind) = SS;
      SDtime(iind) = SD;      
      if(tempevl)
        TUtime(iind) = TU;
        TNtime(iind) = TN;
        TStime(iind) = TS;
        TDtime(iind) = TD;
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