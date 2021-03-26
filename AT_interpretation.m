%% Simple AX interpretation code
% 2021/01/29 Calculates
clear all; close all; 
global H w w0 FT KF CHCl KW wHCl ST KS S k AT BT KB
% How to extract data will probably have to change, or be removed
% completely 
if pwd ~'/Users/May-Linn/Dropbox (Personal)/Matlab/AXdata';
    cd '/Users/May-Linn/Dropbox (Personal)/Matlab/AXdata'
end
FileName = '20210203 SW-1_sim.csv';
    [w0,S,emf0,t0,CNaOH,CHCl,emf1,t1,w1,emf2,t2,w2] = Extract_AXdata(FileName);
% Check sample type, I have set options for running "fake" seawater or
% chloride solutions vs real seawater, this might be relevant for people
% who use chloride/bicarbonate solutions to calibrate their acid, but I
% assume we'll need another way of setting this up. 
    if contains(FileName,'NaCl') == 1
        Sol = 1;
    elseif contains(FileName,'KCl') == 1
        Sol = 3;
    elseif contains(FileName,'Cl') == 0
        Sol = 2; 
    end
% This is work to set it up as a function, leaving it here for now in case that's
% relevant 
% other sample properties 
%     if isfield(Titrants,'CHCl')
%     CHCl = Titrants.CHCl; 
%     end
%     if isfield(Titrants,'CNaOH')
%     CNaOH = Titrants.CNaOH; 
%     end
%     if isfield(Titrants,'CT_NaOH')
%     CT_NaOH = Titrants.CT_NaOH; 
%     else
%     CT_NaOH = 0e-6;
%     end
%     CT_sample = 0e-6; 
    w02     = w0 + w1(end); 
%     if isfield(MinorAlks,'BTratio')
%     BTratio = MinorAlks.BTratio; 
%     else
    BTratio = 1;
%     end
%     if isfield(MinorAlks,'SiT')
%     SiT = MinorAlks.SiT; 
%     else
    SiT = 0e-6;
%     end
%     if isfield(MinorAlks,'PT')
%     PT = MinorAlks.PT; 
%     else
    PT = 0e-6;
%     end
    
    [k,KW,K1,K2,ST,FT,BT,KS,KF,KB,KSi,KP1,KP2,KP3] = EqConstants(S,t,Sol,BTratio); % This is a function that calculates equilibrium constants, I haven't included it but will make some edits and share it once I know what the code end result will look like 
    options = optimset('TolX',1e-16,'TolFun',1e-16,'Display','off','Algorithm','levenberg-marquardt');
    Idx     = find (emf > 0.205 & emf < 0.235); %pH 3.5–3 (for Metrohm ecotrode plus), but this part needs refinement to work with different electrodes, some iterative approach to find the right pH interval
%% Interpretation code, first Gran to estimate E0 then nonlinear to calculate AT and refine E0
    [p1,p2,weq] = Gran1(w0,w,emf,Idx);
    E01est      = (emf(Idx) - k.*log(CHCl.*(w(Idx) - weq)./(w0 + w(Idx))));
    H           = exp((emf - mean(E0est))./k); 
    ATest       = -p2/p1*CHCl/w0; 
    % Solver
    x0(1)       = 1;     %f
    x0(2)       = ATest; %AT
    [x]         = fsolve(@interp1, x0, options);
    r1          = 1000000*interp1(x); 
    s_of_sqrs   = r1*r1';
    E0          = mean(E0est) - k*log(x(1));
    H           = exp((emf - E0)./k); 
    AT          = x(2); 
    sprintf('E0 estimated from FWD data is %.6f, and AT %.2f umol/kg.',E0,AT*1e6)
%% Sub-routines 
    function r1 = interp1(x) 
    global w0 CHCl w H KW KF KS ST FT 
        f    = x(1); 
        AT   = x(2);
        Z    = 1 + ST/KS; 
        HSO4 = w0.*ST./(1+(Z*KS)./(f*H)); 
        HF   = w0.*FT./(1+KF./(f*H));
        r1   = w0*AT + HSO4 + HF - w.*CHCl + (w0 + w).*((f*H./Z) - (KW./(f*H)));
    end    
    function [p1,p2,weq] = Gran1(w0,w,emf,Idx)
    global k
        F1          = (w0+w(Idx)).*exp((emf(Idx))./k);
        ft          = fittype('poly1');
        fitresult   = fit(w(Idx)',F1',ft);
        p1          = fitresult.p1;
        p2          = fitresult.p2;
        weq         = -p2/p1;
    end 