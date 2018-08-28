%%% Validation Toolbox V.01
%%% By: Ali Mehran and Amir AghaKouchak
%%%
%%% A simple and easy to use Validation Toolbox (MATLAB source code)that can be used for validation
%%% of gridded data including satellite observations, reanalysis data, and weather and climate model
%%% simulations. In addition to the commonly used categorical indices, the toolbox includes the 
%%% Volumetric Hit Index (VHI), Volumetric False Alarm Ration (VFAR),
%%% Volumetric Missed Index (VMI) and Volumetric Critical Success Index (VCSI).
%%% 
%%% Citations:
%%% AghaKouchak A., Mehran A., 2013, Extended Contingency Table: Performance Metrics for Satellite Observations and Model Simulations, Water Resources Research, doi:10.1002/wrcr.20498.
%%% Mehran A., AghaKouchak A., 2013, Capabilities of Satellite Precipitation datasets to simulated Heavy Precipitation Rates at Different Temporal Accumulations, Hydrological Processes, doi: 10.1002/hyp.9779.
%%% AghaKouchak A., Mehran A., Norouzi H., Behrangi A., 2012, Systematic and Random Error Components in Satellite Precipitation Data Sets, Geophysical Research Letters, L09406, 39, doi:10.1029/2012GL051592.
%%% AghaKouchak A., Behrangi A., Sorooshian S., Hsu K., Amitai E., 2011, Evaluation of satellite-retrieved extreme precipitation rates across the Central United States, Journal of Geophysical Research, 116, D02115, doi:10.1029/2010JD014741.
%%%
%%% Link:
%%% http://amir.eng.uci.edu/software.html
%%%
%%% Disclaimer:
%%% This program (hereafter, software) is provided 'as is' without warranty of any kind, either express or implied. The software could include technical or other mistakes, inaccuracies or typographical errors. The use of the software is done at your own discretion and risk and with agreement that you will be solely responsible for any damage and that the authors and their affiliate institutions accept no responsibility for errors or omissions in the software or documentation. In no event shall the authors or their affiliate institutions be liable to you or any third parties for any special, indirect or consequential damages of any kind, or any damages whatsoever.
%%%------------------------------------------------------------
clc
clear all
load sample_simulation
load sample_observation
sim = sample_simulation; obs = sample_observation; 
trr=0.0;  %Filters out values samller than trr (threshold), e.g., negative values
% Can also be used to filter out values below a certain quantile (e.g., 75th percentile).
minsamp = 50;   %Filters out pixels in which the total number of observations are less than minsamp to avoid unreliable statistics
sim(sim<trr) = 0.0; sim(sim==-999) = nan; %Replace -999 with the value that corresponds to "no data"
obs(obs<trr) = 0.0; obs(obs==-999) = nan; 
oii = size(sim,1); ojj = size(sim,2); okk = size(obs,3); template = nan(1,1,okk); result = [];     
bias = sum(sim,3)./sum(obs,3);
for ii=1:oii
    for j=1:ojj 
        temp = obs(ii,j,:); temp = temp(temp>trr); 
        if size(temp,3)<minsamp
            sim(ii,j,:)=template;
        end
        observed = obs(ii,j,:); simulated = sim(ii,j,:);
        [biasmap, hit_bias, NoHit, NoFalse, NoMiss, SumMiss, SumFalse, sumSimhit] = ValidationFunction(simulated, observed,trr);
        result=[result;  NoHit  NoFalse   NoMiss  NoFalse/(NoHit+NoFalse)   NoHit/(NoMiss+NoHit) SumMiss  SumFalse  sumSimhit biasmap];
    end
end
Lat=48.875:-.25:25.125;  Lon=-124.875:0.25:-67.125;

% Plot Bias
par9 = rot90(reshape(result(:,9),ojj,oii)); 
subplot(3,3,1)
pcolor(1:ojj,1:oii,par9); set(gca,'xticklabel','','yticklabel','') ;caxis([0 2]);
title('Bias','FontWeight','bold'); colorbar
% Plot Probability of Detection (POD)
par5 = rot90(reshape(result(:,5),ojj,oii)); 
subplot(3,3,2)
pcolor(1:ojj,1:oii,par5);set(gca,'xticklabel','','yticklabel','') ;caxis([0 1]);
title('POD','FontWeight','bold');colorbar
% plot Volumetric Hit Index (VHI)
par6 = rot90(reshape(result(:,6),ojj,oii)); par7(:,:) = rot90(reshape(result(:,7),ojj,oii));
par8 = rot90(reshape(result(:,8),ojj,oii)); vhit = par8./(par6+par8);vhit(vhit<0)=nan;
subplot(3,3,3) 
pcolor(1:ojj,1:oii,vhit);set(gca,'xticklabel','','yticklabel','') ;caxis([0 1]);
title('VHI','FontWeight','bold');colorbar
% Plot False Alarm Ratio (FAR)
par4 = rot90(reshape(result(:,4),ojj,oii)); 
subplot(3,3,4)
pcolor(1:ojj,1:oii,par4);
set(gca,'xticklabel','','yticklabel','') ;caxis([0 1]);
title('FAR','FontWeight','bold');colorbar
% Plot Volumetric False Alarm Ration (VFAR)
par7 = rot90(reshape(result(:,7),ojj,oii));
vfalse = par7./(par7+par8);vfalse(vfalse<0)=nan;
subplot(3,3,5) 
pcolor(1:ojj,1:oii,vfalse);set(gca,'xticklabel','','yticklabel','') ;caxis([0 1]);
title('VFAR','FontWeight','bold');colorbar
% Plot Categorical Miss
par5 = rot90(reshape(result(:,5),ojj,oii)); 
subplot(3,3,6)
pcolor(1:ojj,1:oii,1-par5);set(gca,'xticklabel','','yticklabel','') ;caxis([0 1]);
title('Categorical Miss','FontWeight','bold');colorbar
% plot Volumetric Miss Index (VMI)
vmissed = par6./(par8+par6);vmissed(vmissed<0)=nan;
subplot(3,3,7) 
pcolor(1:ojj,1:oii,vmissed);set(gca,'xticklabel','','yticklabel','') ;caxis([0 1]);
title('VMI','FontWeight','bold');colorbar
% % Plot Critical Success Index (CSI)
par1 = rot90(reshape(result(:,1),ojj,oii));
par2 = rot90(reshape(result(:,2),ojj,oii));
par3 = rot90(reshape(result(:,3),ojj,oii));
subplot(3,3,8) 
pcolor(1:ojj,1:oii,par1./(par1+par2+par3));set(gca,'xticklabel','','yticklabel','') ;caxis([0 1]);
title('CSI','FontWeight','bold');colorbar
% % Plot Volumetric Critical Success Index (VCSI)
subplot(3,3,9) 
pcolor(1:ojj,1:oii,par8./(par8+par7+par6));set(gca,'xticklabel','','yticklabel','') ;caxis([0 1]);
title('VCSI','FontWeight','bold');colorbar
