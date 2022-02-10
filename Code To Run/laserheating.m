clear
close all
%parpool
%clc
%% Constants for material and laser
film_name='PMMA'; %Must match column heading in materialprops.csv EXACTLY
subs_name='CarbonSteel'; %Substrate name. Must match column heading EXACTLY
laser_type='NdYAG'; %Enter 'KrF' or 'NdYAG'
[rhos,rhov,ks,cps,cpv,delh,T1,T2,refl,alphafilm,subrefl,alphasub] = matpropget(film_name,subs_name,laser_type);

tamb=25; %Ambient Temperature (Deg. C)

subthick=1e-5; %Substrate thickness
filmthick=1e-5; %PMMA thickness

fluence=0.05*10000; %.10 J/cm^2. Enter quantity as W/m^2
tp=50e-9; %pulse width
I0=fluence/tp; %Peak irradiance W/m^2 
freq=10000; %10 kHz
numpulses=5; %Number of pulses you wish to apply to system

x=[-subthick:1e-7:-2e-8,-1e-8:1e-9:filmthick]; %Specify desired points for nodes as a 1xn array
%% Start a timer
tic
delh=delh/(T2-T1);
% Initialize Phase Change Pts array pcp as global variable to facilitate
% use in functions
global pcp
pcp=[]; %pcp stores x coordinate of points that have been vaporized (undergone phase change)
vaptime=[]; %vaptime stores the time at which points in pcp vaporize

%% Create the array of times where we want solutions.
% Want to record Temp at several times during pulse
% But less frequently when pulsing not active
realtime=[];
for i=1:numpulses
    realtime=[realtime,(i-1)/freq:tp/20:(i-1)/freq+4*tp,(i-1)/freq+4*tp+0.1/freq:0.1/freq:i/freq];
end

%Plot Irradiance vs Time
figure
for i=1:length(realtime)-1
    Ft(i)=FofT(1,tp,freq,i,realtime);
end
plot(realtime(1:end-1),I0*Ft,'LineWidth',2)
xlabel('Time (s)')
ylabel('Laser Irradiance (W/m^2)')
title('Laser Irradiance over Time')
drawnow
    

%% Initiate the solver


sols=25+zeros(1,length(x)); %Initializes array sols to store solutions at all of our timesteps. First row = IC
global stepsol; %Use global variables to "trick" pdepe into using temp. output from previous step
global x_vector; %as initial condition for new step. See function @theic for details
x_vector=x;
stepsol=sols; %Solution at end of initial "timestep"
for i=1:length(realtime)-1 %For all times in array t
    %The next line uses matlab's PDE solver pdepe. For more info type "help
    %pdepe"
    t=[0,(realtime(i+1)-realtime(i))/2,realtime(i+1)-realtime(i)];
    sol=pdepe(0,...
            @(x,t,u,DuDx)thepde(x,t,u,DuDx,T2,T1,cps,cpv,delh,ks,rhos,rhov,subrefl,alphafilm,refl,I0,tp,freq,i,realtime,pcp,filmthick,alphasub),...
            @theic,...
            @(xl,ul,xr,ur,t)thebcs(xl,ul,xr,ur,t,I0,refl,tp,freq,i,realtime,alphafilm,pcp,filmthick),...
            x,t); %This line solves the PDE
    stepsol=sol(end,:); %This is the final solution at time t. It is passed as the IC for next solver iteration.
    sols=[sols;stepsol]; %Record final solution

    for j=1:length(x) %Loop to check if any elements have vaporized
        if stepsol(j) >=T2-5 %If element j in stepsol is near vap. threshold T2
            if ismember(x(j),pcp)==0 %If element j is not already in pcp
                pcp=[pcp;x(j)]; %Then x(j) is added to the list of vaporized coordinates
                vaptime=[vaptime;realtime(i+1)]; %Record the time at which this point vaporized
            end
        end
    end
    if rem(i,floor(length(realtime)/20))==0 %This if statement tracks progress for user's convenience
        c=clock;
        disp(['Timestep ',num2str(i),' of ', num2str(length(realtime)),' complete at ', num2str(c(4)),':',num2str(c(5)),' hrs'])
        save('workspace.mat') %Save the workspace for backup just in case
    elseif i==1
        disp('First Timestep Complete!')
    end
end


toc

%% Plot Results
% 3D surface plot
figure
surf(x,realtime,sols,'EdgeColor','none','LineStyle','none') %Surface plot of solution over time
xlabel('x (m)')
ylabel('time (s)')
zlabel('Temp (C)')
print('-dpng','3DResults.png')

%Ablation vs. Time
figure
plot(vaptime,filmthick-pcp,'.','MarkerSize',15) %Ablation vs. Time
xlabel('time (s)')
ylabel('Ablation (m)')
title('Ablation vs. Time')
print('-dpng','Ablation.png')

%Plot temp. Distribution during pulses
for kk=0:1/freq:(numpulses-1)/freq
    tlist=[0 25e-9 50e-9 75e-9 100e-9 150e-9 200e-9]+kk;
    counter=0;
    figure
    hold on
    for ll=1:length(tlist)
       counter=counter+1;
       plot(x,sols(realtime==tlist(ll),:),'-','LineWidth',2,'DisplayName',strcat(num2str(tlist(ll)/1e-9),' ns')) %,colors{counter}
    end
    legend('show','Location','northwest')
    %xlim([-2e-6 2e-6])
    xlabel('z (m)')
    ylabel('Temperature (^{o}C)')
    title(['Temperature at Pulse ',num2str(kk*freq+1)])
    fname = ['TempDistribution',num2str(kk*freq+1),'.png']; % filename: JPEG file 
    print( '-dpng', fname ); % print figure: JPEG file 
end

%Compute and plot thermal energy density at x=0
ind0=find(abs(x)==min(abs(x))); %index where x=0 (or is closest to 0)
t0=sols(:,ind0); %Array of temperatures at x=0 over all timesteps
Cfun= @(u) (cps)*smhv(T2-u)+smhv(u - T1).*smhv(T2 - u).*delh+cpv.*smhv(u-T2); %Function for C
e=zeros(1,length(t0));
for i=1:length(t0)
    e(i)=rhos*integral(Cfun,25,t0(i)); %Integrate to find increase in thermal energy
end
figure
plot(realtime,e/(100^3),'LineWidth',2)
xlabel('Time (s)')
ylabel('\Delta e (J/cm^3)')
title('Increase in Thermal Energy Density')
print('-dpng','EnergyDensity.png')

%% Functions for execution of pdepe
function [c,f,s]=thepde(x,t,u,DuDx,T2,T1,cps,cpv,delh,ks,rhos,rhov,subrefl,alphafilm,refl,I0,tp,freq,i,realtime,pcp,filmthick,alphasub) %defines the governing PDE
    c=heaviside(x).*calcC1D(x,u,T2,T1,cps,cpv,delh).*calcrho1D(x,u,T2,T1,rhos,rhov)+heaviside(-x).*(7860*434); %c*rho for film and substrate
    f=calcK1D(x,u,T2,ks)*DuDx.*heaviside(x)+52*heaviside(-x)*DuDx; 
    s=heaviside(x).*subrefl*alphafilm*(1-refl)*I0.*FofT(t,tp,freq,i,realtime).*exp(-alphafilm*x).*exp(-alphafilm*min([pcp;filmthick]))+heaviside(-x).*(1-subrefl)*alphasub*(1-refl)*I0.*FofT(t,tp,freq,i,realtime).*exp(alphasub*x).*exp(-alphafilm*min([pcp;filmthick])); %Source heating term due to reflected light at interface
end

function u0=theic(x) %Use previous solution as IC for next step.
    %For details, see 
    %https://www.mathworks.com/matlabcentral/answers/299432-how-to-use-pdepe-with-arbitrary-initial-condition
    global stepsol
    global x_vector
    [R,C] = find(x_vector==x,1,'first');
    u0=stepsol(R,C); %Apply the previous solution as the new initial condition
end

function [pl,ql,pr,qr] = thebcs(xl,ul,xr,ur,t,I0,refl,tp,freq,i,realtime,alphafilm,pcp,filmthick) %Boundary Conditions
    pl=ul;
    ql=25; 
    pr=-I0*(1-refl).*FofT(t,tp,freq,i,realtime)*(1-exp(-alphafilm*min([pcp;filmthick]))); %Heat flux BC on right edge
    qr=1;
end
