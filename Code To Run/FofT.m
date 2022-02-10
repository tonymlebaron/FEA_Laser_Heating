function [Ft] = FofT(timenow,tp,freq,i,realtime)
    timenow=(realtime(i)+realtime(i+1))/2; %Average time at this timestep
    timenow=timenow-(1/freq)*floor(timenow/(1/freq)); %Adjusts for pulse frequency
    if timenow < 2*tp
        Ft=(-1/tp)*abs(timenow-tp)+1;
    elseif timenow >=2*tp
        Ft=0;
    else
        Ft=0;
        disp('Error computing Ft, assuming 0')
    end
end

