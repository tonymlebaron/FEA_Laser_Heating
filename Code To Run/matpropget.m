function [rhos,rhov,ks,cps,cpv,delh,T1,T2,refl,alphafilm,subrefl,alphasub,subcp,subk,subrho] = matpropget(film_name,subs_name,laser_type)
%materialprops.csv has stored material properties for mineral oil, lubricating
%oil, carbon steel, stainless steel, and copper. All SI units. This
%function retrieves those properties.

A=readtable('materialprops.csv','ReadRowNames',true); %Import CSV with material props

%% Lookup Material Thermal Properties
rhos=A{'Density',film_name};
rhov=A{'rhov',film_name};
ks=A{'Thermal Conductivity',film_name};
cps=A{'Specific Heat',film_name};
cpv=A{'cpv',film_name}
delh=A{'LatentHeat',film_name};
T1=A{'T1',film_name};
T2=A{'T2',film_name};

subrho=A{'Density',subs_name}; %kg/m^3
subcp=A{'Specific Heat',subs_name}; %J/kg*k
subk=A{'Thermal Conductivity',subs_name}; %W/m*k

if length(laser_type)==5 %If this is the NdYAG laser
    if laser_type=='NdYAG'
    refl=A{'Reflectance NdYAG',film_name};
    alphafilm=A{'Absorption Coeff NdYAG',film_name};
    subrefl=A{'Reflectance NdYAG',subs_name};
    alphasub=A{'Absorption Coeff NdYAG',subs_name};
    %disp(num2str(subabsc))
    end
elseif length(laser_type)==3
    if laser_type=='KrF'
    refl=A{'Reflectance KrF',film_name};
    alphafilm=A{'Absorption Coeff KrF',film_name};
    subrefl=A{'Reflectance KrF',subs_name};
    alphasub=A{'Absorption Coeff KrF',subs_name};    
    end
end




end
