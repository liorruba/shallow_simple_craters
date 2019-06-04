function [T] = sphericalCraterTemperature(solarConstant, depthToDiameter, albedo, emissivity, solarElevationAngle)
%
% This function calculates the shadow temperature a spherical crater would
% have according to the model built by Ingersoll et al., 1992.
%
% Input:
% solarConstant -           the solar constant at the modeled crater.
% depthToDiameter -         the crater depth to diameter. This can either
%                           be a matrix or a vector.
% albedo -                  the surface visible albedo.
% emissivity -              the surface IR emissivity.
% solarElevationAngle -     the solar elevation angle measured in degress. 
%                           This can either be a matric or a vector. 
%
% Output:
% T -                       The equilibrium shadow temperature of the
%                           surface. If either depthToDiameter or
%                           solarElevationAngle is a vector, T is also a
%                           vector. If both are vectors, T is a matrix.
%
% Example:
% T = sphericalCraterTemperature(1367, 0.2, 0.136, 0.95, 1.6)
%
% Written by Lior Rubanenko, University of California Los Angeles,
% Last edited May 2016.
%

def.solarConstant = 1367;
def.depthToDiameter = 0.2;
def.albedo = 0.136;
def.emissivity = 0.95;
def.solarElevationAngle = 1.6;

if nargin == 0
    solarConstant = def.solarConstant;
    depthToDiameter = def.depthToDiameter;
    albedo = def.albedo;
    emissivity = def.emissivity; 
    solarElevationAngle = def.solarElevationAngle;

elseif nargin < 2
    depthToDiameter = def.depthToDiameter;
    albedo = def.albedo;
    emissivity = def.emissivity; 
    solarElevationAngle = def.solarElevationAngle;
    
elseif nargin < 3
    albedo = def.albedo;
    emissivity = def.emissivity; 
    solarElevationAngle = def.solarElevationAngle;
    
elseif nargin < 4
    emissivity = def.emissivity; 
    solarElevationAngle = def.solarElevationAngle;
    
elseif nargin < 5
    solarElevationAngle = def.solarElevationAngle;

elseif nargin > 5
    error('Too many input arguments.');
end   

if size(depthToDiameter) ~= size(solarElevationAngle)
    if length(depthToDiameter) > 1
        if length(solarElevationAngle) > 1
            [solarElevationAngle, depthToDiameter] = meshgrid(solarElevationAngle, depthToDiameter);
        end
    end
end

solarElevationAngle(solarElevationAngle<0)=0;

f=1./(1 + (1./depthToDiameter).^2 ./ 4);
solarFlux = (solarConstant .* sind(solarElevationAngle) .* f * (1 - albedo) ./ (1 - albedo .* f)) .* (1 + albedo .* (1 - f) ./ emissivity);

T = nthroot(solarFlux ./ 5.67e-8, 4);
