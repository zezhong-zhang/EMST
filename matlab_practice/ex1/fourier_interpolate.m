function [arrayout] = fourier_interpolate(arrayin,newgrid,fftin,fftout)
%Fourier interpolation (and downsampling) routine. Interpolates by 
%zero padding and down samples by cropping in Fourier space. 
%fftin and fftout are optional arguments. Pass fftin=1
%if arrayin is already in Fourier space (assumed default), pass fftout=1 
%if you want the  output in Fourier space (not the default).

if(exist('fftin','var'))
    %Transforms into Fourier space or not as requested by user.
    if(fftin)
       arrayin_ = arrayin;
    else
       arrayin_ = fft2(arrayin);
    end
else
    %Default is to assume that the array is in real space
    arrayin_ = fft2(arrayin);
end

%Get size of oldgrid
oldgrid = size(arrayin);

%Initialize new array
arrayout = zeros(newgrid);

%Work out whether original array is odd or even
yodd = mod(oldgrid(1),2);
xodd = mod(oldgrid(2),2);

%length of list of elements from arrayin that will go into newarray
ylength = floor(min(newgrid(1),oldgrid(1))/2); %in y
xlength = floor(min(newgrid(2),oldgrid(2))/2); %in x

%Put elements of old array into new array
arrayout(mod([-ylength:ylength-1+yodd],newgrid(1))+1,...
         mod([-xlength:xlength-1+xodd],newgrid(2))+1)...
       = arrayin_(mod([-ylength:ylength-1+yodd],oldgrid(1))+1,...
                  mod([-xlength:xlength-1+xodd],oldgrid(2))+1);


if(exist('fftout','var'))
    %Inverse Fourier transform to real space if requested by user
    if(~fftout)
       arrayout = ifft2(arrayout);
    end
else
    %Default is to assume that the output is in real space
    arrayout = ifft2(arrayout);
end     

%Since a fft will multiply the magnitude by a factor (#pixels)
%and an ifft will divide by this number, we have to multiply by a factor
%(#pixelsout)/(#pixelsin) to maintain the correct scaling of the
%output

arrayout = arrayout* (prod(newgrid)/prod(oldgrid));

%If input is real then output should be too
if(isreal(arrayin))
    arrayout = abs(arrayout);
end

end