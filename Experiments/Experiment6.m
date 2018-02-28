
resolution=[];
spectralgap=[];
for i=5:21;
    [res,spg]=Experiment5(i);
    resolution=[resolution, res];
    spectralgap=[spectralgap, spg];
end
plot(resolution,spectralgap)