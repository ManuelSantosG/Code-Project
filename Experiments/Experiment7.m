
resolution=[];
spectralgap=[];
for i=5:21;
    depth=i;
    Experiment4
    resolution=[resolution, length(P)];
    spectralgap=[spectralgap, 1-abs(lambda(2,2))];
end
plot(resolution, spectralgap)