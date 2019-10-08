function BW = createLFBMask(RGB)

hsvimage = rgb2hsv(RGB);
hImage = hsvimage(:,:,1);
vImage = hsvimage(:,:,3);

BW = vImage > 0.65 & vImage < 0.72 & hImage > 0.63 & hImage < 0.7;

end
