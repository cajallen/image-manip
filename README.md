# Compiling
`make` should work.
Use build/main for executable, and sample/image.ext for images.

# Upscaling/Downscaling
Left is Gaussian, Middle is linear, Right is point.
These are for the files named downsampling.png, and upsampling.png.

The gaussian filtering methods I used make it look nicer than linear for easily arguably
both of them. 

# Difficulties
I have gamma-correct blurs, which is cool. I've had a very busy couple of weeks and didn't have
time to get started on this project early, so I didn't get the rotation finished.
For the same reason, my scale code is absolutely dreadful. I had a hard time understanding
both of the sampling methods, I think primarily because any change in coordinate system
is typically pretty confusing. I wish I had more time to figure these out and clean up my code,
not that it was the classes fault, just me personally.