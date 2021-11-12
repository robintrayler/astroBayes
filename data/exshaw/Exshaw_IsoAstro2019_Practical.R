# Notes on running astrochron at WWCC IsoAstro Workshop 2019
# Friday, June 7, 2019

# Today I will use R and astrochron and it will change my life

# use this to clear the working environment if you'd like

# I have to load up astrochron every time you restart R
library(astrochron)

# I want to read my comma separated value (CSV) data file (Site_1262_a.csv)
# that contains color data (red/green) from site 1262 (Walvis Ridge)
# We are going to use a function from astrochron called 'read()'
jura = read.csv(file = './data/exshaw/juraciso.csv')# or you can use <- instead of =
# rund = read()

# What is inside jura? let's see the raw data
# jura

# this takes the file you are going to specify and put it into the variable object "jura"
# it opens a wizard to open up your data file and does a bunch of stuff on this data
# including sorting, removing empty entries, and averaging
# ?read

# what if I want to read in tab-delimined data?
# you can now see all of the arguments that you could include in the command rather than
# accepting the defaults

# what is jura? it is something called a data.frame
# this is how R stores data for your analysis, which can store all different types of data
# together, e.g. text, boolean, numeric, complex variables
# class(jura)

# let's plot the data
plot(jura)  # this plot doesn't look so good with points

# let's plot the data with a line
plot(jura,type="l")
# plot(rund,type="l")

# let's plot the data with lines and points
plot(jura,type="b")

# let's plot the data with lines and points and labels
plot(jura,type="b",xlab="Depth (mcd)",ylab="a",cex=0.75)


# # let's get the data between the ELMO and jura and put it into a new data.frame 'jura'
# # we'll use the function 'iso' to isolate the data with a graphical interface
# jura = iso(jura)
#
# # if you want to learn more about 'iso', type ?iso
# ?iso
#
# # for this demo, let's make sure we all have the same data segment
# jura = iso(dat=jura,xmin=117.62,xmax=139.22)

# what is the original sampling?
# the function 'strats' will allow you to evaluate this...
strats(jura)
# strats(rund)

# now I have the data I want to analyze
# let's put the data on an even sample grid using interpolation
# most fourier methods require data with a constant sampling interval
# we'll use the function 'linterp' to do this
jura5cm=linterp(jura)
# rund5cm <- linterp(rund)

# by default your data set will be interpolated to the median sampling interval
# if you want to change this use the option 'dt'
jura5cm=linterp(jura, dt=0.05)
# rund5cm <- linterp(rund, dt=0.05)

# more info on linterp
# ?linterp

# let's calculate a periodogram, the simplest type of spectral analysis
periodogram(jura5cm)
# periodogram(rund5cm)

# let's look at the options for the function periodogram
# ?periodogram

# change scaling of periodogram to linear rather than logarithmic
periodogram(jura5cm,pl=2)

# Back to the geology, what are plausible sedimentation rates?
# Based on biostratigraphy (1.2 Ma duration through 2.5 m black shale),
# the nominal sedimentation rate is 2 m/Ma or 0.002 mm/a.

# we need to think about smaple resolution and what's reasonable given
# plausible sedimentation rate...
# if sed rate is 2 mm/a, then precession is around 40 mm = 4 cm, or 25 cycles/m
# and obliquity is around 80 mm or 8 cm, or 12.5 cycles/m
# and short eccentricity is around 200 mm or 20 cm, or 5 cycles/m
# and long eccentricity is around 800 mm or 80 cm, or 1.2 cycles/m
# if sed rate is > 1 cm/ka then the peak in the periodogram moves to the left
# assume minimum sedimentation rate ~0.5 mm/kyr,
# therefore we should examine out to 100 cycles/m
periodogram(jura5cm)

# now let's see why interpreting the periodogram is a BAD IDEA…
# see examples in handout

# let's investigate our site data using the multitaper method
mtm(jura5cm)

# let's plot this on a linear power scale
mtm(jura5cm, pl=2)

# let's also remove the linear trend too
mtm(jura5cm,pl=2,detrend=T)

# no matter what we do we can't seem to see a 'significant' eccentricity signal
# over the red noise, assuming a sed rate of 1.18 cm/ka
# let's reinterpolate

jura5cm = linterp(jura,dt=0.05)

mtm(jura5cm,pl=2,detrend=T)
mtm(rund5cm,pl=2,detrend=T)

# let's ues a better approach to noise modeling: lowspec
lowspec(jura5cm,pl=2,detrend=T)

# now we'll test the Milankovitch hypothesis!
# output the frequencies from lowspec

freq = lowspec(jura5cm,pl=2,detrend=T,padfac=10,sigID=T,siglevel=0.95,
               output=2)

### Conduct ASM testing on Site 1262 a* data
# set Rayleigh frequency in cycles/m
rayleigh=0.04640371
# set Nyquist frequency in cycles/m
nyquist=10

# this target was derived using the following:
#e=getLaskar("la11")
#e=iso(e,xmin=53000,xmax=56000)
#eha(e,win=3000,sigID=T,fmin=0,fmax=0.02,pad=30000)
#p=etp(tmin=53000,tmax=56000,dt=1,eWt=0,oWt=0,pWt=1,esinw=T)
#eha(p,win=3000,sigID=T,fmin=0.02,fmax=0.07,pad=30000)

target=c(1/405.4,1/124.5,1/95.2,1/23,1/21.8,1/18.6)

asm(freq,target=target,rayleigh=rayleigh,nyquist=nyquist,sedmin=0.5,sedmax=2,
    numsed=100,iter=100000)

# let's do evolutive harmonic analysis

eha(jura5cm,tbw=2,win=1,step=0.05,pad=1000)
eha(rund5cm,tbw=2,win=1,step=0.05,pad=1000)

# learn more about astrochron? includes an example demos and help!
?astrochron

timeOpt(jura5cm)





