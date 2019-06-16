
install.packages("readMzXmlData")
library("readMzXmlData")
#A list of spectra and metadata.
#[[1]]spectrum$mass: A vector of calculated mass.
#[[1]]spectrum$intensity: A vector of intensity values.
#[[1]]metaData: A list of metaData depending on read spectrum.

library("MALDIquant")
s<- createMassSpectrum(mass= 1:10, intensity= 11:20, metaData=list(name="Example spectrum"))
##get intensity
intensity(s)
##get mass
mass(s)
## get metadata
metaData(s)
##trim spectrum
install.packages("MALDIquant")
library("MALDIquant")
data("fiedler2009subset", package="MALDIquant")
## choose only spectrum 1
s1 <- fiedler2009subset[[1]]
## preprocessing
## sqrt transform (for variance stabilization)
s2 <- transformIntensity(s1, method="sqrt")
## 21 point Savitzky-Golay-Filter for smoothing spectra
## (maybe you have to adjust the halfWindowSize;
## you could use a simple moving average instead)
## see ?smoothIntensity
s3 <- smoothIntensity(s2, method="SavitzkyGolay", halfWindowSize=10)
## remove baseline
## (maybe you have to adjust iterations to your spectra; high resolution
## spectra need a much lower iteration number (halfWindowSize, for some other
## baseline estimation algorithms)
## see ?removeBaseline, ?estimateBaseline
s4 <- removeBaseline(s3, method="SNIP", iterations=100)
## run peak detection
## (maybe you need to adjust halfWindowSize [decreasing it for high resolution
## spectra] and SNR [a higher value increase the True-Positive-Rate but decrease
## sensitivity])
## see ?detectPeaks, ?estimateNoise
p <- detectPeaks(s4, method="MAD", halfWindowSize=20, SNR=2)
## produce some plots
par(mfrow=c(2,3))

xlim <- range(mass(s1)) # use same xlim on all plots for better comparison
plot(s1, main="1: raw", sub="", xlim=xlim)
plot(s2, main="2: variance stabilization", sub="", xlim=xlim)
plot(s3, main="3: smoothing", sub="", xlim=xlim)
plot(s4, main="4: baseline correction", sub="", xlim=xlim)
plot(s4, main="5: peak detection", sub="", xlim=xlim)
points(p)
## label top 20 peaks
top20 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:20]
top20 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:20]
labelPeaks(p, index=top20, underline=TRUE)
plot(p, main="6: peak plot", sub="", xlim=xlim)
labelPeaks(p, index=top20, underline=TRUE)
par(mfrow=c(1,1))



#In-silico Peptide Fragmentation
#The fragment ions of a peptide can be computed following the rules proposed in [4]. Beside
the b and y ions the FUN argument of fragmentIon defines which ions are computed. the
default ions beeing computed are defined in the function defaultIon. The are no limits for
defining other forms of fragment ions for ETD (c and z ions) CID (b and y ions)
defaultIon

library(protViz)
function (b, y)
{
  Hydrogen <- 1.007825
  Oxygen <- 15.994915
  Nitrogen <- 14.003074
  c <- b + (Nitrogen + (3 * Hydrogen))
  z <- y - (Nitrogen + (3 * Hydrogen))
  return(cbind(b, y, c, z))}

peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')
pim<-parentIonMass(peptides)
fi<-fragmentIon(peptides)
par(mfrow=c(3,1))
Hydrogen<-1.007825
(fi.HTLNQIDSVK.1<-fragmentIon('HTLNQIDSVK'))[[1]]
(fi.HTLNQIDSVK.2<-(fi.HTLNQIDSVK.1[[1]] + Hydrogen) / 2)
#Peptide Sequence - Fragment Ion Matching
peptideSequence<-'HTLNQIDSVK'
spec<-list(scans=1138,title="178: (rt=22.3807) [20080816_23_fetuin_160.RAW]",rtinseconds=1342.8402, charge=2, mZ=c(195.139940, 221.211970, 239.251780, 290.221750,316.300770, 333.300050, 352.258420, 448.384360, 466.348830,496.207570, 509.565910, 538.458310, 547.253380, 556.173940,560.358050, 569.122080, 594.435500, 689.536940, 707.624790,803.509240, 804.528220, 822.528020, 891.631250, 909.544400,916.631600, 973.702160, 990.594520, 999.430580, 1008.583600,1017.692500, 1027.605900,1017.692500, 1027.605900,1593, 531.8, 520.4, 976.4, 410.5, 2756, 2279, 5819, 2.679e+05,1267, 1542, 979.2, 9577, 3283, 9441, 1520, 1310, 1.8e+04,587.5, 2685, 671.7, 3734, 8266, 3309))
fi<-fragmentIon(peptideSequence)
n<-nchar(peptideSequence)
by.mZ<-c(fi[[1]]$b, fi[[1]]$y)
by.mZ
by.label<-c(paste("b",1:n,sep=''), paste("y",n:1,sep=''))
by.label
idx<-findNN(by.mZ, spec$mZ)

# For the assignment of a canditate
peptide an in-silico fragment ion spectra fi is computed. The function findNN determines
for each fragment ion the closesed peak in the MS2. If the difference between the in-silico
mass and the measured mass is inside the 'accuracy' mass window of the mass spec device
the in-silico fragment ion is considered as potential hit.

i<-1:5
findNN_(3.5, i)
i<-1:6
findNN_(3.5, i)
