# Read a list of (run, minField, maxField) from a text file and
# download all 6 photoField*fits files in directories
# workingDir/rootName/run
# Libraries
import os
import urllib.request
import numpy


def fetchSDSSphotoField(outdir, run, camcol):
    try:
        infile = "http://data.sdss3.org/sas/dr9/boss/photoObj/301/%d\
/photoField-%06d-%d.fits" % (
            run, run, camcol)
        # Download the fits file
        outfile = outdir + "photoField-%06d-%d.fits" % (run, camcol)
        print("retrieving:", outfile)
        urllib.request.urlretrieve(infile, outfile)
    except:
        print("some problem with run=%06d camCol=%d" % (run, camcol))
        print('infile=', infile)
        print('infile=', outfile)
        os.remove(outfile)
        pass
    # if the file doesn't exist, urllib still makes an (almost) empty file,
    # remove it...
    statinfo = os.stat(outfile)
    filesize = statinfo.st_size
    if filesize < 300:
        os.remove(outfile)
        return 0
    else:
        print("downloaded photoField-%06d-%d.fits" % (run, camcol))
        return 1


# Read in file with SDSS run data
objlist = numpy.loadtxt('data/Stripe82RunList.dat')
numberAll = 0
numberOK = 0
numberBad = 0
rootName = "SDSSdata"

for line in objlist:
    run = int(line[0])
    numberAll += 1
    outdir = "%s/%d/" % (rootName, run)
    try:
        os.stat(outdir)
    except:
        os.makedirs(outdir)

    for camcol in range(1, 7):
        print('retrieving file for camcol:', camcol)
        try:
            fetchSDSSphotoField(outdir, run, camcol)
            numberOK += 1
        except:
            print('problem with run', run)
            numberBad += 1

print("tried %d runs (times 6 columns)" % (numberAll))
print("success for %d files" % (numberOK))
print("failed on %d files" % (numberBad))
