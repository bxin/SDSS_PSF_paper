#!/usr/bin/env python

class sdssinst(object):

    def __init__(self):
        self.nBand = 5
        self.nCamcol = 6
        # show which band
        self.band = {}
        self.band[0] = "u"
        self.band[1] = "g"
        self.band[2] = "r"
        self.band[3] = "i"
        self.band[4] = "z"

