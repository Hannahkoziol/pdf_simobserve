#newest update from Dan from July 5

Phase_Center_GC = 'GALACTIC 0.0deg +0.0deg'

simobserve( project         = 'compact',
            skymodel        = 'flux10_B1snap_000.fits',
            indirection     = Phase_Center_GC,
            setpointings    = True,
            integration     = '10s',
            antennalist     = 'sma.cmzoom.compact.cfg',
            incenter        = '230GHz',
            inwidth         = '50MHz',
            inbright        = '',
            totaltime       = '6h',
            graphics        = 'both',
            verbose         = True,
            overwrite       = True,
            mapsize         = '2.1arcmin',
            incell          = '0.49arcsec' )

simobserve( project         = 'subcompact',
            skymodel        = 'flux10_B1snap_000.fits',
            indirection     = Phase_Center_GC,
            setpointings    = True,
            integration     = '10s',
            antennalist     = 'sma.cmzoom.subcompact.cfg',
            incenter        = '230GHz',
            inwidth         = '50MHz',
            inbright        = '',
            totaltime       = '6h',
            graphics        = 'both',
            verbose         = True,
            overwrite       = True,
            mapsize         = '2.1arcmin',
            incell          = '0.49arcsec' )

concat(vis       = ['compact/compact.sma.cmzoom.compact.noisy.ms',
                    'subcompact/subcompact.sma.cmzoom.subcompact.noisy.ms'],
       concatvis = 'complete.ms')


tclean( vis                      = 'complete.ms',
        imagename                = 'sim_temp.cont.tclean',
        specmode                 = 'mfs',
        deconvolver              = 'hogbom',
        imsize                   = 256,
        weighting                = 'briggs',
        robust                   = 0.5,
        niter                    = 100000,
        interactive              = False,
        gridder                  = 'mosaic',
        cell                     = '0.49arcsec',
        phasecenter              = Phase_Center_GC,
        threshold                = '1.5mJy',
        usemask                  = "auto-multithresh",
        sidelobethreshold        = 1.0,
        noisethreshold           = 3.0,
        lownoisethreshold        = 2.0,
        minbeamfrac              = 0.1,
        growiterations           = 75,
        negativethreshold        = 0.0,
        pbcor                    = True )

#exportfits(imagename = 'sim_temp.cont.tclean.image.pbcor',
#           fitsimage = 'sim_temp.cont.tclean.pbcor.fits',
#           dropdeg   = True)
