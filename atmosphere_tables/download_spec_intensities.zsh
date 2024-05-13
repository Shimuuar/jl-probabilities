for Teff in {2300..4900..100} ; do
    fname="lte0${Teff}-4.00-0.0.PHOENIX-ACES-AGSS-COND-SPECINT-2011.fits"

    curl ftp://phoenix.astro.physik.uni-goettingen.de/SpecIntFITS/PHOENIX-ACES-AGSS-COND-SPECINT-2011/Z-0.0/$fname > phoenix_spec_intensities/$fname
done
