for Teff in {2300..5000..100} ; do
    fname="lte0${Teff}-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"

    curl ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/$fname > phoenix/$fname
done

curl ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS//WAVE_PHOENIX-ACES-AGSS-COND-2011.fits > phoenix/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits