import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

class SpectralWeighting:
    def __init__(self):
        self.float_nan = float('nan')
        self.exponent = 2

    def set_exponent(self, e):
        self.exponent = e
        return True

    def splat(self, infile, outfile, bchan=0, echan=0):
        try:
            # Open input file
            infits = fits.open(infile, mode='readonly')
        except Exception as e:
            print(f"Error opening input file: {e}")
            return False
        
        # Create output file
        try:
            outfits = fits.HDUList([fits.PrimaryHDU()])
            outfits.writeto(outfile, overwrite=True)
            outfits = fits.open(outfile, mode='update')
        except Exception as e:
            print(f"Error creating output file: {e}")
            return False

        # Clone header
        outfits[0].header = infits[0].header.copy()
        outfits[0].header['NAXIS'] = 2
        outfits[0].header['NAXIS1'] = infits[0].header['NAXIS1']
        outfits[0].header['NAXIS2'] = infits[0].header['NAXIS2']
        outfits[0].data = np.zeros((infits[0].header['NAXIS2'], infits[0].header['NAXIS1']), dtype=np.float32)

        # Determine channel range
        if echan == 0:
            echan = infits[0].header['NAXIS3']
        else:
            echan -= 1
        if bchan > 0:
            bchan -= 1

        print(f"Splat cube from channel {bchan + 1} to channel {echan + 1}")

        # Initialize arrays
        nx = infits[0].header['NAXIS1']
        ny = infits[0].header['NAXIS2']
        sum_image = np.zeros((ny, nx), dtype=np.float32)
        wsum_image = np.zeros((ny, nx), dtype=np.float32)
        splat_image = np.full((ny, nx), self.float_nan, dtype=np.float32)

        # Main loop
        for k in range(bchan, echan):
            print(f"Working on channel: {k + 1} of {echan + 1}", end='\r')

            # Read the channel image
            image = infits[0].data[k, :, :]

            # Process the image
            mask = np.isfinite(image)
            sum_image[mask] += image[mask]
            wsum_image[mask] += 1.0

        # Calculate the splat image
        valid = wsum_image > 0
        splat_image[valid] = sum_image[valid] / wsum_image[valid]

        # Write the final splat image
        outfits[0].data = splat_image
        outfits.flush()

        # Close files
        infits.close()
        outfits.close()

        return True


    def gaussian_smoothing(self, infile, outfile, sigma):
        # Compute the smoothing kernel dimensions
        m = 1 + 2 * int(3.0 * sigma)  # It's always odd
        n = 1 + 2 * int(3.0 * sigma)

        # Open input file
        try:
            infits = fits.open(infile, mode='readonly')
        except Exception as e:
            print(f"Error opening input file: {e}")
            return False

        # Create output file
        try:
            outfits = fits.HDUList([fits.PrimaryHDU()])
            outfits.writeto(outfile, overwrite=True)
            outfits = fits.open(outfile, mode='update')
        except Exception as e:
            print(f"Error creating output file: {e}")
            return False

        # Clone header
        outfits[0].header = infits[0].header.copy()
        outfits[0].header['NAXIS1'] = infits[0].header['NAXIS1']
        outfits[0].header['NAXIS2'] = infits[0].header['NAXIS2']
        outfits[0].header['NAXIS3'] = infits[0].header['NAXIS3']

        nx = infits[0].header['NAXIS1']
        ny = infits[0].header['NAXIS2']
        nz = infits[0].header['NAXIS3']

        for k in range(nz):
            print(f"Working on channel: {k + 1} of {nz}", end='\r')
            image = infits[0].data[k, :, :]

            # Apply Gaussian smoothing 
            smooth_image = gaussian_filter(image, sigma=sigma, mode='constant', cval=self.float_nan)

            # Handle NaN regions (preserve NaNs in the original data)
            mask = np.isfinite(image)
            smooth_image[~mask] = self.float_nan

            if k == 0:
                outfits[0].data = np.zeros((nz, ny, nx), dtype=np.float32)

            outfits[0].data[k, :, :] = smooth_image

        outfits.flush()
        infits.close()
        outfits.close()

        return True


    def get_plane_sigma_image(self, image, nx, ny, sigma_image, m):
        for j in range(ny):
            for i in range(nx):
                values = []

                # Extract nearby pixel values
                for ii in range(m):
                    for jj in range(m):
                        xpix = i + ii - int(m / 2.0 + 0.5)
                        ypix = j + jj - int(m / 2.0 + 0.5)

                        # If inside image
                        if 0 <= xpix < nx and 0 <= ypix < ny:
                            value = image[xpix + nx * ypix]
                            if np.isfinite(value):
                                values.append(value)

                # Calculate statistic of nearby pixel values
                if values:
                    sigma_image[i + nx * j] = np.std(values)  # Standard Deviation
                else:
                    sigma_image[i + nx * j] = self.float_nan

        return True


    def local_weights(self, infile, outfile, size, bchan=0, echan=0, sigma=0):
        # Smooth input image if sigma > 0
        smooth_infile = infile
        if sigma > 0:
            smooth_infile = "smoothed_" + infile
            if not self.gaussian_smoothing(infile, smooth_infile, sigma):
                return False

        # Open input file
        try:
            infits = fits.open(smooth_infile, mode='readonly')
        except Exception as e:
            print(f"Error opening input file: {e}")
            return False

        print(f"Creating file {outfile}")

        # Create output file
        try:
            outfits = fits.HDUList([fits.PrimaryHDU()])
            outfits.writeto(outfile, overwrite=True)
            outfits = fits.open(outfile, mode='update')
        except Exception as e:
            print(f"Error creating output file: {e}")
            return False

        # Clone header
        outfits[0].header = infits[0].header.copy()
        outfits[0].header['NAXIS1'] = infits[0].header['NAXIS1']
        outfits[0].header['NAXIS2'] = infits[0].header['NAXIS2']
        outfits[0].header['NAXIS3'] = infits[0].header['NAXIS3']
        outfits[0].data = np.zeros((infits[0].header['NAXIS3'], infits[0].header['NAXIS2'], infits[0].header['NAXIS1']), dtype=np.float32)

        nx = infits[0].header['NAXIS1']
        ny = infits[0].header['NAXIS2']
        nz = infits[0].header['NAXIS3']

        m = 1 + 2 * size  # It's always odd

        # Select channel range for weights' cube
        if echan == 0:
            echan = nz
        else:
            echan -= 1
        if bchan > 0:
            bchan -= 1

        print(f"Weights' cube from channel {bchan + 1} to channel {echan + 1}")
    
        mini_splat_image = np.full(nx*ny, self.float_nan, dtype=np.float32)
        sum_image = np.zeros(ny*nx, dtype=np.float32)
        wsum_image = np.zeros(ny*nx, dtype=np.float32)

        for k in range(bchan, echan):
            print(f"Working on channel: {k + 1} of {echan + 1}", end='\r')
            sum_image.fill(0.0)
            wsum_image.fill(0.0)

            # Compute a mini splat of the nearby channels to highlight faint broad-band RFI
            k1 = max(bchan, k - m)
            k2 = min(echan, k + m)

            for kk in range(k1, k2):
                image = infits[0].data[kk, :, :].flatten()
                for i in range(nx*ny):
                    if np.isfinite(image[i]):
                        sum_image[i] += image[i]
                        wsum_image[i] += 1.0

                for i in range(nx*ny): 
                    if np.isfinite(sum_image[i]) and wsum_image[i] > 0:
                        mini_splat_image[i]=sum_image[i]/wsum_image[i]
                
            sigma_image = np.zeros_like(mini_splat_image)
            self.get_plane_sigma_image(mini_splat_image, nx, ny, sigma_image, m)

            # Compute weights as the inverse of the standard deviation raised to the exponent
            with np.errstate(divide='ignore', invalid='ignore'):  
                weights = 1.0 / np.power(sigma_image, self.exponent)
                weights[~np.isfinite(weights)] = self.float_nan
            
            for i in range(nx):
                for j in range(ny):
                    outfits[0].data[k, j, i] = weights[i+j*nx]
        
        outfits.flush()

        infits.close()
        outfits.close()

        return True


    def weighted_splat(self, infile, outfile, size, bchan=0, echan=0, sigma=0):
        # Open input file
        try:
            infits = fits.open(infile, mode='readonly')
        except Exception as e:
            print(f"Error opening input file: {e}")
            return False

        # Create output file
        try:
            outfits = fits.HDUList([fits.PrimaryHDU()])
            outfits.writeto(outfile, overwrite=True)
            outfits = fits.open(outfile, mode='update')
        except Exception as e:
            print(f"Error creating output file: {e}")
            return False

        # Clone header
        outfits[0].header = infits[0].header.copy()
        outfits[0].header['NAXIS'] = 2
        outfits[0].header['NAXIS1'] = infits[0].header['NAXIS1']
        outfits[0].header['NAXIS2'] = infits[0].header['NAXIS2']
        outfits[0].data = np.zeros((infits[0].header['NAXIS2'], infits[0].header['NAXIS1']), dtype=np.float32)

        # Calculate weights cube
        print("Calculating local weights")
        weights_file = "weights_" + infile
        if not self.local_weights(infile, weights_file, size, bchan, echan, sigma):
            return False

        # Open weights file
        try:
            winfits = fits.open(weights_file, mode='readonly')
        except Exception as e:
            print(f"Error opening weights file: {e}")
            return False

        # Select channel range for splat
        if echan == 0:
            echan = infits[0].header['NAXIS3']
        else:
            echan -= 1
        if bchan > 0:
            bchan -= 1

        print(f"Splat cube from channel {bchan + 1} to channel {echan + 1}")

        nx = infits[0].header['NAXIS1']
        ny = infits[0].header['NAXIS2']
        sum_image = np.zeros(ny*nx, dtype=np.float32)
        wsum_image = np.zeros(ny*nx, dtype=np.float32)
        splat_image = np.full(ny*nx, self.float_nan, dtype=np.float32)

        for k in range(bchan, echan):
            print(f"Working on channel: {k + 1} of {echan + 1}", end='\r')

            # Read image and weights for the current channel
            image = infits[0].data[k, :, :].flatten()
            wimage = winfits[0].data[k, :, :].flatten()

            # Update sum and weighted sum
            for i in range(nx*ny):
                if np.isfinite(image[i]):
                    sum_image[i] += image[i] * wimage[i]
                    wsum_image[i] += wimage[i]

            for i in range(nx*ny): 
                if np.isfinite(sum_image[i]) and wsum_image[i] > 0:
                    splat_image[i] = sum_image[i] / wsum_image[i]

        # Write the final splat image to the output file
        for i in range(nx):
            for j in range(ny):
                outfits[0].data[j, i] = splat_image[i+j*nx]
        outfits.flush()

        # Close files
        infits.close()
        winfits.close()
        outfits.close()

        return True


    def weighted_merge(self, filelist, outfile, size, bchan=0, echan=0, sigma=0):
        # Read the list of input files
        print("Starting weighted merge")
        try:
            with open(filelist, 'r') as f:
                files = [line.strip() for line in f if line.strip() and not line.startswith('#')]
        except Exception as e:
            print(f"Cannot open filelist: {e}")
            return False

        if not files:
            print("No input files found in filelist")
            return False

        print(f"Processing {len(files)} files")

        # Produce weights cubes if size > 0
        # If size=0, no weight is applied
        if size > 0:
            for filename in files:
                if not self.local_weights(filename, "weights_" + filename, size, bchan, echan, sigma):
                    return False

        # Open input FITS files and weights files if needed
        vinfits = []
        vwinfits = []

        for filename in files:
            try:
                vinfits.append(fits.open(filename, mode='readonly'))
                if size > 0:
                    vwinfits.append(fits.open("weights_" + filename, mode='readonly'))
                else:
                    vwinfits.append(None)
            except Exception as e:
                print(f"Error opening input or weight files: {e}")
                return False

        # Create output FITS file
        try:
            outfits = fits.HDUList([fits.PrimaryHDU()])
            outfits.writeto(outfile, overwrite=True)
            outfits = fits.open(outfile, mode='update')
        except Exception as e:
            print(f"Error creating output file: {e}")
            return False

        # Clone header
        outfits[0].header = vinfits[0][0].header.copy()
        outfits[0].header['NAXIS1'] = vinfits[0][0].header['NAXIS1']
        outfits[0].header['NAXIS2'] = vinfits[0][0].header['NAXIS2']
        outfits[0].header['NAXIS3'] = vinfits[0][0].header['NAXIS3']
        infits = vinfits[0]
        outfits[0].data = np.zeros((infits[0].header['NAXIS3'], infits[0].header['NAXIS2'], infits[0].header['NAXIS1']), dtype=np.float32)

        # Determine the channel range
        if echan == 0:
            echan = vinfits[0][0].header['NAXIS3']
        else:
            echan -= 1
        if bchan > 0:
            bchan -= 1

        print(f"Merge cube from channel {bchan + 1} to channel {echan + 1}")

        nx = vinfits[0][0].header['NAXIS1']
        ny = vinfits[0][0].header['NAXIS2']

        for k in range(bchan, echan):
            print(f"Working on channel: {k + 1} of {echan + 1}", end='\r')

            sum_image = np.zeros(ny*nx, dtype=np.float32)
            wsum_image = np.zeros(ny*nx, dtype=np.float32)
            merge_image = np.full(ny*nx, self.float_nan, dtype=np.float32)

            # Loop over files
            for i, vinfit in enumerate(vinfits):
                image = vinfit[0].data[k, :, :].flatten()

                if size > 0:
                    wimage = vwinfits[i][0].data[k, :, :].flatten()
                else:
                    wimage = np.ones_like(image)
                
                for i in range(nx*ny):
                    if np.isfinite(image[i]):
                        sum_image[i] += image[i] * wimage[i]
                        wsum_image[i] += wimage[i]
              
                for i in range(nx*ny): 
                    if np.isfinite(sum_image[i]) and wsum_image[i] > 0:
                        merge_image[i] = sum_image[i] / wsum_image[i]

            for i in range(nx):
                for j in range(ny):
                    outfits[0].data[k, j, i] = merge_image[i+j*nx]

        outfits.flush()

        # Close files
        for vinfit in vinfits:
            vinfit.close()
        for wfit in vwinfits:
            if wfit:
                wfit.close()

        outfits.close()

        return True