import sys
from spectral_weighting import SpectralWeighting 


def main():
    if len(sys.argv) >= 2:
        option = sys.argv[1]
    else:
        option = '-help'

    sigma = 0.0

    # Check smooth option
    if option == '-smooth' and len(sys.argv) >= 3:
        sigma = float(sys.argv[2])
        sys.argv = [sys.argv[0]] + sys.argv[3:]
        option = sys.argv[1]

    print(f"Option: {option}")

    # Help
    if len(sys.argv) <= 2 or option == '-help':
        print("Type one of the following:")
        print("python main.py -splat input.fits output.fits [bchan echan]")
        print("python main.py [-smooth sigma] -weighted_splat input.fits output.fits size [bchan echan]")
        print("python main.py -merge filelist.txt output.fits [bchan echan]")
        print("python main.py [-smooth sigma] -weighted_merge exponent filelist.txt output.fits size [bchan echan]")
        print("python main.py [-smooth sigma] -local_weights input.fits output.fits size")
        print("python main.py -gaussian_smooth sigma input.fits output.fits")
        return

    sp = SpectralWeighting()

    # Splat
    if option == '-splat' and len(sys.argv) >= 4:
        if len(sys.argv) == 4 and not sp.splat(sys.argv[2], sys.argv[3]):
            sys.exit(1)
        if len(sys.argv) == 6 and not sp.splat(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])):
            sys.exit(1)

    # Gaussian Smoothing
    elif option == '-gaussian_smooth' and len(sys.argv) == 5:
        if not sp.gaussian_smoothing(sys.argv[3], sys.argv[4], float(sys.argv[2])):
            sys.exit(1)

    # Local Noise
    elif option == '-local_noise' and len(sys.argv) == 5:
        if not sp.local_noise(sys.argv[2], sys.argv[3], int(sys.argv[4])):
            sys.exit(1)

    # Local Weights
    elif option == '-local_weights' and len(sys.argv) == 5:
        if not sp.local_weights(sys.argv[2], sys.argv[3], int(sys.argv[4]), sigma=sigma):
            sys.exit(1)

    # Weighted Splat
    elif option == '-weighted_splat':
        if len(sys.argv) == 5 and not sp.weighted_splat(sys.argv[2], sys.argv[3], int(sys.argv[4]), 0, 0, sigma):
            sys.exit(1)
        if len(sys.argv) == 7 and not sp.weighted_splat(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sigma):
            sys.exit(1)

    # Weighted Merge
    elif option == '-weighted_merge':
        print(len(sys.argv))
        print(sp.set_exponent(float(sys.argv[2])))
        if len(sys.argv) == 6 and (not sp.set_exponent(float(sys.argv[2])) or not sp.weighted_merge(sys.argv[3], sys.argv[4], int(sys.argv[5]), 0, 0, sigma)):
            sys.exit(1)
        if len(sys.argv) == 8 and (not sp.set_exponent(float(sys.argv[2])) or not sp.weighted_merge(sys.argv[3], sys.argv[4], int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sigma)):
            sys.exit(1)

    # Merge
    elif option == '-merge':
        if len(sys.argv) == 4 and not sp.weighted_merge(sys.argv[2], sys.argv[3], 0.0):
            sys.exit(1)
        if len(sys.argv) == 6 and not sp.weighted_merge(sys.argv[2], sys.argv[3], 0.0, int(sys.argv[4]), int(sys.argv[5])):
            sys.exit(1)


if __name__ == "__main__":
    main()
