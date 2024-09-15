There are two main methods:
the weighted splat and the weighted merge

The weighted splat method consists in:
1. Open an input fits file, and create an output fits file;
2. Compute the local weights spectral cube;
3. Loop over all spectral channels c:
   a) Read the local weights image
   b) Compute the weighted average
   c) Create the final splat image
4. Close input and output fits files;

The weighted merge method consists in:
1. Load fits files list;
2. Open input fits files, and create an output fits file;
3. For each input spectral cube compute the local weights spectral cube;
4. Loop over all spectral channels c:
   a) Read the local weights image
   b) Compute the weighted average
   c) Create the final merge image and add it to the output spectral cube
5. Close input and output fits files;

Functionalities: Select channel range for splat.
Particular Features: ”Mini Splat”: to eliminate broadband, but very weak RFI, we implemented a mini splat.
