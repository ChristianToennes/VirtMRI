# Javascript Kiss FFT #

This is a simple JavaScript project that implements 2D FFT and IFFT functions.
For a usage example, see the `benchmark()` function in `index.html`.

The implementation is based on the [Kiss FFT](http://kissfft.sf.net) library,
which has been compiled to JavaScript using
[emscripten](https://github.com/kripken/emscripten). This is the first time
I've used Kiss or emscripten, so expect bugs. Works for me™ in Firefox 44 and
Chrome 48.

## Try it out ##

    $ git clone https://github.com/frederikhermans/js-kiss-fft2
    $ python -m SimpleHTTPServer

Open [http://localhost:8000](http://localhost:8000) in your browser, click the
"Benchmark" button.

## Build it yourself ##

To build the library yourself, you need `emcc` on your $PATH. Just run:

    $ make

## TODO ##

* Simplify memory handling
* Compare performance to other JavaScript 2D FFT libraries
