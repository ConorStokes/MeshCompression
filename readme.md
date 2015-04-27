# Vertex Cache Optimised Mesh Compression

This is a proof of concept library expanding on previously presented index buffer compression mechanisms and adding compression for vertices as well.

## How does it work?

Basically, it expands on the index buffer compression using a parallelogram predictor and a kind of universal code (related to exponential golomb) along with an exponential moving average to give fast adaptive compression for vertex attributes. A more complete version of this is outlined in [this blog post](http://conorstokes.github.io/2015/04/28/adding-vertex-compression-to-index-buffer-compression/).
